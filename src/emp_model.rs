#![allow(dead_code)]
use rand::distr::{Distribution, weighted::WeightedIndex};
use rand::prelude::*;

pub struct EmpHistMode {
    pub prob: Vec<f64>,
    pub edge: Vec<u32>,
    pub median: u32,
}

impl EmpHistMode {
    pub fn new(prob: Vec<f64>, edge: Vec<u32>, median: u32) -> Self {
        assert_eq!(
            edge.len(),
            prob.len() + 1,
            "Lengths of prob and edge are inconsistent!"
        );
        Self { prob, edge, median }
    }

    pub fn from_lengths(lengths: &mut [u32], bin_width: u32) -> Self {
        let max_length = *lengths.iter().max().unwrap_or(&1000);
        let num_bins = (max_length / bin_width) + 1;
        let mut bin_counts = vec![0usize; num_bins as usize];

        for length in lengths.iter() {
            let bin_idx = (*length / bin_width) as usize;
            bin_counts[bin_idx] += 1;
        }

        let total: usize = bin_counts.iter().sum();
        let prob: Vec<f64> = bin_counts
            .iter()
            .map(|&count| count as f64 / total as f64)
            .collect();

        let edge: Vec<u32> = (0..=num_bins).map(|i| i * bin_width).collect();

        let median = lengths.select_nth_unstable(lengths.len() / 2).1;

        let (edge, prob) = Self::truncate_outliers(edge, prob);
        Self::new(prob, edge, *median)
    }

    pub fn sample(&self, mut rng: &mut StdRng, target_median: Option<u32>) -> u32 {
        let distribution = WeightedIndex::new(&self.prob).unwrap();
        let bin_idx = distribution.sample(&mut rng);
        let start = self.edge[bin_idx];
        let end = self.edge[bin_idx + 1];
        if let Some(target_median) = target_median {
            let value = rng.random_range(start..end);
            (value as f64 * target_median as f64 / self.median as f64) as u32
        } else {
            rng.random_range(start..end)
        }
    }

    pub fn samples(
        &self,
        num_samples: usize,
        target_median: Option<u32>,
        mut seed: u64,
    ) -> Vec<u32> {
        if seed == 0 {
            seed = rand::rng().random();
        }
        let mut rng = StdRng::seed_from_u64(seed);
        let distribution = WeightedIndex::new(&self.prob).unwrap();

        let mut samples: Vec<u32> = (0..num_samples)
            .map(|_| {
                let bin_idx = distribution.sample(&mut rng);
                let start = self.edge[bin_idx];
                let end = self.edge[bin_idx + 1];
                rng.random_range(start..end)
            })
            .collect();

        if let Some(target_median) = target_median {
            Self::rescale_to_median(&mut samples, target_median as f64);
        }
        samples
    }

    fn truncate_outliers(edge: Vec<u32>, prob: Vec<f64>) -> (Vec<u32>, Vec<f64>) {
        let peak_index = prob
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);

        let mut start = peak_index;
        let mut end = peak_index;

        let mut zero_span = 0;
        while start > 0 && zero_span < 10 {
            if prob[start] > 0.0 {
                zero_span = 0;
            } else {
                zero_span += 1;
            }
            start -= 1;
        }

        zero_span = 0;
        while end + 1 < prob.len() && zero_span < 10 {
            if prob[end] > 0.0 {
                zero_span = 0;
            } else {
                zero_span += 1;
            }
            end += 1;
        }

        let mut trimmed_probs = prob[start..=end].to_vec();
        trimmed_probs[0] += prob[..start].iter().sum::<f64>();
        trimmed_probs[end - start] += prob[end + 1..].iter().sum::<f64>();

        let mut trimmed_edges = edge[start..=end + 1].to_vec();
        trimmed_edges[0] = edge[0];
        trimmed_edges[end - start + 1] = edge[edge.len() - 1];

        (trimmed_edges, trimmed_probs)
    }

    fn rescale_to_median(samples: &mut [u32], target_median: f64) {
        let mid = samples.len() / 2;
        let current_median = *samples.select_nth_unstable(mid).1 as f64;
        let scale = target_median / current_median;

        for value in samples.iter_mut() {
            *value = ((*value as f64 * scale).round().max(1.0)) as u32;
        }
    }
}
