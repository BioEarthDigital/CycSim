use crossbeam_channel::bounded;
use edlib_rs::edlibrs::*;
use kseq::parse_path;
use log::info;
use mimalloc::MiMalloc;
use nohash_hasher::IntMap;
use parking_lot::RwLock;
use rand::{
    distr::{
        Bernoulli, Uniform,
        weighted::{Weight, WeightedIndex},
    },
    prelude::*,
};
use rayon::prelude::*;
use rust_htslib::bam::{
    self, HeaderView, Read as Read_, Record, Writer, ext::BamRecordExtensions, record::*,
};
use rustc_hash::{FxBuildHasher, FxHashMap as HashMap};
use std::{
    sync::{Arc, Mutex},
    thread,
};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

mod mlib;
use mlib::{FileIO, resource_str};

mod option;
use option::Option as Opt;

mod mode_io;
use mode_io::{dump_file, load_file};

mod emp_model;
use emp_model::EmpHistMode;

mod sim;
use sim::{sim_aligned_seq, sim_unaligned_seq};

mod kmer;
use kmer::{iter2kmer, seq2kmer};

mod qual;

mod reads;
use reads::{SimReads, create_bam_header};

const COUNTER_BITS: usize = 10;
const UNALIGNED_ERROR_RATE: f64 = 0.4;
const TRANSIT_WIN: usize = 10;
const DELETION_CIGAR: usize = 0b01; //same as edlib cigar
const INSERTION_CIGAR: usize = 0b10; //same as edlib cigar
const MISMATCH_CIGAR: usize = 0b11; //same as edlib cigar
const KMER_ERROR_BIN: f32 = 0.01;
const LEN_ERROR_BIN: f32 = 1000.;

fn read_seqs_with_tid(path: &str, header: &HeaderView) -> Vec<Vec<u8>> {
    let mut seqs = vec![Vec::new(); header.target_count() as usize];

    let mut records = parse_path(path).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        if let Some(tid) = header.tid(record.head().as_bytes()) {
            seqs[tid as usize] = record.seq().to_ascii_uppercase().as_bytes().to_owned();
        }
    }
    seqs
}

fn read_seqs_with_names(path: &str) -> Vec<(String, Vec<u8>)> {
    let mut seqs = Vec::new();

    let mut records = parse_path(path).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        seqs.push((
            record.head().to_owned(),
            record.seq().to_ascii_uppercase().as_bytes().to_owned(),
        ));
    }
    seqs
}

fn get_read_aligned_range(cigar: &CigarStringView, mut end: u32) -> (usize, usize) {
    let mut start = 0;

    if let Some(Cigar::SoftClip(len)) = cigar.iter().next() {
        start = *len;
    } else if let Some(Cigar::HardClip(len)) = cigar.iter().next() {
        start = *len;
    }

    if let Some(Cigar::SoftClip(len)) = cigar.iter().next_back() {
        end -= *len;
    } else if let Some(Cigar::HardClip(len)) = cigar.iter().next_back() {
        end -= *len;
    }

    (start as usize, end as usize)
}

fn parse_cigar_span(cigar: &str) -> (usize, usize, usize, usize) {
    let (mut left_clip, mut right_clip, mut read_span, mut ref_span) = (0, 0, 0, 0);
    let mut num = String::new();
    let mut is_first = true;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            let len = num
                .parse::<usize>()
                .unwrap_or_else(|_| panic!("Invalid number before {ch}"));

            match ch {
                'S' | 'H' => {
                    if is_first {
                        left_clip = len;
                    } else {
                        right_clip = len;
                    }
                }
                'M' | 'X' | '=' => {
                    read_span += len;
                    ref_span += len;
                }
                'I' => {
                    read_span += len;
                }
                'D' | 'N' => {
                    ref_span += len;
                }
                _ => {
                    panic!("Unrecognized CIGAR operation: '{ch}'");
                }
            }
            is_first = false;
            num.clear();
        }
    }

    assert!(
        num.is_empty(),
        "CIGAR string ended unexpectedly with trailing digits"
    );

    (left_clip, right_clip, read_span, ref_span)
}

type ErrorPatternType = u32;

#[derive(Clone, Eq, Hash, PartialEq)]
struct ErrorPattern {
    ref_skip: usize,  // skipped sequence length in reference
    seq_sym: Vec<u8>, // substituted/inserted bases
}

impl ErrorPattern {
    const NUMBER_MAX: usize = u16::MAX as usize;

    fn new(seq_sym: &[u8], ref_skip: usize) -> Self {
        Self {
            seq_sym: seq_sym.to_owned(),
            ref_skip,
        }
    }

    fn error_type(&self) -> usize {
        if self.is_match() {
            0
        } else if self.is_deletion() {
            DELETION_CIGAR
        } else if self.is_insertion() {
            INSERTION_CIGAR
        } else {
            MISMATCH_CIGAR
        }
    }

    fn is_match(&self) -> bool {
        self.ref_skip == 0 && self.seq_sym.is_empty()
    }

    fn is_mismatch(&self) -> bool {
        self.ref_skip > 0 && !self.seq_sym.is_empty()
    }

    fn is_deletion(&self) -> bool {
        self.ref_skip > self.seq_sym.len()
    }

    fn is_insertion(&self) -> bool {
        self.ref_skip < self.seq_sym.len()
    }

    /// Parses an encoded error number.
    /// Returns: (type, length)
    /// - type: 0 => mismatch, 1 => DELETION, 2 => INSERTION
    fn parse_error_num(num: usize) -> (usize, usize) {
        (num & 0b11, num >> 2)
    }

    /// Encodes the error into two u16-like values.
    /// Returns: (mismatch, indel)
    fn to_error_num(&self) -> (usize, usize) {
        let max_val = (1 << 14) - 1;
        let seq_len = self.seq_sym.len().min(max_val);
        let ref_len = self.ref_skip.min(max_val);

        if seq_len == ref_len {
            ((ref_len << 2) | MISMATCH_CIGAR, 0)
        } else if seq_len > ref_len {
            (
                (ref_len << 2) | MISMATCH_CIGAR,
                ((seq_len - ref_len) << 2) | INSERTION_CIGAR,
            )
        } else {
            (
                (seq_len << 2) | MISMATCH_CIGAR,
                ((ref_len - seq_len) << 2) | DELETION_CIGAR,
            )
        }
    }
}

struct ErrorPatternInterner {
    map: RwLock<HashMap<ErrorPattern, ErrorPatternType>>,
    vec: RwLock<Vec<ErrorPattern>>,
}

impl ErrorPatternInterner {
    fn new() -> Self {
        Self {
            map: RwLock::new(HashMap::with_capacity_and_hasher(65536, FxBuildHasher)),
            vec: RwLock::new(Vec::with_capacity(65536)),
        }
    }

    fn get_or_insert(&self, pattern: ErrorPattern) -> ErrorPatternType {
        // Try read lock first
        if let Some(&sym) = self.map.read().get(&pattern) {
            return sym;
        }

        // Upgrade to write lock
        let mut map = self.map.write();
        if let Some(&sym) = map.get(&pattern) {
            return sym; // double-checked locking
        }

        let mut vec = self.vec.write();
        let sym = ErrorPatternType::try_from(vec.len())
            .expect("exceeds the maximum range of ErrorPatternType");
        map.insert(pattern.clone(), sym);
        vec.push(pattern);
        sym
    }

    // fn resolve(&self, sym: ErrorPatternType) -> ErrorPattern {
    //     self.vec.read()[sym as usize].clone()
    // }
}

fn extend_op_to_match(op: &[u8]) -> (usize, usize, usize) {
    let (mut i, mut ref_offset, mut read_offset) = (0, 0, 0);
    while i < op.len() && op[i] != 0 {
        if op[i] == 1 {
            // deletion
            ref_offset += 1;
        } else if op[i] == 2 {
            //insertion
            read_offset += 1;
        } else {
            //mismatch
            ref_offset += 1;
            read_offset += 1;
        }
        i += 1;
    }

    (i, ref_offset, read_offset)
}

fn trim_alignment_to_long_match(
    ops: &[u8],
    min_len: usize,
) -> Option<((usize, usize), (usize, usize), (usize, usize))> {
    let (mut ref_left_offset, mut ref_right_offset) = (0, 0);
    let (mut read_left_offset, mut read_right_offset) = (0, 0);
    let (mut ops_left_offset, mut ops_right_offset) = (0, 0);

    let mut match_count = 0;
    for op in ops.iter().copied() {
        if op == 0 {
            match_count += 1;
            ref_left_offset += 1;
            read_left_offset += 1;
            ops_left_offset += 1;
            if match_count >= min_len {
                ops_left_offset -= min_len;
                ref_left_offset -= min_len;
                read_left_offset -= min_len;
                break;
            }
        } else {
            match_count = 0;
            ops_left_offset += 1;
            match op {
                1 => ref_left_offset += 1,
                2 => read_left_offset += 1,
                3 => {
                    ref_left_offset += 1;
                    read_left_offset += 1;
                }
                _ => unreachable!("Not a valid op {op}!"),
            }
        }
    }

    if ops_left_offset > ops.len() - min_len {
        return None;
    }

    match_count = 0;
    for op in ops.iter().rev().copied() {
        if op == 0 {
            match_count += 1;
            ref_right_offset += 1;
            read_right_offset += 1;
            ops_right_offset += 1;
            if match_count >= min_len {
                ops_right_offset -= min_len;
                ref_right_offset -= min_len;
                read_right_offset -= min_len;
                break;
            }
        } else {
            match_count = 0;
            ops_right_offset += 1;
            match op {
                1 => ref_right_offset += 1,
                2 => read_right_offset += 1,
                3 => {
                    ref_right_offset += 1;
                    read_right_offset += 1;
                }
                _ => unreachable!("Not a valid op {op}!"),
            }
        }
    }

    if ops_right_offset + ops_left_offset >= ops.len() {
        return None;
    }

    Some((
        (ref_left_offset, ref_right_offset),
        (read_left_offset, read_right_offset),
        (ops_left_offset, ops_right_offset),
    ))
}

fn uniq_count(slice: &[u8]) -> impl Iterator<Item = (usize, usize)> {
    let mut iter = slice.iter().cloned().peekable();
    std::iter::from_fn(move || {
        let first = iter.next()?;
        let mut count = 1;
        while let Some(peeked) = iter.peek() {
            if peeked == &first {
                iter.next();
                count += 1;
            } else {
                break;
            }
        }
        Some((count, first as usize))
    })
}

fn analyze_cigar_ops(
    ops: &[u8],
    error_transit: &mut [[[u32; 5]; 5]],
    match_len: &mut [[[u32; 62]; 5]],
) {
    let mut last_op = 4; //start cigar by default
    let mut last_error = 4; //start cigar by default
    let mut last_len = 0;

    for (count, op) in uniq_count(ops) {
        let bin1 = (last_len / TRANSIT_WIN).min(20);
        let bin2 = (count / TRANSIT_WIN).min(60);
        if op == 0 {
            match_len[bin1][last_op][bin2 + 1] += 1;
        } else {
            if last_op == 0 {
                error_transit[bin1 + 1][last_error][op] += 1;
            } else {
                match_len[bin1][last_op][0] += 1;
                error_transit[0][last_error][op] += 1;
            }
            last_error = op;
        }
        last_op = op;
        last_len = count;
    }
}

fn sum_3d<const M: usize, const K: usize>(total: &mut [[[u32; K]; M]], part: &[[[u32; K]; M]]) {
    for (i, bin) in total.iter_mut().enumerate() {
        for (j, op) in bin.iter_mut().enumerate() {
            for (k, v) in op.iter_mut().enumerate() {
                *v = v.checked_add(part[i][j][k]).expect("u32 overflow");
            }
        }
    }
}

fn build_distribution<const N: usize>(
    data: &[[[u32; N]; 5]],
    skip_indexes: Vec<usize>,
) -> Vec<Vec<WeightedIndex<u32>>> {
    data.iter()
        .map(|bin| {
            bin.iter()
                .map(|weights| {
                    WeightedIndex::new(weights).unwrap_or_else(|_| {
                        WeightedIndex::new(
                            std::iter::repeat(1)
                                .enumerate()
                                .map(|(p, x)| if skip_indexes.contains(&p) { 0 } else { x })
                                .take(N),
                        )
                        .unwrap()
                    })
                })
                .collect()
        })
        .collect()
}

struct KmerDist {
    keys: Vec<u32>,
    pos_dist: WeightedIndex<u32>,
    pos_error_dist: WeightedIndex<u32>,
    pos_has: ErrorPresenceFlags,
    neg_dist: WeightedIndex<u32>,
    neg_error_dist: WeightedIndex<u32>,
    neg_has: ErrorPresenceFlags,
    error_rate: f32,
    match_index: usize,
}

#[derive(Default)]
struct ErrorPresenceFlags {
    mismatch: bool,
    deletion: bool,
    insertion: bool,
}

impl KmerDist {
    const ALL: usize = 0;

    fn new(map: &IntMap<u32, [u32; 2]>, errors: &[ErrorPattern]) -> Self {
        let mut keys = Vec::with_capacity(map.len());
        let mut pos_weights = Vec::with_capacity(map.len());
        let mut neg_weights = Vec::with_capacity(map.len());

        let mut pos_flags = ErrorPresenceFlags::default();
        let mut neg_flags = ErrorPresenceFlags::default();

        let (mut match_count, mut total_count) = (0, 0);
        let mut match_index = map.len();

        for (p, (&k, &[pos, neg])) in map.iter().enumerate() {
            keys.push(k);
            pos_weights.push(pos);
            neg_weights.push(neg);

            let e = &errors[k as usize];

            if e.is_match() {
                match_index = p;
                match_count += pos + neg;
            }
            total_count += (pos + neg) as usize * e.ref_skip.max(e.seq_sym.len()).max(1);

            if pos > 0 {
                pos_flags.mismatch |= e.is_mismatch();
                pos_flags.deletion |= e.is_deletion();
                pos_flags.insertion |= e.is_insertion();
            }
            if neg > 0 {
                neg_flags.mismatch |= e.is_mismatch();
                neg_flags.deletion |= e.is_deletion();
                neg_flags.insertion |= e.is_insertion();
            }
        }

        KmerDist {
            keys,
            pos_dist: WeightedIndex::new(&pos_weights).unwrap_or_else(|_| Self::default_dist()),
            pos_error_dist: WeightedIndex::new(Self::pseudo_zero_chain(&pos_weights, match_index))
                .unwrap_or_else(|_| Self::default_dist()),
            pos_has: pos_flags,
            neg_dist: WeightedIndex::new(&neg_weights).unwrap_or_else(|_| Self::default_dist()),
            neg_error_dist: WeightedIndex::new(Self::pseudo_zero_chain(&neg_weights, match_index))
                .unwrap_or_else(|_| Self::default_dist()),
            neg_has: neg_flags,
            error_rate: 1. - ((match_count as f64) / (total_count as f64)) as f32,
            match_index,
        }
    }

    fn sample<'a>(
        &self,
        is_pos: bool,
        error_type: usize,
        rng: &mut StdRng,
        errors: &'a [ErrorPattern],
        temperature: f64,
    ) -> Option<(&'a ErrorPattern, f64)> {
        // probability of correct base
        let (dist, err_dist, flags) = if is_pos {
            (&self.pos_dist, &self.pos_error_dist, &self.pos_has)
        } else {
            (&self.neg_dist, &self.neg_error_dist, &self.neg_has)
        };

        match error_type {
            DELETION_CIGAR => {
                if flags.deletion {
                    self.sample_filtered(err_dist, rng, errors, |e| e.is_deletion())
                        .map(|x| (x, 0.))
                } else {
                    None
                }
            }
            INSERTION_CIGAR => {
                if flags.insertion {
                    self.sample_filtered(err_dist, rng, errors, |e| e.is_insertion())
                        .map(|x| (x, 0.))
                } else {
                    None
                }
            }
            MISMATCH_CIGAR => {
                if flags.mismatch {
                    self.sample_filtered(err_dist, rng, errors, |e| e.is_mismatch())
                        .map(|x| (x, 0.))
                } else {
                    None
                }
            }
            Self::ALL => {
                let (idx, pro) = if temperature != 1. {
                    let dist = Self::update_index_weight_with_temperature(
                        dist,
                        temperature,
                        self.match_index,
                    );
                    (
                        dist.sample(rng),
                        dist.weight(self.match_index).unwrap_or(0) as f64
                            / dist.total_weight() as f64,
                    )
                } else {
                    (
                        dist.sample(rng),
                        dist.weight(self.match_index).unwrap_or(0) as f64
                            / dist.total_weight() as f64,
                    )
                };
                self.keys.get(idx).map(|x| (&errors[*x as usize], pro))
            }
            _ => panic!("Unknown error_type value: {error_type}"),
        }
    }

    fn sample_filtered<'a, F, W>(
        &self,
        dist: &WeightedIndex<W>,
        rng: &mut dyn RngCore,
        errors: &'a [ErrorPattern],
        pred: F,
    ) -> Option<&'a ErrorPattern>
    where
        W: Weight + rand::distr::uniform::SampleUniform + std::cmp::PartialOrd,
        F: Fn(&ErrorPattern) -> bool,
    {
        loop {
            let idx = dist.sample(rng);
            let err = &errors[self.keys[idx] as usize];
            if pred(err) {
                return Some(err);
            }
        }
    }

    fn pseudo_zero_chain(weights: &[u32], zero_index: usize) -> impl Iterator<Item = &u32> {
        weights
            .iter()
            .enumerate()
            .map(move |(p, x)| if p == zero_index { &0 } else { x })
    }

    fn default_dist() -> WeightedIndex<u32> {
        WeightedIndex::new(std::iter::once(&1)).unwrap()
    }

    fn update_index_weight_with_temperature(
        dist: &WeightedIndex<u32>,
        temperature: f64,
        index: usize,
    ) -> WeightedIndex<u32> {
        let mut new_dist = dist.to_owned();
        if let Some(weight) = dist.weight(index) {
            let new_weight = ((weight as f64 * temperature) as u32).clamp(1, u32::MAX / 2);
            let _ = new_dist.update_weights(&[(index, &new_weight)]);
            new_dist
        } else {
            new_dist
        }
    }
}

fn build_kmer_features<const N: usize>(
    features: &[Vec<IntMap<u32, [u32; 2]>>; N],
    errors: &[ErrorPattern],
) -> [Vec<KmerDist>; N] {
    let mut results: [Vec<_>; N] = std::array::from_fn(|_| Vec::new());
    results.par_iter_mut().enumerate().for_each(|(i, slot)| {
        let vecs = &features[i];
        *slot = vecs.iter().map(|map| KmerDist::new(map, errors)).collect();
    });

    results
}

fn build_kmer_identity_distribution(
    errors_rate: &[(f32, f32)],
    use_len: bool,
    global_error_rate: bool,
) -> Vec<Option<(Vec<u32>, WeightedIndex<u32>)>> {
    let rounded_key = |x: f32| (x * 10000.0).round() as u32;

    let kmer_max_error = errors_rate
        .iter()
        .map(|(err, _)| *err)
        .fold(0.0_f32, f32::max);

    let bin = if use_len {
        LEN_ERROR_BIN
    } else {
        KMER_ERROR_BIN
    };

    let bin_num = if global_error_rate { 1 } else { ((kmer_max_error / bin).ceil() as usize) + 1 };

    let mut binned = vec![HashMap::default(); bin_num];

    for &(kmer_error, identity) in errors_rate {
        let bin_idx = if global_error_rate { 0 } else { (kmer_error / bin).floor() as usize };
        let rounded_identity = rounded_key(identity);
        *binned[bin_idx].entry(rounded_identity).or_insert(0) += 1;
    }

    let mut identities_buffer = Vec::new();
    let mut dist_buffer = Vec::new();
    let mut identity_samplers = Vec::with_capacity(bin_num);

    for map in binned.into_iter() {
        // println!("{:?}", map);
        if map.len() < 3 {
            identity_samplers.push(None);
        } else {
            identities_buffer.clear();
            dist_buffer.clear();
            for (k, v) in map.into_iter() {
                identities_buffer.push(k);
                dist_buffer.push(v);
            }

            let weights_dist = WeightedIndex::new(&dist_buffer).unwrap();
            identity_samplers.push(Some((identities_buffer.to_owned(), weights_dist)));
        }
    }
    identity_samplers
}

fn read2kmer(
    read_seq: &[u8],
    ref_seq: &[u8],
    k: usize,
    error_transit: &mut [[[u32; 5]; 5]],
    match_len: &mut [[[u32; 62]; 5]],
) -> (impl Iterator<Item = (u32, ErrorPattern)>, f32) {
    assert!(k <= 16);
    let (left_half_k, right_half_k) = (k >> 1, k - (k >> 1));
    let config = EdlibAlignConfigRs {
        task: EdlibAlignTaskRs::EDLIB_TASK_PATH,
        ..Default::default()
    };
    let align_res = edlibAlignRs(ref_seq, read_seq, &config);
    let mut kmer_positions = Vec::with_capacity(ref_seq.len());

    // let mut read_aln_seq = Vec::with_capacity(read_seq.len());
    // let mut ref_aln_seq = Vec::with_capacity(ref_seq.len());
    // let mut ref_seq_i = 0;
    // let mut read_seq_i = 0;
    // if let Some(ref ops) = align_res.alignment {
    //  for op in ops.into_iter().copied() {
    //      if op == 0 || op == 3{
    //          ref_aln_seq.push(ref_seq[ref_seq_i]);
    //          ref_seq_i += 1;
    //          read_aln_seq.push(read_seq[read_seq_i]);
    //          read_seq_i += 1;
    //      }else if op == 1{
    //          ref_aln_seq.push(ref_seq[ref_seq_i]);
    //          ref_seq_i += 1;
    //          read_aln_seq.push(b'-');
    //      }else{
    //          ref_aln_seq.push(b'-');
    //          read_aln_seq.push(read_seq[read_seq_i]);
    //          read_seq_i += 1;
    //      }
    //  }
    // }
    // assert!(ref_seq_i == ref_seq.len() && read_seq_i == read_seq.len());
    // println!("read:{}\nref :{}\n", &String::from_utf8_lossy(&read_aln_seq), &String::from_utf8_lossy(&ref_aln_seq));

    let mut identity = 0.;
    if let Some(ops) = align_res.alignment {
        identity = ops.iter().filter(|x| **x == 0).count() as f32 / ops.len() as f32;
        if let Some((
            (ref_left_offset, ref_right_offset),
            (read_left_offset, _read_right_offset),
            (alignment_left_offset, alignment_right_offset),
        )) = trim_alignment_to_long_match(&ops, 8)
        {
            let mut ref_seq_i = ref_left_offset;
            let mut read_seq_i = read_left_offset;
            let ops = &ops[alignment_left_offset..ops.len() - alignment_right_offset];
            analyze_cigar_ops(ops, error_transit, match_len);

            let mut i = 0;
            while i < ops.len() {
                i += 1;
                read_seq_i += 1;
                ref_seq_i += 1;
                let (offset, ref_offset, read_offset) = extend_op_to_match(&ops[i..]);

                if left_half_k + ref_left_offset <= ref_seq_i - 1
                    && ref_seq_i <= ref_seq.len() - right_half_k - ref_right_offset
                {
                    kmer_positions.push((ref_seq_i - 1, ref_offset, read_seq_i, read_offset));
                }

                i += offset;
                ref_seq_i += ref_offset;
                read_seq_i += read_offset;
            }

            // assert!(
            //     ref_seq_i + ref_right_offset == ref_seq.len()
            //         && read_seq_i + read_right_offset == read_seq.len(),
            //     "{ref_seq_i} vs {} {read_seq_i} vs {}",
            //     ref_seq.len(),
            //     read_seq.len()
            // );
        }
    }

    (
        kmer_positions
            .into_iter()
            .map(move |(ref_seq_i, ref_offset, read_seq_i, read_offset)| {
                // println!("{} {} {ref_offset}", String::from_utf8_lossy(&ref_seq[ref_seq_i - left_half_k .. ref_seq_i + right_half_k]),
                // String::from_utf8_lossy(&read_seq[read_seq_i..(read_seq_i+read_offset)]));
                (
                    seq2kmer(&ref_seq[ref_seq_i - left_half_k..ref_seq_i + right_half_k]),
                    ErrorPattern::new(
                        &read_seq[read_seq_i..(read_seq_i + read_offset)],
                        ref_offset,
                    ),
                )
            }),
        identity,
    )
}

fn stat_error_rate(
    features: &[Vec<IntMap<u32, [u32; 2]>>],
    errors: &[ErrorPattern],
) -> (Vec<(usize, u64)>, f64) {
    let errors_num: Vec<(usize, usize)> = errors.iter().map(|x| x.to_error_num()).collect();

    let global_count = Mutex::new(vec![0u64; ErrorPattern::NUMBER_MAX]);
    features.par_iter().for_each(|vecs| {
        let mut local_count = vec![0u64; ErrorPattern::NUMBER_MAX];

        for ((n1, n2), (c1, c2)) in vecs.iter().flat_map(|x| {
            x.iter()
                .map(|(k, v)| (errors_num[*k as usize], (v[0], v[1])))
        }) {
            let n1_len = ErrorPattern::parse_error_num(n1).1;

            if n1_len > 0 {
                local_count[n1] += (c1 + c2) as u64;
            }
            if n2 > 0 {
                local_count[n2] += (c1 + c2) as u64;
            }

            if n1_len == 0 && n2 == 0 {
                local_count[0] += (c1 + c2) as u64;
            }
        }

        let mut global = global_count.lock().unwrap();
        for (i, val) in local_count.iter().enumerate() {
            global[i] += *val;
        }
    });

    let locked = global_count.lock().unwrap();
    let total = locked.iter().sum::<u64>();
    //skip positions where no errors occurred
    (
        locked
            .iter()
            .copied()
            .enumerate()
            .skip(1)
            .filter(|(_, x)| *x != 0)
            .collect(),
        1. - locked[0] as f64 / total as f64,
    )
}

fn parse_sa_tag(sa_tag: &str, header: &HeaderView) -> Vec<(usize, usize, u32, usize, usize)> {
    sa_tag
        .trim_end_matches(';')
        .split(';')
        .filter_map(|s| {
            let fields: Vec<&str> = s.split(',').collect();
            if fields.len() != 6 {
                return None;
            }
            let rid = header
                .tid(fields[0].as_bytes())
                .expect("Faield map header to tid!"); //safe
            let pos = fields[1].parse::<usize>().ok()? - 1;
            let std = fields[2];
            let (left_clip, right_clip, read_span, ref_span) = parse_cigar_span(fields[3]);
            if std == "+" {
                Some((left_clip, left_clip + read_span, rid, pos, pos + ref_span))
            } else {
                Some((right_clip, right_clip + read_span, rid, pos, pos + ref_span))
            }
        })
        .collect()
}

fn analyze_chimeric_feature(
    mut regions: Vec<(usize, usize, u32, usize, usize)>,
    min_span: usize,
) -> (usize, usize) {
    regions.sort();

    //self_chimeric is formed by + and - strands of DNA.
    let (mut dif_chimeric, mut self_chimeric) = (0, 0);
    for (p, &(read_start2, _, ref_name2, ref_start2, ref_end2)) in
        regions.iter().enumerate().skip(1)
    {
        let (_, read_end1, ref_name1, ref_start1, ref_end1) = regions[p - 1];
        let read_span = read_end1.abs_diff(read_start2);
        if read_span > min_span {
            continue;
        }

        if ref_name1 != ref_name2 {
            dif_chimeric += 1;
            continue;
        }

        let ref_start_span = ref_start1.abs_diff(ref_start2);
        let ref_end_span = ref_end1.abs_diff(ref_end2);

        if ref_start_span <= min_span && ref_end_span <= min_span {
            self_chimeric += 1;
            continue;
        }

        let ref_span = ref_end1.abs_diff(ref_start2);
        if ref_span > min_span {
            dif_chimeric += 1;
        }
    }

    (dif_chimeric, self_chimeric)
}

//TODO, re-use buffer
fn estimate_kmer_identity(seq: &[u8], kmer_dists: &[Vec<KmerDist>], ksize: usize) -> f32 {
    let mut buffer = Vec::with_capacity(seq.len());
    for kmer in iter2kmer(seq.iter().copied(), ksize).map(|x| x as usize) {
        let error_rate =
            kmer_dists[kmer & ((1 << COUNTER_BITS) - 1)][kmer >> COUNTER_BITS].error_rate;
        if error_rate > 0. {
            buffer.push(error_rate);
        }
    }

    buffer.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let cutoff = buffer.len() * 30 / 100;
    if cutoff > 0 {
        buffer.iter().take(cutoff).sum::<f32>() / cutoff as f32
    } else {
        0.0
    }
}

fn estimate_kmer_identitys<const N: usize>(
    kmer_features: &[Vec<IntMap<u32, [u32; 2]>>; N],
    errors: &[ErrorPattern],
    aln_regions: Vec<(i32, u32, u32)>,
    reads_identitys: Vec<f32>,
    ref_seqs: &[Vec<u8>],
    ksize: usize,
    use_len: bool,
) -> Vec<(f32, f32)> {
    if use_len {
        return aln_regions
            .iter()
            .map(|x| (x.2 - x.1) as f32)
            .zip(reads_identitys)
            .collect();
    }

    let kmer_dists = &build_kmer_features(kmer_features, errors);
    let mut kmers_errors = vec![0.; aln_regions.len()];
    kmers_errors
        .par_iter_mut()
        .enumerate()
        .for_each(|(p, ide)| {
            let (tid, start, end) = &aln_regions[p];
            *ide = estimate_kmer_identity(
                &ref_seqs[*tid as usize][*start as usize..*end as usize],
                kmer_dists,
                ksize,
            );
        });

    // for (kide, ide) in kmers_errors.iter().zip(reads_identitys.iter()){
    //     println!("{kide}\t{ide}");
    // }
    // panic!("");

    kmers_errors.into_iter().zip(reads_identitys).collect()
}

fn main() {
    let opt = &Opt::from_args();
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.thread)
        .build_global()
        .unwrap();

    if opt.cmd == "rhq" {
        let mut bam = bam::Reader::from_path(&opt.bam).unwrap();
        let header = bam::Header::from_template(bam.header());
        let mut out = bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap();
        bam.set_threads(opt.thread)
            .expect("Failed enable multi-threaded BAM reading!");
        out.set_threads(opt.thread)
            .expect("Failed enable multi-threaded BAM writing!");

        let mut alns = HashMap::default();
        let mut record = Record::new();
        while let Some(ret) = bam.read(&mut record) {
            ret.expect("Failed parse BAM!");

            if record.is_unmapped() || record.mapq() > 0 {
                continue;
            }

            let is_primary = !record.is_secondary() && !record.is_supplementary();

            let read_len = if is_primary {
                record.seq_len()
            } else {
                record.seq_len_from_cigar(true)
            } as u32;
            let (read_start, read_end) = get_read_aligned_range(&record.cigar(), read_len);
            if read_end - read_start < (read_len as f32 * opt.map_fra) as usize {
                continue;
            }

            let as_score = match record.aux(b"AS") {
                Ok(Aux::U8(v)) => v as i64,
                Ok(Aux::U16(v)) => v as i64,
                Ok(Aux::U32(v)) => v as i64,
                Ok(Aux::I8(v)) => v as i64,
                Ok(Aux::I16(v)) => v as i64,
                Ok(Aux::I32(v)) => v as i64,
                Ok(other) => panic!("unexpected AS tag type: {other:?}"),
                Err(_) => i64::MIN,
            };

            let aln = alns.entry(record.qname().to_owned()).or_insert(Vec::new());
            aln.push((
                is_primary,
                as_score,
                record.tid() as u32,
                record.reference_start() as u32,
                record.reference_end() as u32,
            ));
        }

        let ref_seqs = &read_seqs_with_tid(&opt.genome, bam.header());

        //Primary and secondary alignments share score and reference
        alns.retain(|_, v| {
            if let Some(primary_index) = v.iter().position(|item| item.0) {
                let (_, score1, tid1, start1, end1) = v[primary_index];
                let seq1 = &ref_seqs[tid1 as usize][start1 as usize..end1 as usize];

                let (dup_count, dif_count) = v
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != primary_index)
                    .fold(
                        (0, 0),
                        |(dup, dif), (_, &(_, score2, tid2, start2, end2))| {
                            if score2 > score1 {
                                (dup, dif + 1)
                            } else if score2 == score1 {
                                let seq2 = &ref_seqs[tid2 as usize][start2 as usize..end2 as usize];
                                if seq1 == seq2 {
                                    (dup + 1, dif)
                                } else {
                                    (dup, dif + 1)
                                }
                            } else {
                                (dup, dif)
                            }
                        },
                    );

                if dup_count > 0 && dif_count == 0 {
                    *v = vec![v[primary_index]];
                    true
                } else {
                    false
                }
            } else {
                false
            }
        });

        bam = bam::Reader::from_path(&opt.bam).unwrap();
        while let Some(ret) = bam.read(&mut record) {
            ret.expect("Failed parse BAM!");
            if let Some(aln) = alns.get(record.qname()) {
                if record.tid() as u32 == aln[0].2
                    && record.reference_start() as u32 == aln[0].3
                    && record.reference_end() as u32 == aln[0].4
                {
                    record.set_mapq(61);
                }
            }
            out.write(&record).unwrap();
        }
    } else if opt.cmd == "sim" {
        let ref_seqs = &read_seqs_with_names(&opt.genome);
        let (
            read_type,
            ksize,
            strand_fra,
            chimeric_fra,
            kmer_features,
            errors,
            error_transits,
            _match_lens,
            emp_aligned_len,
            emp_unaligned_left_len,
            emp_unaligned_right_len,
            errors_rate,
        ) = &load_file(&opt.mode, opt.thread);

        let (global_errors, global_has_error) = &stat_error_rate(kmer_features, errors);
        let depth = if opt.base > 0 {
            opt.base as f64 / ref_seqs.iter().map(|x| x.1.len()).sum::<usize>() as f64
        } else {
            opt.depth as f64
        };
        info!("Model loading completed, kmer error rate: {global_has_error:.4}");

        // output_kmer_features(kmer_features, errors);

        let average_error = *global_has_error;
        let strand_dist = &Bernoulli::new(*strand_fra).unwrap();
        let chimeric_dist = &WeightedIndex::new(chimeric_fra).unwrap();
        let global_error_dist =
            &WeightedIndex::new(global_errors.iter().map(|item| item.1)).unwrap();
        let global_has_error = &Bernoulli::new(*global_has_error).unwrap();
        let has_error = &Bernoulli::new(UNALIGNED_ERROR_RATE).unwrap();
        let random_error = &Bernoulli::new(opt.noise).unwrap();
        let error_transits = &build_distribution(error_transits, vec![0, 4]);
        // let match_lens = &build_distribution(match_lens, vec![]);
        let kmer_features = &build_kmer_features(kmer_features, errors);
        let kmer_identitys = &build_kmer_identity_distribution(errors_rate, read_type == "hifi", opt.global_error_rate);

        thread::scope(|work| {
            let (ou_s, ou_r) = bounded(opt.thread + 2);
            for i in 0..opt.thread as u64 {
                let ou_s = ou_s.clone();
                work.spawn(move || {
                    let seed = if opt.seed == 0 {
                        rand::rng().random()
                    } else {
                        opt.seed << 12 | i
                    };
                    let mut rng = StdRng::seed_from_u64(seed);

                    let (mut seq_buffer, mut cigar_buffer, mut qual_buffer) =
                        (Vec::new(), Vec::new(), Vec::new());
                    let (mut bam_index, mut bam_buffer) = (1, String::new());
                    let mut sim_reads = SimReads::new();
                    for (tid, (ref_name, ref_seq)) in ref_seqs.iter().enumerate() {
                        let ref_len = ref_seq.len();
                        let total_bases = (ref_len as f64 * depth) as usize / opt.thread + 1;
                        let mut sim_bases = 0;
                        let pos_dist = Uniform::new(0, ref_len).unwrap();
                        let mut chimeric_reads = Vec::with_capacity(4);
                        let (mut chimeric, mut aligned_len, mut unaligned_left_len) = (0, 0, 0);
                        let (mut unaligned_right_len, mut start_pos, mut is_positive) = (0, 0, true);
                        while sim_bases <= total_bases {
                            if chimeric != 4 || chimeric_reads.is_empty(){//random chimeric
                                aligned_len = emp_aligned_len.sample(
                                    &mut rng,
                                    if opt.median_len == 0 {
                                        None
                                    } else {
                                        Some(opt.median_len)
                                    },
                                );
                                unaligned_left_len = emp_unaligned_left_len.sample(&mut rng, None);
                                unaligned_right_len = emp_unaligned_right_len.sample(&mut rng, None);
                                start_pos = pos_dist.sample(&mut rng) + (ksize >> 1);
                                if chimeric == 0 && !opt.disable_chimeric {
                                    chimeric = chimeric_dist.sample(&mut rng).min(3);
                                }
                                is_positive = strand_dist.sample(&mut rng);
                            }else{//self chimeric
                                is_positive = !is_positive;
                                chimeric = 2;
                            }

                            let total_len = (unaligned_left_len + aligned_len + unaligned_right_len) as usize;
                            if start_pos + total_len >= ref_len {
                                continue;
                            }

                            let (seq_start, cigar_start) = (seq_buffer.len(), cigar_buffer.len());

                            let mut end_pos = start_pos;
                            if !opt.disable_unaligned
                                && end_pos + (unaligned_left_len as usize) < ref_len
                            {
                                end_pos += sim_unaligned_seq(
                                    &mut rng,
                                    unaligned_left_len as usize,
                                    &ref_seq[end_pos..],
                                    global_errors,
                                    global_error_dist,
                                    has_error,
                                    &mut seq_buffer,
                                    &mut cigar_buffer,
                                    &mut qual_buffer,
                                );
                            }

                            if end_pos + (aligned_len as usize) < ref_len {
                                end_pos += sim_aligned_seq(
                                    &mut rng,
                                    aligned_len as usize,
                                    &ref_seq[end_pos - (ksize >> 1)..],
                                    global_errors,
                                    global_error_dist,
                                    global_has_error,
                                    random_error,
                                    error_transits,
                                    *ksize,
                                    kmer_features,
                                    errors,
                                    is_positive,
                                    kmer_identitys,
                                    opt.temperature,
                                    average_error,
                                    read_type == "hifi",
                                    opt.global_error_rate,
                                    &mut seq_buffer,
                                    &mut cigar_buffer,
                                    &mut qual_buffer,
                                ) - (ksize >> 1);
                            }

                            if !opt.disable_unaligned
                                && end_pos + (unaligned_right_len as usize) < ref_len
                            {
                                end_pos += sim_unaligned_seq(
                                    &mut rng,
                                    unaligned_right_len as usize,
                                    &ref_seq[end_pos..],
                                    global_errors,
                                    global_error_dist,
                                    has_error,
                                    &mut seq_buffer,
                                    &mut cigar_buffer,
                                    &mut qual_buffer,
                                );
                            }

                            if seq_buffer.len() < opt.min_len as usize
                                || aligned_len < opt.min_len
                                || seq_buffer.len() > opt.max_len as usize
                                || aligned_len > opt.max_len
                            {   
                                seq_buffer.clear();
                                cigar_buffer.clear();
                                qual_buffer.clear();
                                chimeric_reads.clear();
                                continue;
                            }

                            // println!(">{is_positive} {start_pos}\n{}\n{}\n{}",
                            // std::str::from_utf8(&seq_buffer).unwrap(),
                            // std::str::from_utf8(&cigar_buffer).unwrap(),
                            // std::str::from_utf8(&qual_buffer.iter().map(|x| x+33).collect::<Vec<u8>>()).unwrap());
                            if chimeric > 0 {
                                chimeric_reads.push((
                                    seq_start,
                                    cigar_start,
                                    ref_name.as_str(),
                                    start_pos,
                                    end_pos,
                                    unaligned_left_len,
                                    aligned_len,
                                    unaligned_right_len,
                                    tid as i32,
                                    is_positive,
                                    chimeric,
                                    i,
                                    bam_index),
                                );

                                if chimeric_reads.len() >= chimeric {
                                    sim_reads.push_chimeric(&chimeric_reads, ref_seqs, &mut seq_buffer, &mut cigar_buffer, &mut qual_buffer, &mut bam_buffer);
                                    sim_bases += seq_buffer.len();
                                    seq_buffer.clear();
                                    cigar_buffer.clear();
                                    qual_buffer.clear();
                                    chimeric_reads.clear();
                                    chimeric = 0;
                                }
                            }else{
                                sim_reads.push(
                                    (&seq_buffer,
                                    &cigar_buffer,
                                    &qual_buffer,
                                    ref_name,
                                    start_pos,
                                    end_pos,
                                    unaligned_left_len,
                                    aligned_len,
                                    unaligned_right_len,
                                    tid as i32,
                                    is_positive,
                                    chimeric,
                                    i,
                                    bam_index),
                                    &mut bam_buffer,
                                );
                                sim_bases += seq_buffer.len();
                                seq_buffer.clear();
                                cigar_buffer.clear();
                                qual_buffer.clear();
                            }

                            bam_index += 1;
                            if sim_reads.is_full() {
                                ou_s.send(sim_reads).unwrap();
                                sim_reads = SimReads::new();
                            }
                        }
                    }
                    ou_s.send(sim_reads).unwrap();
                });
            }
            drop(ou_s);

            work.spawn(move || {
                let (mut count, mut bases) = (0, 0);
                let header = create_bam_header(ref_seqs);
                let mut writer = Writer::from_stdout(&header, bam::Format::Bam).unwrap();
                writer.set_threads(opt.thread.min(4)).unwrap();
                while let Ok(sim_reads) = ou_r.recv() {
                    for read in sim_reads.iter() {
                        writer.write(read).expect("Failed to write read!");
                        count += 1;
                        bases += read.seq_len();
                        if count % 10000 == 0 {
                            info!("Simulated reads: {count}, bases: {bases} bp");
                        }
                    }
                }
                info!("Simulated reads: {count}, bases: {bases} bp");
            });
        })
    } else {
        let mut bam = bam::Reader::from_path(&opt.bam).unwrap();
        let header = bam.header();

        let ref_seqs = &read_seqs_with_tid(&opt.genome, header);

        thread::scope(|work| {
            let (in_s, in_r) = bounded(1);
            let result1 = work.spawn(move || {
                let mut strands_stat = [0, 0]; //+ 
                let mut chimeric_count = vec![0; 5];
                let mut aligned_len = Vec::with_capacity(65536);
                let mut unaligned_len = [Vec::with_capacity(65536), Vec::with_capacity(65536)];
                let mut read_seqs = Vec::with_capacity(opt.batch_size + 65536);
                let mut read_seqs_info = Vec::with_capacity(opt.batch_size / 16384 + 1);

                let mut reads_idx = 0;
                let mut total_bases = 0;
                let mut reads_aln_regions = Vec::with_capacity(opt.batch_size / 16384 + 1);

                let mut record = Record::new();
                let header = bam.header().to_owned();
                while let Some(ret) = bam.read(&mut record) {
                    ret.expect("Failed parse BAM!");
                    if record.is_secondary()
                        || record.is_supplementary()
                        || record.is_unmapped()
                        || record.mapq() < opt.mapq
                    {
                        continue;
                    }

                    let ref_start = record.reference_start() as usize;
                    let ref_end = record.reference_end() as usize;

                    //here we don't consider strand
                    let (read_start, read_end) =
                        get_read_aligned_range(&record.cigar(), record.seq_len() as u32);
                    aligned_len.push((read_end - read_start) as u32);

                    if let Ok(sa_aux) = record.aux(b"SA") {
                        if let Aux::String(sa_str) = sa_aux {
                            let mut aligned_regions = parse_sa_tag(sa_str, &header);
                            aligned_regions.push((
                                read_start,
                                read_end,
                                record.tid() as u32,
                                ref_start,
                                ref_end,
                            ));
                            let (dif_chimeric, self_chimeric) =
                                analyze_chimeric_feature(aligned_regions, 30);
                            if dif_chimeric != 0 {
                                chimeric_count[dif_chimeric.min(3)] += 1;
                            }
                            chimeric_count[4] += self_chimeric;
                        }
                    } else {
                        //don't record clip len if record has supplementary alignments
                        //read_start, read_end need to reverse for "-" strand
                        if record.is_reverse() {
                            unaligned_len[0].push((record.seq_len() - read_end) as u32);
                            unaligned_len[1].push(read_start as u32);
                        } else {
                            unaligned_len[0].push(read_start as u32);
                            unaligned_len[1].push((record.seq_len() - read_end) as u32);
                        }
                        chimeric_count[0] += 1;
                    }

                    if read_end - read_start
                        < std::cmp::max(
                            opt.map_len,
                            (record.seq_len() as f32 * opt.map_fra) as usize,
                        )
                    {
                        continue;
                    }

                    let std = if record.is_reverse() {
                        strands_stat[1] += 1;
                        1
                    } else {
                        strands_stat[0] += 1;
                        0
                    };

                    let read_seq = record.seq();
                    let read_seq_offset = read_seqs.len();
                    read_seqs.extend((read_start..read_end).map(|i| read_seq[i]));
                    read_seqs_info.push((
                        read_seq_offset,
                        read_seqs.len(),
                        record.tid() as usize,
                        ref_start,
                        ref_end,
                        std,
                    ));

                    reads_aln_regions.push((record.tid(), ref_start as u32, ref_end as u32));
                    reads_idx += 1;

                    total_bases += read_seq.len();
                    if read_seqs.len() >= opt.batch_size {
                        info!("Processing {total_bases} bases!");
                        in_s.send((
                            reads_idx - read_seqs_info.len(),
                            read_seqs.to_owned(),
                            read_seqs_info.to_owned(),
                        ))
                        .expect("Failed to send data.");
                        read_seqs.clear();
                        read_seqs_info.clear();
                    }
                }

                if !read_seqs.is_empty() {
                    info!("Processing {total_bases} bases!");
                    in_s.send((reads_idx - read_seqs_info.len(), read_seqs, read_seqs_info))
                        .expect("Failed to send data.");
                }

                let total = chimeric_count.iter().sum::<usize>() as f64;
                let chimeric_fra = chimeric_count
                    .into_iter()
                    .map(|x| x as f64 / total)
                    .collect::<Vec<f64>>();
                (
                    strands_stat[0] as f64 / strands_stat.iter().sum::<usize>() as f64,
                    chimeric_fra,
                    EmpHistMode::from_lengths(&mut aligned_len, 500),
                    EmpHistMode::from_lengths(&mut unaligned_len[0], 1),
                    EmpHistMode::from_lengths(&mut unaligned_len[1], 1),
                    reads_aln_regions,
                )
            });

            let result2 = work.spawn(move || {
                let mut kmer_features: [Vec<IntMap<u32, [u32; 2]>>; 1 << COUNTER_BITS] =
                    std::array::from_fn(|_| {
                        vec![IntMap::default(); 1 << (2 * opt.k - COUNTER_BITS)]
                    });
                let mut match_lens = vec![[[0u32; 62]; 5]; 21];
                let mut error_transits = vec![[[0u32; 5]; 5]; 22];
                let mut reads_identitys = vec![0.; 65536];
                let interner = Arc::new(ErrorPatternInterner::new());

                while let Ok((reads_idx, read_seqs, seqs_info)) = in_r.recv() {
                    let (read_seqs, seqs_info) = (&read_seqs, &seqs_info);
                    let step = seqs_info.len().div_ceil(opt.thread);

                    let result = thread::scope(|kmer| {
                        let mut handles = Vec::with_capacity(opt.thread);
                        for i in 0..opt.thread {
                            let interner = interner.clone();
                            let (s, e) = (i * step, ((i + 1) * step).min(seqs_info.len()));
                            if s >= e {
                                continue;
                            }

                            let handle = kmer.spawn(move || {
                                let mut kmers: [Vec<(u32, u32, u32)>; 1 << COUNTER_BITS] =
                                    std::array::from_fn(|_| {
                                        Vec::with_capacity(
                                            (opt.batch_size >> COUNTER_BITS) / opt.thread + 1024,
                                        )
                                    });
                                let mut match_len = vec![[[0u32; 62]; 5]; 21];
                                let mut error_transit = vec![[[0u32; 5]; 5]; 22];
                                let mut identitys = Vec::with_capacity(e - s);

                                for (p, (read_start, read_end, ref_id, ref_start, ref_end, std)) in
                                    seqs_info[s..e].iter().copied().enumerate()
                                {
                                    let (read_kmers, identity) = read2kmer(
                                        &read_seqs[read_start..read_end],
                                        &ref_seqs[ref_id][ref_start..ref_end],
                                        opt.k,
                                        &mut error_transit,
                                        &mut match_len,
                                    );

                                    for (kmer, error) in read_kmers {
                                        let error_sym = interner.get_or_insert(error);
                                        kmers[kmer as usize & ((1 << COUNTER_BITS) - 1)].push((
                                            kmer >> COUNTER_BITS,
                                            error_sym,
                                            std,
                                        ));
                                    }
                                    identitys.push((reads_idx + s + p, identity));
                                }
                                (kmers, error_transit, match_len, identitys)
                            });
                            handles.push(handle);
                        }
                        handles
                            .into_iter()
                            .map(|h| h.join().unwrap())
                            .collect::<Vec<_>>()
                    });

                    //(kmer_errors, error_transit, match_len)
                    kmer_features
                        .par_iter_mut()
                        .enumerate()
                        .for_each(|(p, vec)| {
                            for (kmer, error_sym, std) in result.iter().flat_map(|x| &x.0[p]) {
                                vec[*kmer as usize].entry(*error_sym).or_insert([0, 0])
                                    [*std as usize] += 1;
                            }
                        });

                    for (error_transit, match_len, identitys) in
                        result.into_iter().map(|x| (x.1, x.2, x.3))
                    {
                        sum_3d(&mut match_lens, &match_len);
                        sum_3d(&mut error_transits, &error_transit);

                        let max_index = identitys[identitys.len() - 1].0;
                        if reads_identitys.len() <= max_index {
                            reads_identitys.resize(max_index + 1, 0.);
                        }

                        for (index, identity) in identitys.into_iter() {
                            reads_identitys[index] = identity;
                        }
                    }
                }

                (
                    kmer_features,
                    interner,
                    error_transits,
                    match_lens,
                    reads_identitys,
                )
            });

            let (
                strand_fra,
                chimeric_fra,
                emp_aligned_len,
                emp_unaligned_left_len,
                emp_unaligned_right_len,
                reads_aln_regions,
            ) = result1.join().unwrap();
            let (kmer_features, interner, error_transits, match_lens, reads_identitys) =
                result2.join().unwrap();
            let kmer_identitys = estimate_kmer_identitys(
                &kmer_features,
                &interner.vec.read(),
                reads_aln_regions,
                reads_identitys,
                ref_seqs,
                opt.k,
                opt.read_type == "hifi",
            );

            dump_file(
                &opt.read_type,
                opt.k,
                opt.thread as u32,
                strand_fra,
                &chimeric_fra,
                &kmer_features,
                &interner,
                &error_transits,
                &match_lens,
                &emp_aligned_len,
                &emp_unaligned_left_len,
                &emp_unaligned_right_len,
                &kmer_identitys,
            );
        });
    }

    info!("\n{}", resource_str());
}
