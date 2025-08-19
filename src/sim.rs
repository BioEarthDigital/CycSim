use crate::{
    COUNTER_BITS, DELETION_CIGAR, ErrorPattern, INSERTION_CIGAR, KMER_ERROR_BIN, KmerDist,
    LEN_ERROR_BIN, MISMATCH_CIGAR, TRANSIT_WIN, UNALIGNED_ERROR_RATE, estimate_kmer_identity,
    iter2kmer, qual::identity_to_qual,
};
use rand::{
    distr::{Bernoulli, Uniform, weighted::WeightedIndex},
    prelude::*,
};
use rand_distr::Distribution;
use std::sync::LazyLock;

static BASES: &[u8] = b"ATGC";
static UNIFORM_BASE: LazyLock<Uniform<usize>> =
    LazyLock::new(|| Uniform::new(0, BASES.len()).expect("Failed to create Uniform distribution"));

fn iter_random_base(rng: &mut StdRng) -> impl Iterator<Item = u8> {
    std::iter::from_fn(move || Some(BASES[UNIFORM_BASE.sample(rng)]))
}

fn random_error_pattern(rng: &mut StdRng, ref_skip_seq: &[u8], seq_sym_len: usize) -> ErrorPattern {
    let mut i = 0;
    let mut seq_sym = Vec::with_capacity(seq_sym_len);
    if seq_sym_len > 0 {
        for base in iter_random_base(rng) {
            if i < ref_skip_seq.len() {
                if base != ref_skip_seq[i] {
                    seq_sym.push(base);
                    i += 1;
                }
            } else {
                seq_sym.push(base);
            }

            if seq_sym.len() >= seq_sym_len {
                break;
            }
        }
    }

    ErrorPattern::new(&seq_sym, ref_skip_seq.len())
}

fn extend_seq_cigar_with_random_error(
    seq: &mut Vec<u8>,
    cigar: &mut Vec<u8>,
    ref_seq: &[u8],
    rng: &mut StdRng,
    error_type: usize,
    error_len: usize,
) -> (usize, usize) {
    let (mut i, mut l) = (0, 0);
    let error_len = error_len.min(ref_seq.len());
    match error_type {
        MISMATCH_CIGAR => {
            //mismatch
            let mut count = 0;
            for base in iter_random_base(rng) {
                if base != ref_seq[i] {
                    seq.push(base);
                    cigar.push(b'X');
                    count += 1;
                    i += 1;
                    l += 1;
                    if count >= error_len {
                        break;
                    }
                }
            }
        }
        DELETION_CIGAR => {
            //DELETION
            i += error_len;
            cigar.extend(std::iter::repeat_n(b'D', error_len));
        }
        INSERTION_CIGAR => {
            //INSERTION
            seq.extend(iter_random_base(rng).take(error_len));
            cigar.extend(std::iter::repeat_n(b'I', error_len));
            l += error_len;
        }
        _ => panic!("Unknown error type: {error_type}"),
    }
    (i, l)
}

fn extend_seq_cigar_with_error_pattern(
    seq: &mut Vec<u8>,
    cigar: &mut Vec<u8>,
    error: &ErrorPattern,
) -> usize {
    seq.extend(&error.seq_sym);

    if error.is_deletion() {
        if !error.seq_sym.is_empty() {
            cigar.extend(std::iter::repeat_n(b'X', error.seq_sym.len()));
        }
        cigar.extend(std::iter::repeat_n(
            b'D',
            error.ref_skip - error.seq_sym.len(),
        ));
        error.ref_skip
    } else if error.is_insertion() {
        if error.ref_skip > 0 {
            cigar.extend(std::iter::repeat_n(b'X', error.ref_skip));
        }
        cigar.extend(std::iter::repeat_n(
            b'I',
            error.seq_sym.len() - error.ref_skip,
        ));
        error.seq_sym.len()
    } else {
        assert_eq!(error.seq_sym.len(), error.ref_skip);
        cigar.extend(std::iter::repeat_n(b'X', error.seq_sym.len()));
        error.seq_sym.len()
    }
}

// pub fn sim_unaligned_seq(
//     rng: &mut StdRng,
//     len: usize,
//     seq: &mut Vec<u8>,
//     cigar: &mut Vec<u8>,
//     qual: &mut Vec<u8>,
// ) -> usize {
//     seq.extend(iter_random_base(rng).take(len));
//     cigar.extend(std::iter::repeat_n(b'S', len));
//     qual.extend(std::iter::repeat_n(MIN_QUAL, len));
//     0
// }

pub fn sim_unaligned_seq(
    rng: &mut StdRng,
    len: usize,
    ref_seq: &[u8],
    errors: &[(usize, u64)],
    error_dist: &WeightedIndex<u64>,
    has_error: &Bernoulli,
    seq: &mut Vec<u8>,
    cigar: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) -> usize {
    let mut l = 0;
    let mut i = 0;
    while l < len && i < ref_seq.len() {
        if ref_seq[i] == b'N' {
            break;
        }
        if has_error.sample(rng) {
            let idx = error_dist.sample(rng);
            let (error_type, error_len) = ErrorPattern::parse_error_num(errors[idx].0);
            let (ref_skip, extend_len) = extend_seq_cigar_with_random_error(
                seq,
                cigar,
                &ref_seq[i..],
                rng,
                error_type,
                error_len,
            );
            i += ref_skip; //may exceed than ref_seq length
            l += extend_len;
        } else {
            seq.push(ref_seq[i]);
            cigar.push(b'=');
            i += 1;
            l += 1;
        }
    }
    let qs = identity_to_qual(UNALIGNED_ERROR_RATE);
    qual.extend(std::iter::repeat_n(qs, seq.len() - qual.len()));

    i
}

fn sim_aligned_seq_with_temperature(
    rng: &mut StdRng,
    len: usize,
    ref_seq: &[u8],
    global_errors: &[(usize, u64)],
    global_error_dist: &WeightedIndex<u64>,
    global_has_error: &Bernoulli,
    random_error: &Bernoulli,
    error_transits: &[Vec<WeightedIndex<u32>>],
    ksize: usize,
    kmer_dists: &[Vec<KmerDist>],
    errors: &[ErrorPattern],
    is_positive: bool,
    end_temp: f64,
    average_error: f64,
    seq: &mut Vec<u8>,
    cigar: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) -> (usize, usize) {
    //no error in the leftmost half of the kmer
    let mut i = ksize >> 1;
    let mut match_count = 0; //match cigar count
    let mut l = 0;

    let mut last_op = 4; //start cigar by default
    let mut last_error = 4; //start cigar by default
    let mut last_len = 0;
    let mut del_len = 0;
    let start_temp = end_temp * 0.93;
    let bias_len = 5000;

    for kmer in iter2kmer(ref_seq.iter().copied(), ksize).map(|x| x as usize) {
        if del_len > 0 {
            del_len -= 1;
            continue;
        }

        //the error rate decreases as the length of the read increases
        let temperature = if l < bias_len {
            start_temp + (end_temp - start_temp) * (i as f64 / bias_len as f64)
        } else {
            end_temp
        };

        //every incorrect kmer starts with a match base
        seq.push(ref_seq[i]);
        cigar.push(b'=');
        i += 1;
        l += 1;
        match_count += 1;

        if last_op != 0 {
            last_len = 0;
            last_op = 0;
        }
        last_len += 1;

        let kmer_dist = &kmer_dists[kmer & ((1 << COUNTER_BITS) - 1)][kmer >> COUNTER_BITS];
        let mut error_buffer = None;
        let (mut error, pro) = kmer_dist
            .sample(is_positive, KmerDist::ALL, rng, errors, temperature)
            .unwrap_or_else(|| {
                //the kmer does not exist in the database
                let has_error = global_has_error.sample(rng);
                if has_error {
                    let idx = global_error_dist.sample(rng);
                    let (error_type, error_len) =
                        ErrorPattern::parse_error_num(global_errors[idx].0);
                    error_buffer = Some(match error_type {
                        DELETION_CIGAR => random_error_pattern(rng, &ref_seq[i..i + error_len], 0),
                        INSERTION_CIGAR => random_error_pattern(rng, &ref_seq[i..i], error_len),
                        MISMATCH_CIGAR => {
                            random_error_pattern(rng, &ref_seq[i..i + error_len], error_len)
                        }
                        _ => panic!("Unknown error_type value: {error_type}"),
                    });
                } else {
                    error_buffer = Some(random_error_pattern(rng, &[], 0));
                }
                (error_buffer.as_ref().unwrap(), 1. - average_error)
            });

        if !error.is_match() {
            let bin1 = (last_len / TRANSIT_WIN).min(20);
            let error_transit = if last_op == 0 {
                error_transits[bin1 + 1][last_error].sample(rng)
            } else {
                error_transits[0][last_error].sample(rng)
            };

            if random_error.sample(rng) {
                let idx = global_error_dist.sample(rng);
                let (mut error_type, mut error_len) =
                    ErrorPattern::parse_error_num(global_errors[idx].0);

                if error_transit == last_error {
                    while error_type != error_transit {
                        let idx = global_error_dist.sample(rng);
                        let (et, el) = ErrorPattern::parse_error_num(global_errors[idx].0);
                        error_type = et;
                        error_len = el;
                    }
                }
                let (ref_skip, extend_len) = extend_seq_cigar_with_random_error(
                    seq,
                    cigar,
                    &ref_seq[i..],
                    rng,
                    error_type,
                    error_len,
                );
                i += ref_skip;
                l += extend_len;
                del_len = ref_skip;

                if last_op != error_type {
                    last_len = 0;
                    last_op = error_type;
                }
                last_len += error_len;
                last_error = error_type;
            } else {
                //use HMM mode
                if error_transit == last_error {
                    //a error may belong to both mismatch and indel.
                    if (error_transit == MISMATCH_CIGAR && !error.is_mismatch())
                        || ((error_transit == DELETION_CIGAR || error_transit == INSERTION_CIGAR)
                            && error.error_type() != error_transit)
                    {
                        //deletion
                        if let Some((e, _)) =
                            kmer_dist.sample(is_positive, error_transit, rng, errors, temperature)
                        {
                            error = e;
                        }
                    }
                }

                if i + error.ref_skip >= ref_seq.len() {
                    let qs = identity_to_qual(pro);
                    qual.extend(std::iter::repeat_n(qs, seq.len() - qual.len()));
                    break;
                }

                let op_len = extend_seq_cigar_with_error_pattern(seq, cigar, error);
                i += error.ref_skip;
                l += error.seq_sym.len();
                del_len = error.ref_skip;

                if last_op != error.error_type() {
                    last_len = 0;
                    last_op = error.error_type();
                }
                last_len += op_len;
                last_error = error.error_type();
            }
        }

        let qs = identity_to_qual(pro);
        qual.extend(std::iter::repeat_n(qs, seq.len() - qual.len()));

        if l >= len {
            break;
        }
    }

    (i, match_count)
}

pub fn sim_aligned_seq(
    rng: &mut StdRng,
    len: usize,
    ref_seq: &[u8],
    global_errors: &[(usize, u64)],
    global_error_dist: &WeightedIndex<u64>,
    global_has_error: &Bernoulli,
    random_error: &Bernoulli,
    error_transits: &[Vec<WeightedIndex<u32>>],
    ksize: usize,
    kmer_dists: &[Vec<KmerDist>],
    errors: &[ErrorPattern],
    is_positive: bool,
    kmer_identitys: &[Option<(Vec<u32>, WeightedIndex<u32>)>],
    opt_temperature: f64,
    average_error: f64,
    use_len: bool,
    seq: &mut Vec<u8>,
    cigar: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) -> usize {
    let bin_idx = if use_len {
        (len as f32 / LEN_ERROR_BIN).floor() as usize
    } else {
        let kmer_identity = estimate_kmer_identity(&ref_seq[..(len + ksize).min(ref_seq.len())], kmer_dists, ksize);
        (kmer_identity / KMER_ERROR_BIN).floor() as usize
    };
    let expected_identity = kmer_identitys
        .get(bin_idx)
        .and_then(|opt| opt.as_ref())
        .map(|(vals, dist)| vals[dist.sample(rng)] as f64 / 10000.0)
        .unwrap_or(1. - average_error);

    let (seq_start, cigar_start) = (seq.len(), cigar.len());

    let mut low = 0.01;
    let mut high = 3.0;
    let mut temperature = if opt_temperature == 0. {
        // if expected_identity >= 1. {
        //     seq.extend(&ref_seq[..len]);
        //     cigar.extend(std::iter::repeat_n(b'=', len));
        //     qual.extend(std::iter::repeat_n(MAX_QUAL, len));
        //     return len;
        // }
        2.0
    } else {
        opt_temperature
    };

    let (max_error_tolerate, max_sim_count) = if average_error < 0.01 {
        (0.0005, 100)
    } else {
        (0.001, 10)
    };

    let mut sim_count = 0;
    loop {
        seq.truncate(seq_start);
        cigar.truncate(cigar_start);
        qual.truncate(seq_start);
        let (i, match_count) = sim_aligned_seq_with_temperature(
            rng,
            len,
            ref_seq,
            global_errors,
            global_error_dist,
            global_has_error,
            random_error,
            error_transits,
            ksize,
            kmer_dists,
            errors,
            is_positive,
            temperature,
            average_error,
            seq,
            cigar,
            qual,
        );

        let seq_identity = match_count as f64 / cigar[cigar_start..].len() as f64;
        if expected_identity == 0.
            || (expected_identity - seq_identity).abs() <= max_error_tolerate
            || opt_temperature != 0.
        {
            // println!("{sim_count}\t{seq_identity}\t{expected_identity}\t{temperature}\t{low}\t{high}\t{}", seq_identity - expected_identity);
            return i;
        }

        if seq_identity < expected_identity {
            if high - temperature < 0.1 {
                high += 1.;
            }
            low = temperature;
        } else {
            if temperature - low < 0.1 {
                low = (low - 1.0).max(0.00001);
            }
            high = temperature;
        }
        temperature = (low + high) / 2.0;

        sim_count += 1;
        if sim_count >= max_sim_count {
            return i;
        }
    }
}
