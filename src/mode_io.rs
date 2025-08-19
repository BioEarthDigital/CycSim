#![allow(dead_code)]
use crate::{COUNTER_BITS, EmpHistMode, ErrorPattern, ErrorPatternInterner, FileIO, IntMap};
use byteorder::{ByteOrder, LittleEndian, ReadBytesExt, WriteBytesExt};
use crossbeam_channel::bounded;
use std::{
    io::{self, BufWriter, Read, Write},
    thread,
};
use zstd::stream::{Decoder, Encoder};

const CFSIM_MAGIC: &[u8] = "CF02\u{1}".as_bytes();

fn dump_kmer_features<W: Write>(w: &mut W, data: &[Vec<IntMap<u32, [u32; 2]>>]) {
    w.write_u64::<LittleEndian>(data.len() as u64).unwrap();
    for vec in data.iter() {
        w.write_u64::<LittleEndian>(vec.len() as u64).unwrap();
        w.write_u64::<LittleEndian>(vec.iter().map(|x| x.len()).sum::<usize>() as u64)
            .unwrap(); //len of map in this vec
        for map in vec.iter() {
            w.write_u16::<LittleEndian>(
                map.len().try_into().expect("Error feature map too large!"),
            )
            .unwrap();
            for (k, v) in map.iter() {
                w.write_u32::<LittleEndian>(*k).unwrap();
                w.write_u32::<LittleEndian>(v[0]).unwrap();
                w.write_u32::<LittleEndian>(v[1]).unwrap();
            }
        }
    }
}

fn load_kmer_feature(len: usize, data: &[u8]) -> Vec<IntMap<u32, [u32; 2]>> {
    let mut vec = vec![IntMap::default(); len];
    let (mut j, mut vec_i) = (0, 0);
    while j < data.len() {
        let map_len = LittleEndian::read_u16(&data[j..j + 2]) as usize;
        let map = &mut vec[vec_i];
        j += 2;
        vec_i += 1;
        map.reserve(map_len);
        for _ in 0..map_len {
            let k = LittleEndian::read_u32(&data[j..j + 4]);
            let v1 = LittleEndian::read_u32(&data[j + 4..j + 8]);
            let v2 = LittleEndian::read_u32(&data[j + 8..j + 12]);
            j += 12;
            map.insert(k, [v1, v2]);
        }
    }
    vec
}

pub fn output_kmer_features(features: &[Vec<IntMap<u32, [u32; 2]>>], errors: &[ErrorPattern]) {
    for (p, maps) in features.iter().enumerate() {
        for (kmer, (sym, count)) in maps
            .iter()
            .enumerate()
            .flat_map(|(pos, x)| x.iter().map(move |y| (pos, y)))
        {
            let error = &errors[*sym as usize];
            println!(
                "{}\t{}\t{}\t{}\t{}",
                kmer << COUNTER_BITS | p,
                &String::from_utf8_lossy(&error.seq_sym),
                error.ref_skip,
                count[0],
                count[1]
            );
        }
    }
}

fn dump_errors<W: Write>(w: &mut W, data: &ErrorPatternInterner) {
    let vec = &data.vec.read();
    w.write_u64::<LittleEndian>(vec.len() as u64).unwrap();
    for error in vec.iter() {
        w.write_u16::<LittleEndian>(
            error
                .ref_skip
                .try_into()
                .expect("Error ref_skip too large!"),
        )
        .unwrap();
        w.write_u16::<LittleEndian>(
            error
                .seq_sym
                .len()
                .try_into()
                .expect("Error seq_sym too large!"),
        )
        .unwrap();
    }
    for error in vec.iter() {
        w.write_all(&error.seq_sym).unwrap();
    }
}

fn load_errors<R: Read>(r: &mut R) -> Vec<ErrorPattern> {
    let mut buffer = [0; 8];
    r.read_exact(&mut buffer).unwrap();
    let count = LittleEndian::read_u64(&buffer) as usize;
    let mut meta = Vec::with_capacity(count);
    let mut total_seq_len = 0;

    for _ in 0..count {
        r.read_exact(&mut buffer[..4]).unwrap();
        let ref_skip = LittleEndian::read_u16(&buffer[..2]) as usize;
        let seq_len = LittleEndian::read_u16(&buffer[2..4]) as usize;
        total_seq_len += seq_len;
        meta.push((ref_skip, seq_len));
    }

    let mut buf = vec![0u8; total_seq_len];
    r.read_exact(&mut buf).unwrap();

    let mut errors = Vec::with_capacity(count);
    let mut offset = 0;
    for (ref_skip, seq_len) in meta {
        let end = offset + seq_len;
        errors.push(ErrorPattern::new(&buf[offset..end], ref_skip));
        offset = end;
    }

    errors
}

fn dump_error_transits<W: Write>(w: &mut W, data: &[[[u32; 5]; 5]]) {
    assert_eq!(data.len(), 22);
    for v in data.iter().flatten().flatten() {
        w.write_u32::<LittleEndian>(*v).unwrap();
    }
}

fn load_error_transits<R: Read>(r: &mut R) -> Vec<[[u32; 5]; 5]> {
    let mut buf = vec![0u8; 5 * 5 * 22 * 4];
    r.read_exact(&mut buf).unwrap();
    let mut error_transits = vec![[[0u32; 5]; 5]; 22];
    for (i, v) in error_transits.iter_mut().flatten().flatten().enumerate() {
        *v = LittleEndian::read_u32(&buf[i * 4..(i + 1) * 4]);
    }
    error_transits
}

fn dump_match_lens<W: Write>(w: &mut W, data: &[[[u32; 62]; 5]]) {
    assert_eq!(data.len(), 21);
    for v in data.iter().flatten().flatten() {
        w.write_u32::<LittleEndian>(*v).unwrap();
    }
}

fn load_match_lens<R: Read>(r: &mut R) -> Vec<[[u32; 62]; 5]> {
    let mut buf = vec![0u8; 62 * 5 * 21 * 4];
    r.read_exact(&mut buf).unwrap();
    let mut match_lens = vec![[[0u32; 62]; 5]; 21];
    for (i, v) in match_lens.iter_mut().flatten().flatten().enumerate() {
        *v = LittleEndian::read_u32(&buf[i * 4..(i + 1) * 4]);
    }
    match_lens
}

fn dump_histogram<W: Write>(w: &mut W, data: &EmpHistMode) {
    w.write_u32::<LittleEndian>(data.median).unwrap();
    w.write_u32::<LittleEndian>(data.edge.len() as u32).unwrap();
    for v in data.edge.iter() {
        w.write_u32::<LittleEndian>(*v).unwrap();
    }
    w.write_u32::<LittleEndian>(data.prob.len() as u32).unwrap();
    for v in data.prob.iter() {
        w.write_f64::<LittleEndian>(*v).unwrap();
    }
}

fn load_histogram<R: Read>(r: &mut R) -> EmpHistMode {
    let median = r.read_u32::<LittleEndian>().unwrap();
    let edge_len = r.read_u32::<LittleEndian>().unwrap() as usize;
    let mut edge = Vec::with_capacity(edge_len);
    let mut buf = vec![0u8; edge_len * 4];
    r.read_exact(&mut buf).unwrap();
    for i in 0..edge_len {
        let v = LittleEndian::read_u32(&buf[i * 4..(i + 1) * 4]);
        edge.push(v);
    }

    let prob_len = r.read_u32::<LittleEndian>().unwrap() as usize;
    let mut prob = Vec::with_capacity(prob_len);
    let mut buf = vec![0u8; prob_len * 8];
    r.read_exact(&mut buf).unwrap();
    for i in 0..prob_len {
        let v = LittleEndian::read_f64(&buf[i * 8..(i + 1) * 8]);
        prob.push(v);
    }
    EmpHistMode::new(prob, edge, median)
}

fn dump_chimeric_fra<W: Write>(w: &mut W, data: &[f64]) {
    w.write_u32::<LittleEndian>(data.len() as u32).unwrap();
    for v in data.iter() {
        w.write_f64::<LittleEndian>(*v).unwrap();
    }
}

fn load_chimeric_fra<R: Read>(r: &mut R) -> Vec<f64> {
    let len = r.read_u32::<LittleEndian>().unwrap() as usize;
    let mut data = Vec::with_capacity(len);
    let mut buf = vec![0u8; len * 8];
    r.read_exact(&mut buf).unwrap();
    for i in 0..len {
        let v = LittleEndian::read_f64(&buf[i * 8..(i + 1) * 8]);
        data.push(v);
    }
    data
}

fn dump_kmer_identitys<W: Write>(w: &mut W, data: &[(f32, f32)]) {
    w.write_u64::<LittleEndian>(data.len() as u64).unwrap();
    for (v1, v2) in data.iter() {
        w.write_f32::<LittleEndian>(*v1).unwrap();
        w.write_f32::<LittleEndian>(*v2).unwrap();
    }
}

fn load_kmer_identitys<R: Read>(r: &mut R) -> Vec<(f32, f32)> {
    let len = r.read_u64::<LittleEndian>().unwrap() as usize;
    let mut data = Vec::with_capacity(len);
    let mut buf = vec![0_u8; len * 8];
    r.read_exact(&mut buf).unwrap();
    for i in 0..len {
        let v1 = LittleEndian::read_f32(&buf[i * 8..(i * 8 + 4)]);
        let v2 = LittleEndian::read_f32(&buf[(i * 8 + 4)..(i + 1) * 8]);
        data.push((v1, v2));
    }
    data
}

fn dump_read_type<W: Write>(w: &mut W, data: &str) {
    w.write_u64::<LittleEndian>(data.len() as u64).unwrap();
    w.write_all(data.as_bytes()).unwrap();
}

fn load_read_type<R: Read>(r: &mut R) -> String {
    let len = r.read_u64::<LittleEndian>().unwrap();
    let mut buf = vec![0u8; len as usize];
    r.read_exact(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

pub fn dump_file(
    read_type: &str,
    ksize: usize,
    thread: u32,
    strand_fra: f64,
    chimeric_fra: &[f64],
    kmer_features: &[Vec<IntMap<u32, [u32; 2]>>],
    errors: &ErrorPatternInterner,
    error_transits: &[[[u32; 5]; 5]],
    match_lens: &[[[u32; 62]; 5]],
    emp_aligned_len: &EmpHistMode,
    emp_unaligned_left_len: &EmpHistMode,
    emp_unaligned_right_len: &EmpHistMode,
    kmer_identitys: &[(f32, f32)],
) {
    let mut encoder = Encoder::new(BufWriter::new(io::stdout().lock()), 19)
        .expect("Failed to create zstd encoder");
    encoder
        .multithread(thread)
        .expect("Failed to set multithreading");
    let mut w = encoder.auto_finish();

    w.write_all(CFSIM_MAGIC).unwrap();
    w.write_u64::<LittleEndian>(ksize as u64).unwrap();
    w.write_f64::<LittleEndian>(strand_fra).unwrap();
    dump_read_type(&mut w, read_type);
    dump_kmer_features(&mut w, kmer_features);
    dump_errors(&mut w, errors);
    dump_error_transits(&mut w, error_transits);
    dump_match_lens(&mut w, match_lens);
    dump_histogram(&mut w, emp_aligned_len);
    dump_histogram(&mut w, emp_unaligned_left_len);
    dump_histogram(&mut w, emp_unaligned_right_len);
    dump_chimeric_fra(&mut w, chimeric_fra);
    dump_kmer_identitys(&mut w, kmer_identitys);
    w.flush().unwrap();
}

pub fn load_file(
    path: &str,
    thread: usize,
) -> (
    String,
    usize,
    f64,
    Vec<f64>,
    [Vec<IntMap<u32, [u32; 2]>>; 1 << COUNTER_BITS],
    Vec<ErrorPattern>,
    Vec<[[u32; 5]; 5]>,
    Vec<[[u32; 62]; 5]>,
    EmpHistMode,
    EmpHistMode,
    EmpHistMode,
    Vec<(f32, f32)>,
) {
    let mut r = Decoder::new(FileIO::new(path).reader()).expect("Failed to create zstd decoder");
    let mut buffer = [0; 21];
    r.read_exact(&mut buffer).unwrap();
    assert!(
        &buffer[..5] == CFSIM_MAGIC,
        "The input binary mode file is incompatible."
    );
    let ksize = LittleEndian::read_u64(&buffer[5..13]);
    let strand_fra = LittleEndian::read_f64(&buffer[13..21]);
    let read_type = load_read_type(&mut r);
    let kmer_features_len = r.read_u64::<LittleEndian>().unwrap() as usize;

    let (
        kmer_features,
        (
            errors,
            error_transits,
            match_lens,
            emp_aligned_len,
            emp_unaligned_left_len,
            emp_unaligned_right_len,
            chimeric_fra,
            kmer_identitys,
        ),
    ) = thread::scope(|work| {
        let (in_s, in_r) = bounded(thread);
        let (ou_s, ou_r) = bounded(thread);
        let others = work.spawn(move || {
            for i in 0..kmer_features_len {
                let mut buffer = [0; 16];
                r.read_exact(&mut buffer).unwrap();
                let vec_len = LittleEndian::read_u64(&buffer[..8]) as usize;
                let map_len = LittleEndian::read_u64(&buffer[8..]) as usize;
                let mut buffer = vec![0; vec_len * 2 + map_len * 12];
                r.read_exact(&mut buffer).unwrap();
                in_s.send((i, vec_len, buffer)).unwrap();
            }

            let errors = load_errors(&mut r);
            let error_transits = load_error_transits(&mut r);
            let match_lens = load_match_lens(&mut r);
            let emp_reads_len = load_histogram(&mut r);
            let emp_unaligned_left_len = load_histogram(&mut r);
            let emp_unaligned_right_len = load_histogram(&mut r);
            let chimeric_fra = load_chimeric_fra(&mut r);
            let kmer_identitys = load_kmer_identitys(&mut r);
            (
                errors,
                error_transits,
                match_lens,
                emp_reads_len,
                emp_unaligned_left_len,
                emp_unaligned_right_len,
                chimeric_fra,
                kmer_identitys,
            )
        });

        (0..thread).for_each(|_| {
            let in_r = in_r.clone();
            let ou_s = ou_s.clone();
            work.spawn(move || {
                while let Ok((i, len, data)) = in_r.recv() {
                    let vec = load_kmer_feature(len, &data);
                    ou_s.send((i, vec)).unwrap();
                }
            });
        });
        drop(ou_s);

        let kmer_features = work.spawn(move || {
            let mut kmer_features: [Vec<IntMap<u32, [u32; 2]>>; 1 << COUNTER_BITS] =
                std::array::from_fn(|_| Vec::new());
            while let Ok((i, vec)) = ou_r.recv() {
                kmer_features[i] = vec;
            }
            kmer_features
        });
        (kmer_features.join().unwrap(), others.join().unwrap())
    });
    (
        read_type,
        ksize as usize,
        strand_fra,
        chimeric_fra,
        kmer_features,
        errors,
        error_transits,
        match_lens,
        emp_aligned_len,
        emp_unaligned_left_len,
        emp_unaligned_right_len,
        kmer_identitys,
    )
}
