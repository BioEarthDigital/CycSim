#![allow(dead_code)]

use rustc_hash::FxHashMap as HashMap;
use rust_htslib::bam::{
    Record,
    header::{Header, HeaderRecord},
    record::{Aux, Cigar, CigarString},
};
use std::{env, fmt::Write, path::Path};

const READS_CAP: usize = 1000;

type Read<'a> = (
    &'a [u8],     // seq
    &'a [u8],     // cigar_chars
    &'a [u8],     // qual
    &'a str,      // name
    usize,        // start
    usize,        // end
    u32,          // unaligned_left_len
    u32,          // aligned_len
    u32,          // unaligned_right_len
    i32,          // tid
    bool,         // is_positive
    usize,        // chimeric
    u64,          // thread
    usize,        // index
);

type ChiRead<'a> = (
    usize,        // seq_start
    usize,        // cigar_start
    &'a str,      // name
    usize,        // start
    usize,        // end
    u32,          // unaligned_left_len
    u32,          // aligned_len
    u32,          // unaligned_right_len
    i32,          // tid
    bool,         // is_positive
    usize,        // chimeric
    u64,          // thread
    usize,        // index
);

pub struct SimReads {
    pub reads: Vec<Record>,
}

impl SimReads {
    pub fn new() -> Self {
        Self {
            reads: Vec::with_capacity(READS_CAP),
        }
    }

    pub fn push(&mut self, args: Read, qname_buffer: &mut String) {
        let (seq, cigar_chars, qual, name, start, end, unaligned_left_len, 
            aligned_len, unaligned_right_len, tid, is_positive, chimeric, 
            thread, index) = args;

        let flag: u16 = if is_positive { 0 } else { 0x10 };

        let qname = {
            qname_buffer.clear();
            write!(qname_buffer, "{:03}_{index}_{name}_{start}_{end}_{unaligned_left_len}_{aligned_len}_{unaligned_right_len}_{}_{chimeric}", 
                thread, if is_positive {0} else {1}).unwrap();
            qname_buffer.as_bytes()
        };

        let cigars = build_cigar_ops(cigar_chars, 0, 0);

        let mut record = Record::new();
        record.set(qname, Some(&CigarString(cigars)), seq, qual);
        record.set_tid(tid);
        record.set_pos(start as i64);
        record.set_mapq(60);
        record.set_flags(flag);

        self.reads.push(record);
    }

    pub fn push_chimeric(&mut self, args: &[ChiRead], ref_names: &[(String, Vec<u8>)], 
        all_seq: &mut [u8], all_cigar: &mut [u8], all_qual: &mut [u8], qname_buffer: &mut String)
    {
        qname_buffer.clear();
        let mut primary_index = 0;
        for (i, arg) in args.iter().enumerate() {
            let (_, _, name, start, end, unaligned_left_len, aligned_len, unaligned_right_len,
                 _, is_positive, chimeric, thread, index) = arg;
            if i > 0 {
                qname_buffer.push('|');
            }
            write!(qname_buffer, 
                "{:03}_{index}_{name}_{start}_{end}_{unaligned_left_len}_{aligned_len}_{unaligned_right_len}_{}_{chimeric}", 
                thread, if *is_positive {0} else {1}
            ).unwrap();

            if *aligned_len > args[primary_index].6 {
                primary_index = i;
            }
        }

        //htslib qname < 254
        let qname_max_len = qname_buffer.len().min(254);
        let qname = &qname_buffer.as_bytes()[..qname_max_len];

        let mut sa_entries: Vec<String> = Vec::with_capacity(args.len());
        for (i, arg) in args.iter().enumerate() {
            let (seq_start, cigar_start, _, start, _, _, _, _,
                 tid, is_positive, _, _, _) = arg;
            let (seq_end, cigar_end) = if i == args.len() - 1 { 
                (all_seq.len(), all_cigar.len())
            }else {
                let next_args = &args[i+1];
                (next_args.0, next_args.1)
            };//exclude

            let (start_clip, end_clip) = if *is_positive {
                (*seq_start, all_seq.len() - seq_end)
            }else {
                (all_seq.len() - seq_end, *seq_start)
            };

            if !*is_positive {
                rev_comp_inplace(&mut all_seq[*seq_start..seq_end]);
                all_qual[*seq_start..seq_end].reverse();
            }

            let (cigar_str, nm) = build_sa_cigar_ops(&all_cigar[*cigar_start..cigar_end], start_clip, end_clip);
            let strand = if *is_positive { '+' } else { '-' };
            
            // format: rname,pos,strand,CIGAR,mapQ,NM;
            let entry = format!("{},{},{strand},{cigar_str},60,{nm};", 
                ref_names[*tid as usize].0, 
                start + 1);
            sa_entries.push(entry);
        }

        if !args[primary_index].9 {
            rev_comp_inplace(all_seq);
            all_qual.reverse();
        }

        for (i, arg) in args.iter().enumerate() {
            let (seq_start, cigar_start, _, start, _, _ul, _, _,
                 tid, is_positive, _, _, _) = arg;
            let (seq_end, cigar_end) = if i == args.len() - 1 { 
                (all_seq.len(), all_cigar.len())
            }else {
                let next_args = &args[i+1];
                (next_args.0, next_args.1)
            };//exclude

            let (start_clip, end_clip) = if *is_positive {
                (*seq_start, all_seq.len() - seq_end)
            }else {
                (all_seq.len() - seq_end, *seq_start)
            };
            let cigars = build_cigar_ops(&all_cigar[*cigar_start..cigar_end], start_clip, end_clip);

            let flag: u16 = if *is_positive { 0 } else { 0x10 };
            let flag = if i == primary_index { flag } else { flag | 0x800 };
            
            let mut record = Record::new();
            record.set(qname, Some(&CigarString(cigars)), all_seq, all_qual);
            record.set_tid(*tid);
            record.set_pos(*start as i64);
            record.set_mapq(60);
            record.set_flags(flag);

            let mut sa_tag = String::with_capacity(1024);
            for (j, entry) in sa_entries.iter().enumerate() {
                if j != i {
                    sa_tag.push_str(entry);
                }
            }
            record.push_aux(b"SA", Aux::String(&sa_tag)).unwrap();

            self.reads.push(record);
        }
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    pub fn clear(&mut self) {
        self.reads.clear();
    }

    pub fn is_full(&self) -> bool {
        self.reads.len() >= READS_CAP
    }

    pub fn iter(&self) -> impl Iterator<Item = &Record> {
        self.reads.iter()
    }
}

fn rev_comp_inplace(seq: &mut [u8]) {
    const COMP: [u8; 256] = {
        let mut table = [0u8; 256];
        let mut i = 0;
        while i < 256 {
            table[i] = i as u8;
            i += 1;
        }

        table[b'A' as usize] = b'T';
        table[b'T' as usize] = b'A';
        table[b'C' as usize] = b'G';
        table[b'G' as usize] = b'C';
        table[b'a' as usize] = b't';
        table[b't' as usize] = b'a';
        table[b'c' as usize] = b'g';
        table[b'g' as usize] = b'c';
        table
    };

    let len = seq.len();
    for i in 0..(len / 2) {
        let a = seq[i];
        let b = seq[len - 1 - i];
        seq[i] = COMP[b as usize];
        seq[len - 1 - i] = COMP[a as usize];
    }

    if len % 2 == 1 {
        seq[len / 2] = COMP[seq[len / 2] as usize];
    }
}

pub fn reverse_complement_seq(seq: &[u8]) -> Vec<u8> {
    let mut seq = seq.to_owned();
    rev_comp_inplace(&mut seq);
    seq
}

fn to_cigar(op: u8, len: u32) -> Cigar {
    match op {
        b'M' | b'=' => Cigar::Equal(len),
        b'X' => Cigar::Diff(len),
        b'I' => Cigar::Ins(len),
        b'D' => Cigar::Del(len),
        b'S' => Cigar::SoftClip(len),
        b'H' => Cigar::HardClip(len),
        b'N' => Cigar::RefSkip(len),
        _ => panic!("Unsupported CIGAR op: {}", op as char),
    }
}

fn build_cigar_ops(cigar_chars: &[u8], start_clip: usize, end_clip: usize) -> Vec<Cigar> {
    let mut cigars = Vec::with_capacity(1024);

    let mut start = 0;
    while cigar_chars[start] == b'S' {
        start += 1;
    }

    let n = cigar_chars.len();
    let mut end = n;
    while cigar_chars[end - 1] == b'S' {
        end -= 1;
    }

    if start + start_clip > 0 {
        cigars.push(to_cigar(b'S', (start + start_clip) as u32));
    }

    let iter = cigar_chars[start..end].iter();

    let mut iter = iter.peekable();
    while let Some(&op) = iter.next() {
        let mut count = 1;
        while let Some(&&next_op) = iter.peek() {
            if next_op == op {
                iter.next();
                count += 1;
            } else {
                break;
            }
        }
        cigars.push(to_cigar(op, count));
    }

    if end_clip + n - end > 0 {
        cigars.push(to_cigar(b'S', (end_clip + n - end) as u32));
    }

    cigars.shrink_to_fit();
    cigars
}

fn build_sa_cigar_ops(cigar_chars: &[u8], start_clip: usize, end_clip: usize) -> (String, usize) {
    let mut start = 0;
    while cigar_chars[start] == b'S' {
        start += 1;
    }
    let n = cigar_chars.len();
    let mut end = n;
    while cigar_chars[end - 1] == b'S' {
        end -= 1;
    }

    let mut counts: HashMap<u8, usize> = HashMap::default();
    let mut order: Vec<u8> = Vec::with_capacity(8);

    for &op in &cigar_chars[start..end] {
        if !counts.contains_key(&op) {
            order.push(op);
        }
        *counts.entry(op).or_insert(0) += 1;
    }

    let mut nm = 0;
    let mut result = String::with_capacity(64);

    if start + start_clip > 0 {
        write!(result, "{}S", start + start_clip).unwrap();
    }

    for &op in &order {
        let count = counts[&op];
        if op != b'='{
            nm += count;
        }
        write!(result, "{count}{}", op as char).unwrap();
    }

    if end_clip + n - end > 0 {
        write!(result, "{}S", end_clip + n - end).unwrap();
    }

    (result, nm)
}

pub fn create_bam_header(seqs: &[(String, Vec<u8>)]) -> Header {
    let mut header = Header::new();
    let mut hd_record = HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", "1.6");
    hd_record.push_tag(b"SO", "unsorted");
    header.push_record(&hd_record);

    for (name, seq) in seqs {
        let mut sq_record = HeaderRecord::new(b"SQ");
        sq_record.push_tag(b"SN", name);
        sq_record.push_tag(b"LN", seq.len());
        header.push_record(&sq_record);
    }

    add_pg_record(header)
}

pub fn add_pg_record(mut header: Header) -> Header{
    let mut args: Vec<String> = env::args().collect();
    let prog_name = args
        .get(0)
        .and_then(|s| Path::new(s).file_name().map(|os| os.to_string_lossy().to_string()))
        .unwrap_or_else(|| "unknown".to_string());

    if let Some(first) = args.first_mut() {
        *first = prog_name.clone(); // 替换成 basename
    }

    let cmdline = args.join(" ");

    let mut pg_record = HeaderRecord::new(b"PG");
    pg_record.push_tag(b"ID", &prog_name);
    pg_record.push_tag(b"PN", &prog_name);
    pg_record.push_tag(b"VN", env!("VERSION"));
    pg_record.push_tag(b"CL", &cmdline);
    header.push_record(&pg_record);
    header
}