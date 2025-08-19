#![allow(dead_code)]

use rust_htslib::bam::{
    Record,
    header::{Header, HeaderRecord},
    record::{Cigar, CigarString},
};
use std::fmt::Write;

const READS_CAP: usize = 1000;

pub struct SimReads {
    pub reads: Vec<Record>,
}

impl SimReads {
    pub fn new() -> Self {
        Self {
            reads: Vec::with_capacity(READS_CAP),
        }
    }

    pub fn push(
        &mut self,
        seq: &[u8],
        cigar_chars: &[u8],
        qual: &mut Vec<u8>,
        name: &str,
        start: usize,
        end: usize,
        unaligned_left_len: u32,
        aligned_len: u32,
        unaligned_right_len: u32,
        tid: i32,
        is_positive: bool,
        chimeric: usize,
        thread: u64,
        index: usize,
        qname_buffer: &mut String,
    ) {
        let flag: u16 = if is_positive { 0 } else { 0x10 };
        let qname = {
            qname_buffer.clear();
            write!(qname_buffer, "{:03}_{index}_{name}_{start}_{end}_{unaligned_left_len}_{aligned_len}_{unaligned_right_len}_{}_{chimeric}", thread, if is_positive {0} else {1}).unwrap();
            qname_buffer.as_bytes()
        };

        let cigars = build_cigar_ops(cigar_chars);
        if qual.is_empty() {
            qual.resize(seq.len(), 255);
        }

        let mut record = Record::new();
        record.set(qname, Some(&CigarString(cigars)), seq, qual);
        record.set_tid(tid);
        record.set_pos(start as i64);
        record.set_mapq(60);
        record.set_flags(flag);

        self.reads.push(record);
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

fn build_cigar_ops(cigar_chars: &[u8]) -> Vec<Cigar> {
    let mut cigars = Vec::with_capacity(1024);
    let iter = cigar_chars.iter();

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
    cigars.shrink_to_fit();
    cigars
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

    header
}
