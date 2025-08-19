#![allow(dead_code)]
const SEQ_NUM: [u8; 128] = [
    // translate ACGTU-NM to 01233456
    //A, C,  G,  T,   -,  N, M
    65, 67, 71, 84, 45, 78, 77, 4, 4, 4, 4, 4, 4, 4, 4, 4, //0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 6, 5, 4, //64-79
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //80-95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 6, 5, 4, //96-111
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //112-127
];

type Ktype = u32;
pub fn iter2kmer(mut iter: impl Iterator<Item = u8>, ksize: usize) -> impl Iterator<Item = Ktype> {
    assert!(
        ksize <= std::mem::size_of::<Ktype>() * 4,
        "ksize should be <= {}",
        std::mem::size_of::<Ktype>() * 4
    );
    let mut l = 0;
    // let shift = 2 * (ksize - 1);
    let mut mask = (1 << (2 * ksize)) - 1;
    //for ksize == 32 if T is u64 or ksize == 64 if T is u128
    if mask == 0 {
        mask -= 1;
    }
    let mut kmer = 0;
    std::iter::from_fn(move || {
        loop {
            if let Some(c) = iter.next() {
                let c = SEQ_NUM[c as usize] as Ktype;
                if c < 4 {
                    kmer = (kmer << 2 | c) & mask; // forward k-mer
                    l += 1;
                } else {
                    l = 0;
                    return None; //break N base
                }

                if l >= ksize {
                    return Some(kmer);
                }
            } else {
                return None;
            }
        }
    })
}

pub fn seq2kmer(seq: &[u8]) -> Ktype {
    seq.iter()
        .fold(0u32, |acc, &b| (acc << 2) | SEQ_NUM[b as usize] as Ktype)
}

pub fn kmer2seq(mut kmer: Ktype, ksize: usize) -> Vec<u8> {
    let mut seq = vec![0u8; ksize];
    for i in (0..ksize).rev() {
        let base = (kmer & 0b11) as usize;
        seq[i] = SEQ_NUM[base];
        kmer >>= 2;
    }
    seq
}
