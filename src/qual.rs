use std::sync::LazyLock;

pub const MAX_QUAL: u8 = 40;
pub const MIN_QUAL: u8 = 1;

const TABLE_SIZE: usize = 100_001;

static IDENTITY_TO_Q_TABLE: LazyLock<Vec<u8>> = LazyLock::new(|| {
    let mut table = vec![0u8; TABLE_SIZE];
    for (i, v) in table.iter_mut().enumerate().take(TABLE_SIZE) {
        let identity = i as f64 / (TABLE_SIZE as f64 - 1.0);
        let err = 1.0 - identity;
        let q_val = -10.0 * err.log10();
        let q = q_val.round() as u8;
        *v = q.clamp(MIN_QUAL, MAX_QUAL);
    }
    table
});

// should + 33 to ASCII char
pub fn identity_to_qual(identity: f64) -> u8 {
    let idx = (identity.min(1.0) * (TABLE_SIZE as f64 - 1.0)).round() as usize;
    IDENTITY_TO_Q_TABLE[idx]
}
