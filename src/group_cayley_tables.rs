use std::collections::HashSet;

use crate::{hs::HyperGroupoid, hypergroups::HyperGroup};
pub const S3:[u64;36] = [
        0, 1, 2, 3, 4, 5,
        1, 0, 5, 4, 3, 2, 
        2, 4, 0, 5, 1, 3,
        3, 5, 4, 0, 2, 1,
        4, 2, 3, 1, 5, 0,
        5, 3, 1, 2, 0, 4
    ];
