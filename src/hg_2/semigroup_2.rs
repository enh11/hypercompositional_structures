// This are tags for representants of semigroups of order 2 up to  isomorpphism
pub const TAGS_REPRESENTANTS_SEMIGROUPS_2:[u128;5]=[90,102,85,86,105];
// This are the Cayley table of representants of semigroups of order 2 up to 
// isomorphism.  This can be obtained using a GAP library "smallsemi" and are available up to order 8.
// S2_2 and S2_3 are anti isomorphic.

pub const S2_0: [u64; 4] = [
    0,     0,    
    0,     0, ];
pub const S2_1: [u64; 4] = [
    0,     0,     
    0,     1, ];
pub const S2_2: [u64; 4] = [
    0,     0,     
    1,     1, ];
pub const S2_3: [u64; 4] = [
    0,     1,     
    1,     0, ];

pub const SEMIGROUP_2: [&[u64]; 4] = [
    &S2_0, &S2_1, &S2_2,&S2_3
];
