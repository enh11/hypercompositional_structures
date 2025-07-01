
use std::{collections::HashSet, fmt::Binary, fs::File, io::Write};
use std::fmt::Debug;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;

use uint::construct_uint;
construct_uint!{
    pub struct U1024(32);
}
/*THIS IS A WORK in PROGRESS*/
#[derive(Clone, Debug)]
pub struct U1024Range {
    start: U1024,
    end: U1024,
}

impl Iterator for U1024Range {
    type Item = U1024;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let current = self.start;
            self.start += U1024::from(1u8);
            Some(current)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, None) // length might not fit usize
    }
}

impl DoubleEndedIterator for U1024Range {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            self.end -= U1024::from(1u8);
            Some(self.end)
        } else {
            None
        }
    }
}

pub trait U1024RangeExt {
    fn to(self, end: U1024) -> U1024Range;
}

impl U1024RangeExt for U1024 {
    fn to(self, end: U1024) -> U1024Range {
        U1024Range { start: self, end }
    }
}
/*MUST BE IMPROVED ALSO WITH PARALLEL ITERATOR */

impl U1024 {
pub fn u1024_to_binary_vec(&self)->Vec<Vec<u64>>{
    self.0.iter()
        .map(|x| u64_to_binary_vec_u64(x, &64))
        .rev()
        .collect_vec()
}
pub fn from_binary_vec(binary_vec:&Vec<u64>)->Self{
    let length = binary_vec.len();
    (0..length).into_iter()
        .fold(U1024::zero(), 
            |acc,x|
            acc+(U1024::from(binary_vec[x])<<length-1-x))
}
}

/// Converts a Vec<u64> into a HashSet<u64>
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// let v=vec![1,2,0];
/// let v_as_set=vec_to_set(&v);
/// println!("{:?}",v_as_set);
///

pub fn vec_to_set(v:&Vec<u64>)->HashSet<u64>{
    v.iter().map(|x|*x).collect()
}
pub fn power_set (n:&u64)->Vec<Vec<u64>>{
    (0..2u64.pow(*n as u32)).into_iter().map(|i|  get_subset(&i, &n)).collect()
  }
/// Represents the integer $k$ as a subset of the set H={0,1,..,n-1}.
/// There are 2^n different integer representing subsets of H. It will panic if $k is greater then 2^n.
/// The subset S represented by k is given by the binary representation of k: 
/// i is in S if and only if the i-th digit of binary representation of k is a one.
/// Output is Vec<u64>. Use fn vec_to_set to convert it into a HashSet<u64>.
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// use hyperstruc::utilities::get_subset;
/// use std::collections::HashSet;
/// 
/// let cardinality = 4u64;
/// let k=6;
/// let subset_as_vec=get_subset(&k,&cardinality);
/// let subset= vec_to_set(&subset_as_vec);
/// 
/// println!("{:?}",subset);
/// let test_subset:HashSet<u64>= (1..=2).into_iter().collect();
/// assert_eq!(subset,test_subset);
///
/// let cardinality = 3u64;
/// let k=2;
/// let subset_as_vec=get_subset(&k,&cardinality);
/// let subset= vec_to_set(&subset_as_vec);
/// 
/// println!("{:?}",subset);
/// let test_subset:HashSet<u64>= [1].into_iter().collect();
/// assert_eq!(subset,test_subset);
pub fn get_subset(k:&u64,cardinality:&u64)->Vec<u64>{
    let mut subset: Vec<u64> = Vec::new();
    if k>=&2u64.pow(*cardinality as u32){panic!("k can't be grater then 2^n");}
    for i in 0..*cardinality {
        if (k >> i)&1==1{
            subset.push(i as u64);
            }
    }
    subset.iter().map(|x|*x).collect()
}
///
/// Let `k` be the integer which identifies a subset `A` of `H`. This function compute the integer k^c which correspond to the complement of k with respect to H.
/// For example, for `A = {0,2}` we have `k = 5`, whose binary representation is `101`. Then `A^c  = {1}`, which is represented by the integer `2`.
/// 
/// # Example
/// ```
/// use hyperstruc::utilities::get_complement_subset;
/// 
/// let cardinality = 3u64;
/// let subset_representation = 5u64;
/// let expected_complement_representation = 2u64;
/// let complement = get_complement_subset(&subset_representation, &cardinality);
/// 
/// assert_eq!(expected_complement_representation, complement)
///
pub fn get_complement_subset(k:&u64,cardinality:&u64)->u64 {
    // (1<<*cardinality)-1-*k This works as well
    // !k%(1<<cardinality) This works as well 
    k^((1<<cardinality)-1) //This works as well and maybe is more efficient
}

/// Returns the positions of `1`s in the binary representation of a subset indicator.
///
/// This function computes the *characteristic function* (as positions of `1`s) 
/// for a subset represented as an integer `subset_a` within the powerset of a set 
/// of size `cardinality`. It returns the indices of elements that are present in the subset.
/// 
/// # Arguments
/// - `subset_a`: A `u64` value representing a subset of the base set (as a bitmask).
/// - `cardinality`: The size `n` of the base set; the valid range of `subset_a` is `[0, 2^n - 1]`.
/// 
/// # Returns
/// A `Vec<usize>` containing the zero-based indices corresponding to set bits in `subset_a`.
///
/// # Panics
/// Panics if `subset_a >= 2^cardinality`, as it would not represent a valid subset.
///
/// # Example
/// ```
/// use hyperstruc::utilities::support;
///
/// let cardinality = 3u64;
/// let subset_representation = 5u64; // binary: 101
/// let expected = vec![0, 2];        // positions of 1s
///
/// let result = support(&subset_representation, &cardinality);
/// assert_eq!(result, expected);
/// ```
pub fn support(subset_a: &u64, cardinality: &u64) -> Vec<usize> {
    let power_set_cardinality = 1 << cardinality;
    if subset_a >= &power_set_cardinality {
        panic!(
            "{} does not represent a subset of H.\nsubset_a must be between 0 and {}",
            subset_a,
            power_set_cardinality - 1
        );
    }

    (0..*cardinality as usize)
        .filter(|x| (subset_a >> x) & 1 == 1)
        .collect()
}
pub fn cartesian_product(set: &Vec<u64>) -> Vec<(u64, u64)> {
    let mut product: Vec<(u64, u64)> = Vec::new();
    for &a in set {
        for &b in set {
            product.push((a, b));
        }
    }
    product
}
/// Represents a subset $S$ as an integer in H={0,1,..,n-1}.
/// There are 2^n different subsets of H. It will panic if $k is greater then 2^n.
/// Elements in S correspond to indices of occurrences of ones for the binary representation of the integer k which identifies S.
/// Therefore, k equals the sum of power two with exponents in S. 
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// use hyperstruc::utilities::subset_as_u64;
/// use hyperstruc::utilities::get_subset;
/// use std::collections::HashSet;
/// 
/// let subset:HashSet<u64>= (1..=2).into_iter().collect();
/// let subset_as_integer= subset_as_u64(&subset);
/// let test_integer=6u64;
/// assert_eq!(subset_as_integer,test_integer);
/// 
/// let cardinality=4;
/// let k=8;
/// let subset=vec_to_set(&get_subset(&k,&cardinality));
/// let subset_as_integer= subset_as_u64(&subset);
/// assert_eq!(subset_as_integer,k);
/// 
///
pub fn subset_as_u64(k:&HashSet<u64>)->u64{
    k.iter().map(|x|2u64.pow(*x as u32)).sum()
}

pub fn permutaton_matrix_from_permutation(n:&u64,sigma:&Permutation)->DMatrix<u64>{
    let identity: DMatrix<u64>=DMatrix::identity(*n as usize,*n as usize);
    let rows:Vec<Vec<u64>> = identity.row_iter().map(|x|x.iter().map(|z|*z).collect()).collect();

    let x: Vec<u64> =sigma.clone().inverse().apply_slice(rows).concat();
    DMatrix::from_row_slice(*n as usize, *n as usize, &x).transpose()
}
/// Gets binary representation of an integer k with width-numbers-of-bit.
/// # Example
/// ```
/// use hyperstruc::utilities::n_to_binary_vec;/// 
/// let cardinality=4;
/// let k=6;
/// let binary=n_to_binary_vec(&k,&cardinality);
/// let v=vec![0,1,1,0];
/// assert_eq!(v,binary);
///
pub fn n_to_binary_vec(k: &u128, cardinality:&u64) -> Vec<u64> {
	format!("{:0>cardinality$}", format!("{:b}", k), cardinality = *cardinality as usize)
		.chars()
		.map(|x| if x == '1' { 1u64 } else { 0u64 })
		.collect()
}
pub fn u64_to_binary_vec_u64(k: &u64, width: &u32) -> Vec<u64> {
	format!("{:0>width$}", format!("{:b}", k), width = *width as usize)
		.chars()
		.map(|x| if x == '1' { 1u64 } else { 0u64 })
		.collect()
}
pub fn u64_to_set(subset_a:&u64,cardinality:&u64)->HashSet<u64>{
    vec_to_set(&get_subset(subset_a, cardinality))
}
/// Compute U1024 given its vector of 32-binary representation vectors. 
/// 
/// 
/// # Example
/// ```
/// 
/// use hyperstruc::utilities::U1024;
/// use hyperstruc::utilities::binary_to_u1024;
/// let binary_vector:Vec<Vec<u64>>=[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0]].to_vec().iter().map(|x|x.to_vec()).collect();
/// let n = U1024::from_dec_str(&"1234543446698765423435676545654332").unwrap();
/// let x:U1024=binary_to_u1024(&binary_vector);
/// assert_eq!(n,x);
/// 
///
pub fn binary_to_u1024(binary_vec:&Vec<Vec<u64>>)->U1024 {
    let mut bin_to_ret = [0u64;32];
    for i in 0..32{
        bin_to_ret[i]=binary_to_n(&binary_vec[31-i])   
    };
    U1024(bin_to_ret)
} 
impl Binary for U1024 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.u1024_to_binary_vec().concat().iter().map(|x|x.to_string()).collect_vec().concat();
write!(f,"{}",s)  }
}

/// Compute u64 given its binary representation vector.
/// 
/// 
/// # Example
/// ```
/// 
/// use hyperstruc::utilities::binary_to_n;
/// 
/// let binary_vector=[1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0].to_vec();
/// let n :u64=14987979559888969728;
/// let x:u64=binary_to_n(&binary_vector);
/// assert_eq!(n,x);
/// 
/// 
pub fn binary_to_n<T: num_traits::Num>(binary_vec:&Vec<u64>)->T where T: num_traits::Num<FromStrRadixErr: Debug> {
    let s= binary_vec.iter().map(|x|x.to_string()).collect::<String>();
    T::from_str_radix(&s, 2).unwrap()

} 
/// The input value is k in (0..2^n), therefore it represent a subset of H with |H|=n.
/// Any occurrence of 1 in the binary representation of k correspond to an element in the subset S corresponding to k 
/// The input value sigma is a permutation of S_n. We build the corresponding permutation matrix and we make it act on the binary representation of k.   
/// We prefer to inverse and normalize sigma.
/// For example for n=3 and k = 6 = '110', the subset corresponding to k is S_k= {1,2}. Let sigma = (1,2,0), meaning 0->1,1->2,2->0.
/// Then, for k=101, sigmra send the 0th bit into the 1st position, the 1st bit into the 2nd position and the 2nd bit into the 0th position.
/// Therefore, sigma_k=101, i.e., the subset {0,2}, which corresponds to the integer 5.
/// 
/// /// # Example
/// ```
/// use hyperstruc::utilities::representation_permutation_subset;
/// use permutation::Permutation;
/// 
/// let sigma= Permutation::oneline([1,2,0].to_vec());
/// let k =6u128;
/// let sigma_k=representation_permutation_subset(&k, &sigma);
/// let expected_subset_representation =5u128;
/// 
/// assert_eq!(expected_subset_representation,sigma_k)
/// 
/// 
pub fn representation_permutation_subset (k:&u128,sigma:&Permutation)->u128 {

    let binary_k=n_to_binary_vec(&k,&(sigma.len() as u64)).iter().rev().map(|x|*x).collect_vec();
    let x =sigma.apply_slice(&binary_k).iter().rev().map(|x|*x).collect_vec();
    
    binary_to_n(&x)
}
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::utilities::U1024;
/// use hyperstruc::utilities::representing_hypergroupoid_u1024;
/// 
/// let card = 3u64;
/// let tag = U1024::from(22150143u128);
/// 
/// assert!(representing_hypergroupoid_u1024(&tag,&card))
/// 
/// 
pub fn representing_hypergroupoid_u1024(tag:&U1024,cardinality:&u64)->bool{
    let cardinality_u1024=U1024::from(*cardinality);
    let (min,max) = get_min_max_u1024(cardinality);
    let mut clone_tag = tag.clone();
    let power_set = U1024::from(2).pow(cardinality_u1024);
    if tag<&min||tag>&max{return false;}
    while clone_tag>=power_set {
        if clone_tag.trailing_zeros()>=*cardinality as u32 {
            
            return false;
        }
        clone_tag>>=*cardinality;
    }
    true
}
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::utilities::representing_hypergroupoid;
/// 
/// let card = 3u64;
/// let tag = 22150143u128;
/// 
/// assert!(representing_hypergroupoid(&tag,&card))
/// 
/// 
pub fn representing_hypergroupoid(tag:&u128,cardinality:&u64)->bool{
    let (min,max) = get_min_max(cardinality);
    let mut clone_tag = tag.clone();
    let power_set = 2u128.pow(*cardinality as u32);
    if tag<&min||tag>&max{return false;}

    while clone_tag>=power_set {
        if clone_tag.trailing_zeros()>=*cardinality as u32 {
            
            return false;
        }
        clone_tag>>=cardinality;
    }
    true
}
pub fn collect_n_digits(width: &u64,m_representation_hypergroupoid:&u128)->Vec<u64>{
    let binary_n = n_to_binary_vec(m_representation_hypergroupoid, &(*width as u64));
    let out :Vec<u64>=(0..*width)
        .rev()
        .into_iter()
        .map(|i|
            binary_n.iter().rev().nth(i as usize).unwrap()).
            into_iter()
            .map(
                |x| if *x == 1 { 1u64 } 
                    else { 0u64 })
            .collect();
    out
}
pub fn collect_n_digits_u1024(width: &u64,m_representation_hypergroupoid:&U1024)->Vec<u64>{
    let binary_n = m_representation_hypergroupoid.u1024_to_binary_vec().concat();
    let mut bin= binary_n.iter().rev().collect_vec();
        bin.truncate(*width as usize);
        bin.reverse();
    let out :Vec<u64>=(0..*width)
        .rev()
        .into_iter()
        .map(|i|
            bin.iter().rev().nth(i as usize).unwrap()).
            into_iter()
            .map(
                |x| if **x == 1 { 1u64 } 
                    else { 0u64 })
            .collect();
    out
}
pub fn from_tag_u1024_to_vec(tag:&U1024,n:&u64) ->Vec<Vec<u64>>{
    let mut tag = tag.clone();
    let mut tag_vec:Vec<Vec<u64>>=Vec::new();
    let power_set_cardinality = U1024::one()<<*n;
    let mut i = 0;
    loop{
        if tag==U1024::zero() {break;}
        i+=1;
        if tag.trailing_zeros()>=*n as u32{
            panic!("Can't be an hypergroupoid. {:b} zeroes found, step {}",tag,i);
        } else if tag>=power_set_cardinality{
            tag_vec.push(collect_n_digits_u1024(n, &tag));
            tag>>=*n;
        }
        if tag<power_set_cardinality {
            tag_vec.push(collect_n_digits_u1024(n, &tag));
        
            break;}
    }
    tag_vec.iter().map(|x|x.clone()).rev().collect()  
}
pub fn from_tag_to_vec(tag:&u128, n:&u64)->Vec<Vec<u64>>{
    let mut tag = tag.clone();
    let mut tag_vec:Vec<Vec<u64>>=Vec::new();
    let power_set_cardinality = 2u128.pow(*n as u32);

    loop{
        if tag.trailing_zeros()>=*n as u32{
            panic!("Can't be an hypergroupoid. {} zeroes found",tag.trailing_zeros());
        } else {
            tag_vec.push(collect_n_digits(n, &tag));
            tag>>=*n as u128;
        }
        if tag<power_set_cardinality {
            tag_vec.push(collect_n_digits(n, &tag));

            break;}
    }
    tag_vec.iter().map(|x|x.clone()).rev().collect()
}
pub fn write(s:String,name:&str)-> std::io::Result<()> {
    let file_name=format!("{}.txt",name);

    let mut file = File::create(file_name)?;
    
    file.write(&s.as_bytes())?;
    Ok(())
}
pub fn get_min_max(cardinality:&u64)->(u128,u128){
    let size= cardinality.pow(3);
    let square= cardinality.pow(2);
let min: u128 = (0..square).into_iter().fold(0, |acc: u128,x: u64|acc+2u128.pow((cardinality*x) as u32));

let max = 2u128.pow(size as u32)-1;
(min,max)
 }
 pub fn get_min_max_u1024(cardinality:&u64)->(U1024,U1024){
    let size= cardinality.pow(3);
    let square= cardinality.pow(2);
let min: U1024 = (0..square).into_iter().fold(U1024::zero(), |acc: U1024,x: u64|acc+U1024::from(2).pow(U1024::from(cardinality*x)));

let max = U1024::from(2).pow(U1024::from(size))-1;
(min,max)
 }

