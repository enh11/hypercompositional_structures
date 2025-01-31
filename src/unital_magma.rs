use core::panic;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::{collections::HashSet, fmt::Display, vec};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use rand:: Rng;
use crate::{hs::HyperGroupoidMat, utilities::{binary_to_u32, cartesian_product, from_tag_to_vec, get_subset, n_to_binary_vec, ones_positions, permutaton_matrix_from_permutation, representation_permutation_subset, representing_hypergroupoid, subset_as_u32, vec_to_set}};
#[derive(Debug, Clone,PartialEq)]
pub struct UnitalMagma{
    pub h:HyperGroupoidMat,
    pub identity: HashSet<u32>,
}

impl UnitalMagma {
    pub fn new_from_tag(tag:&u128,cardinality:&u32)->Self {
        assert!(representing_hypergroupoid(&mut tag.clone(), &cardinality),"Tag doesn't represent a hypergroupoid!");
        println!("tag is {} cardinality is {}",tag,cardinality);
        let h=HyperGroupoidMat::new_from_tag(&tag, &cardinality);
        println!("h is ok");
        let e = h.collect_scalar_identity();
        if e.len()!=1 {panic!("Not representing a unital magmata. No scalar identity found!")}
        let identity=vec_to_set(&get_subset(&e[0], &cardinality));
        UnitalMagma { h, identity: identity }
    }
}
impl Display for UnitalMagma{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table:DMatrix<String>=DMatrix::from_iterator(self.h.n as usize, self.h.n as usize, 
            self.h.hyper_composition.iter().map(|x|format!("{:?}",vec_to_set(&get_subset(x, &self.h.n)))));
        
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {} Identity is {:?} Size:{}\n", self.h.h, table, self.h.hyper_composition,self.identity, self.h.n )
    }
}

