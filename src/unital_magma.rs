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
    pub identity: u32,
}

impl UnitalMagma {
    pub fn new_from_tag(tag:&u128,cardinality:&u32)->Self {
        assert!(representing_hypergroupoid(&mut tag.clone(), &cardinality),"Tag doesn't represent a hypergroupoid!");
        let h=HyperGroupoidMat::new_from_tag(&tag, &cardinality);
        let e = h.collect_scalar_identity();
        if e.len()!=1 {panic!("Not representing a unital magmata. No scalar identity found!")}
        let identity=e[0];
        UnitalMagma { h, identity: identity }
    }
    pub fn is_left_invertible(&self,x:&u32)->bool{
        /*x=1,2,3... is element in H, and they identify the singleton {i},
        therefore, as a subset of H, it is represented by the u32 given by 2.pow(x)
        */
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let col_x=self.h.hyper_composition.column(x as usize);
        col_x.iter().any(|y|self.identity&y==1) //self.identity&y==1 implies identity is in y. 
    
    }
    pub fn collect_left_invertible(&self)->Vec<u32>{
        self.h.get_singleton().iter().filter(|x|(self.is_left_invertible(x))).map(|x|*x).collect::<Vec<u32>>()
    }
    pub fn is_right_invertible(&self,x:&u32)->bool{
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let row_x=self.h.hyper_composition.column(x as usize);
        row_x.iter().any(|y|self.identity&y==1) //self.identity&y==1 implies identity is in y.    
        }
    pub fn collect_right_invertible(&self)->Vec<u32>{
        self.h.get_singleton().iter().filter(|x|(self.is_right_invertible(x))).map(|x|*x).collect::<Vec<u32>>()
    }
    pub fn is_invertible(&self,x: &u32)->bool{
        self.is_left_invertible(x)&&self.is_right_invertible(x)
    }
    pub fn collect_invertible(&self)->Vec<u32>{
        self.h.get_singleton().iter().filter(|x|(self.is_invertible(x))).map(|x|*x).collect::<Vec<u32>>()
    }
    
    
}
impl Display for UnitalMagma{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table:DMatrix<String>=DMatrix::from_iterator(self.h.n as usize, self.h.n as usize, 
            self.h.hyper_composition.iter().map(|x|format!("{:?}",vec_to_set(&get_subset(x, &self.h.n)))));
        let identity=vec_to_set(&get_subset(&self.identity, &self.h.n));
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {} Identity is {:?} Size:{}\n", self.h.h, table, self.h.hyper_composition,identity, self.h.n )
    }
}

