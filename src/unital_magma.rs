use core::panic;
use std::fmt::Display;
extern crate nalgebra as na;
use nalgebra::DMatrix;
use crate::{hs::HyperGroupoid, utilities::{get_subset, representing_hypergroupoid, vec_to_set}};
#[derive(Debug, Clone,PartialEq)]
pub struct UnitalMagma{
    pub h:HyperGroupoid,
    pub identity: u64,
}

impl UnitalMagma {
    pub fn new_from_tag_u128(tag:&u128,cardinality:&u64)->Self {
        assert!(representing_hypergroupoid(&mut tag.clone(), &(*cardinality)),"Tag doesn't represent a hypergroupoid!");
        let h=HyperGroupoid::new_from_tag_u128(&tag, &cardinality);
        let e = h.collect_scalar_identities();
        if e.len()!=1 {panic!("Not representing a unital magmata. No scalar identity found!")}
        let identity=e[0];
        UnitalMagma { h, identity: identity }
    }
    pub fn is_unital_magma(&self)->bool{
        self.h.collect_identities().len()==1
    }
    pub fn is_left_invertible(&self,x:&u64)->bool{
        /*x=1,2,3... is element in H, and they identify the singleton {i},
        therefore, as a subset of H, it is represented by the u64 given by 2.pow(x)
        */
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let col_x=self.h.hyper_composition.column(x as usize);
        col_x.iter().any(|y|self.identity&y==self.identity) //self.identity&y==self.identity implies identity is in y. 
    
    }
    pub fn collect_left_inverses(&self,x:&u64)->Vec<u64>{
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let col_x=self.h.hyper_composition.column(x as usize);
        self.h.get_singleton().iter().filter(|i|col_x.get(i.trailing_zeros() as usize).unwrap()&self.identity==self.identity).map(|x|*x).collect() 
    }
    pub fn collect_left_invertible(&self)->Vec<u64>{
        self.h.get_singleton().iter().filter(|x|(self.is_left_invertible(x))).map(|x|*x).collect::<Vec<u64>>()
    }
    pub fn is_right_invertible(&self,x:&u64)->bool{
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let row_x=self.h.hyper_composition.row(x as usize);
        row_x.iter().any(|y|self.identity&y==self.identity) //self.identity&y==1 implies identity is in y.    
        }
    pub fn collect_right_inverses(&self,x:&u64)->Vec<u64>{
        if !x.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
        let x=x.trailing_zeros();//the number of trailing_zeros in a power of two is equal to the exponent. Also ilog2() works.
        let row_x=self.h.hyper_composition.row(x as usize);
        self.h.get_singleton().iter().filter(|i|row_x.get(i.trailing_zeros() as usize).unwrap()&self.identity==self.identity).map(|x|*x).collect() 
    }
    pub fn collect_inverses(&self,x:&u64)->Vec<u64> {
        assert!(self.is_invertible(x));
        let left_inverces=self.collect_left_inverses(x);
        let right_inverces=self.collect_right_inverses(x);
        left_inverces.iter().filter(|x|right_inverces.contains(x)).map(|x|*x).collect()
    }
    pub fn collect_right_invertible(&self)->Vec<u64>{
        self.h.get_singleton().iter().filter(|x|(self.is_right_invertible(x))).map(|x|*x).collect::<Vec<u64>>()
    }
    pub fn is_invertible(&self,x: &u64)->bool{
        self.is_left_invertible(x)&&self.is_right_invertible(x)
    }
    pub fn collect_invertible(&self)->Vec<u64>{
        self.h.get_singleton().iter().filter(|x|(self.is_invertible(x))).map(|x|*x).collect::<Vec<u64>>()
    }
    pub fn is_invertible_unital_magma(&self)->bool{
        self.h.get_singleton().iter().all(|x|self.is_invertible(x)&&self.collect_inverses(x).len()==1)
    }
    
}
impl Display for UnitalMagma{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table:DMatrix<String>=DMatrix::from_iterator(self.h.n as usize, self.h.n as usize, 
            self.h.hyper_composition.iter().map(|x|format!("{:?}",vec_to_set(&get_subset(x, &(self.h.n as u64))))));
        let identity=vec_to_set(&get_subset(&self.identity, &(self.h.n as u64)));
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {} Identity is {:?} Size:{}\n", self.h.h, table, self.h.hyper_composition,identity, self.h.n )
    }
}

