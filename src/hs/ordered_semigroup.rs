use std::fmt::Display;
use itertools::Itertools;
use permutation::Permutation;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{binary_relations::relations::Relation, hs::{HyperStructureError, hypergroupoids::HyperGroupoid}, utilities::{U1024, permutation_matrix_from_permutation}};


#[derive(Eq,PartialEq, PartialOrd,Ord,Clone)]
pub struct PreOrderedSemigroup {
    pub semigrp: HyperGroupoid,
    pub order: Relation
}
impl Display for PreOrderedSemigroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}{}", self.semigrp,self.order.zero_one_matrix().0)    }
}
impl PreOrderedSemigroup {
    pub fn new(semigroup:&HyperGroupoid,preorder:&Relation)->Result<Self,HyperStructureError>{
        match semigroup.is_semigroup() {
            true => match preorder.is_pre_order() {
                true => match semigroup.is_relation_compatible(&preorder) {
                    true => Ok(PreOrderedSemigroup{
                                    semigrp: semigroup.clone(),
                                    order: preorder.clone()}),
                    false => Err(HyperStructureError::NotPreOrderedSemigroup),
                    }
                false => Err(HyperStructureError::NotPreOrder),
            },
            false => Err(HyperStructureError::NotAssociative),
        }
    }
pub fn get_integer_tag_u1024(&self)->(U1024,U1024){
    (self.semigrp.get_integer_tag_u1024(),self.order.get_tag())
}
pub fn isomorphic_preordered_semigroup_from_permutation(&self, sigma:&Permutation)->Result<PreOrderedSemigroup, HyperStructureError>{
    let isomorphic_semigroup = self.semigrp.isomorphic_hypergroup_from_permutation(sigma);
    let isomorphic_preorder = self.order.isomorphic_relation_from_permutation(sigma);

    PreOrderedSemigroup::new(&isomorphic_semigroup, &isomorphic_preorder)
}

pub fn collect_isomorphism_class(&self)-> Result<(Self, Vec<Self>), HyperStructureError>{
    let cardinality = self.semigrp.n as usize;
    let permutation_vec:Vec<Vec<usize>> = (0..cardinality).permutations(cardinality ).collect();
    let permutation:Vec<Permutation> = permutation_vec.par_iter()
        .map(|sigma| 
                Permutation::oneline(sigma.clone())
            ).collect();
    let mut isomorphism_classes:Vec<_>=permutation.par_iter()
        .map(|sigma|
                self.isomorphic_preordered_semigroup_from_permutation(&sigma)
            ).filter_map(|s|s.ok()).collect();
   isomorphism_classes.sort();
   isomorphism_classes.dedup();
   let representant_of_class= isomorphism_classes[0].clone();
        
       Ok((representant_of_class,isomorphism_classes))


    /* let sg_class = self.semigrp.collect_isomorphism_class();

    let rel_class = self.order.collect_isomorphism_class();

    // Create the representant
    let representant = PreOrderedSemigroup::new(
        &HyperGroupoid::new_from_tag_u1024(&sg_class.0, &self.semigrp.n),
        &rel_class.0,
    )?; // <- replaces unwrap()

    // Build the class
    let mut class = sg_class
        .1
        .iter()
        .cartesian_product(rel_class.1.iter())
        .filter_map(|(s, r)| {
            PreOrderedSemigroup::new(
                &HyperGroupoid::new_from_tag_u1024(s, &self.semigrp.n),
                r,
            ).ok()     
        })
        .collect_vec();
    class.sort();
    class.dedup();
 let representant = &class[0];
    Ok((representant.clone(), class)) */
}

    
}

