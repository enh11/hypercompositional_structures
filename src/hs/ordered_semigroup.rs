use std::fmt::Display;
use itertools::Itertools;
use permutation::Permutation;
use rayon::{iter::{IntoParallelRefIterator, ParallelIterator}, slice::ParallelSliceMut};

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
pub fn antiisomoprhic(&self)->Self{
    let anti_semigroup = self.semigrp.antiisomoprhic();
    let anti_relation = self.order.antiisomorphic();
    PreOrderedSemigroup { semigrp: anti_semigroup, order:anti_relation }
}
pub fn isomorphic_preordered_semigroup_from_permutation(&self, sigma:&Permutation)->Result<PreOrderedSemigroup, HyperStructureError>{
    let isomorphic_semigroup = self.semigrp.isomorphic_hypergroup_from_permutation(sigma);
    let isomorphic_preorder = self.order.isomorphic_relation_from_permutation(sigma);

    PreOrderedSemigroup::new(&isomorphic_semigroup, &isomorphic_preorder)
}

pub fn collect_antiisomorphism_class(&self)-> Result<(Self, Vec<Self>), HyperStructureError>{
    self.antiisomoprhic().collect_isomorphism_class()
}
/// Anti-isomorphism class = Iso(self) ∪ Iso(anti(self))
pub fn collect_equivalence_classe(&self)-> Result<(Self, Vec<Self>), HyperStructureError>
{
    let iso  = self.collect_isomorphism_class();
    let anti = self.collect_antiisomorphism_class();

    match (iso, anti) {
        (Ok(mut u), Ok(mut v)) => {
            // merge
            u.1.append(&mut v.1);
            u.1.sort();
            u.1.dedup();

            // canonical representative = first element in sorted list
            let representant = u.1[0].clone();
            Ok((representant, u.1))
        }

        (Ok(u), Err(_)) => Ok(u),          // no anti-isomorphic structures
        (Err(_), Ok(u)) => Ok(u),          // no iso structures (rare)
        (Err(e), Err(_)) => Err(e),
    }
}
pub fn collect_isomorphism_class(&self)-> Result<(Self, Vec<Self>), HyperStructureError>{
    let cardinality = self.semigrp.n as usize;
    let permutation_vec:Vec<Vec<usize>> = (0..cardinality).permutations(cardinality ).collect();
    let permutation:Vec<Permutation> = permutation_vec.par_iter()
        .map(|sigma| 
                Permutation::oneline(sigma.clone())
            ).collect();
    let mut isomorphism_classes:Vec<_>=permutation
        .par_iter()
        .map(|sigma|
                self.isomorphic_preordered_semigroup_from_permutation(&sigma)
            )
        .filter_map(|s|s.ok()
        ).collect();
   isomorphism_classes.par_sort();
   isomorphism_classes.dedup();
   let representant_of_class= isomorphism_classes[0].clone();
        
       Ok((representant_of_class,isomorphism_classes))        

}    
pub fn collect_extended_isomorphism_class(&self)-> Result<(Self, Vec<Self>), HyperStructureError>{
    let class = self.collect_isomorphism_class();
    match class {
        Ok(mut class) => {
            let anti_relation = self.order.antiisomorphic();
            let s = PreOrderedSemigroup::new(&self.semigrp, &anti_relation);
            match s {
                Ok(s) => {
                    class.1.append(&mut s.collect_isomorphism_class().unwrap().1);
                    class.1.sort();
                    class.1.dedup();
                    Ok((class.1[0].clone(),class.1))
                },
                Err(e) => Err(e),
            } 
        },
        Err(e) => Err(e)
    }
}
}
    

