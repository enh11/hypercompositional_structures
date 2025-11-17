use std::fmt::Display;
use itertools::Itertools;

use crate::{binary_relations::relations::Relation, hs::{hypergroupoids::HyperGroupoid, hypergroups::HyperStructureError}};


#[derive(Eq,PartialEq, PartialOrd,Ord,Clone)]
pub struct PreOrderedSemigroup {
    pub semigrp: HyperGroupoid,
    pub order: Relation
}
impl Display for PreOrderedSemigroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{},{}", self.semigrp,self.order.zero_one_matrix().0)    }
}
impl PreOrderedSemigroup {
    pub fn new(sg:&HyperGroupoid,r:&Relation)->Result<Self,HyperStructureError>{
        match sg.is_semigroup() {
            true => match r.is_pre_order() {
                true => match sg.is_relation_compatible(&r) {
                    true => Ok(PreOrderedSemigroup{
                                    semigrp: sg.clone(),
                                    order: r.clone()}),
                    false => Err(HyperStructureError::NotPreOrderedSemigroup),
                    }
                false => Err(HyperStructureError::NotPreOrder),
            },
            false => Err(HyperStructureError::NotAssociative),
        }
    }
    pub fn collect_isomorphism_class(&self)->(Self,Vec<Self>) {
        /* Seams to work but must be improved.. Avoid .unwrap()*/
        
        let sg_class= self.semigrp.collect_isomorphism_class();
        println!("There are {} isomorphic semigroup",sg_class.1.len());
        let rel_class = self.order.collect_isomorphism_class();
        
        println!("There are {} isomorphic relation",rel_class.1.len());
        
        let sg_representant = HyperGroupoid::new_from_tag_u1024(&sg_class.0, &self.semigrp.n);
        let rel_representant = rel_class.0;
        let representant = PreOrderedSemigroup::new(&sg_representant, &rel_representant);
       /*  match representant {
            Ok(rep) => println!("A representant for the class {}",rep),
            Err(e) => println!("{}",e),
        } */
        let mut class = sg_class.1.iter().cartesian_product(rel_class.1).map(|(s,r)| {
            
            let sg = HyperGroupoid::new_from_tag_u1024(s, &self.semigrp.n);
            PreOrderedSemigroup::new(&sg, &r)
        }).filter(|s|s.is_ok()).map(|a|a.unwrap()).collect_vec();
        class.sort();
        class.dedup();
        (representant.unwrap(),class)
        


    }
    
}

