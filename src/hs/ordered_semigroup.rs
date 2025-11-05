use std::fmt::Display;
use crate::{binary_relations::relations::Relation, hs::{hypergroupoids::HyperGroupoid, hypergroups::HyperStructureError}};


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
    
}

