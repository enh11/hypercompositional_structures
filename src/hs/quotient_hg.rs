use std::fmt::Display;
use crate::{binary_relations::relations::Relation, hs::{HyperStructureError,hypergroupoids::QuotientHyperGroupoid, hypergroups::HyperGroup, quotient_hg}};
extern crate nalgebra as na;

#[derive(Debug,Clone,PartialEq)]
pub struct QuotientHyperGroup(pub QuotientHyperGroupoid);
impl QuotientHyperGroup {
    pub fn new_from_equivalence_relation(hg:HyperGroup,equivalence:Relation)->Result<Self,HyperStructureError>{
        assert!(equivalence.is_equivalence(),"Error, {:?} is not an equivalence!",equivalence.rel);
        let hs = hg.0;
        match  QuotientHyperGroupoid::new_from_regular_relation(&hs, &equivalence){
            Ok(quotient_hg) => Ok(QuotientHyperGroup(quotient_hg)),
            Err(e) => Err(e),
        }
    }
}
impl Display for QuotientHyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
}