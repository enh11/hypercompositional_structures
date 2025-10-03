use std::fmt::Display;
use crate::{hs::QuotientHyperGroupoid, hypergroups::HyperGroup, binary_relations::relations::Relation};
extern crate nalgebra as na;

#[derive(Debug,Clone,PartialEq)]
pub struct QuotientHyperGroup(pub QuotientHyperGroupoid);
impl QuotientHyperGroup {
    pub fn new_from_equivalence_relation(hg:HyperGroup,equivalence:Relation)->Self{
        assert!(equivalence.is_equivalence(),"Error, {:?} is not an equivalence!",equivalence.rel);
        let hs = hg.0;
        let hs_quotient = QuotientHyperGroupoid::new_from_equivalence_relation(&hs, &equivalence);

        QuotientHyperGroup(hs_quotient)
    }
}
impl Display for QuotientHyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
}