use std::fmt::Display;
use crate::{hs::QuotientHyperGroupoid, hypergroups::HyperGroup, relations::Relation};
extern crate nalgebra as na;

#[derive(Debug,Clone,PartialEq)]
pub struct QuotientHyperGroup(pub QuotientHyperGroupoid);
impl QuotientHyperGroup {
    pub fn new_from_equivalence_relation(hg:HyperGroup,rel:Relation)->Self{
        let hs = hg.0;
        let hs_quotient = QuotientHyperGroupoid::new_from_equivalence_relation(&hs, &rel);
         println!("{}",hs_quotient);
        QuotientHyperGroup(hs_quotient)

    }
}
impl Display for QuotientHyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
}