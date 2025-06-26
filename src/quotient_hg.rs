use std::fmt::{self, Display};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use num_rational::Rational64;
use permutation::Permutation;
use rayon::{iter::{IntoParallelRefIterator, ParallelIterator}};
use crate::{fuzzy::FuzzySubset, hs::{circumference_radius_d_filtered, hg_in_circumference_radius_one, HyperGroupoid, QuotientHyperGroupoid}, hypergroups::HyperGroup, relations::Relation, utilities::{chi_a, get_complement_subset, U1024}};

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