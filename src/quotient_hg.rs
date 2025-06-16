use std::fmt::{self, Display};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use num_rational::Rational64;
use permutation::Permutation;
use rayon::{iter::{IntoParallelRefIterator, ParallelIterator}};
use crate::{fuzzy::FuzzySubset, hs::{circumference_radius_d_filtered, hg_in_circumference_radius_one, HyperGroupoid, QuotientHyperGroupoid}, utilities::{get_complement_subset, chi_a, U1024}};

#[derive(Debug, Clone,PartialEq)]
pub struct QuotientHyperGroup(QuotientHyperGroupoid);