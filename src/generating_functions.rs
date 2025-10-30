use std::collections::HashSet;

use crate::{binary_relations::relations::Relation, hs::hypergroupoids::HyperGroupoid, utilities::{subset_as_u64, support}};

pub fn b_hypercomposition() -> impl Fn(usize, usize) -> u64 {
    Box::new(move |a:usize,b:usize| {1<<a|1<<b})
}
pub fn tropical_hypergroup() -> impl Fn(usize, usize) -> u64 {
    Box::new(move |a: usize,b: usize| {
                    if a!=b {return 1<<a.max(b)}
                    else {return (0..=a).into_iter().fold(0, |acc,x|acc|1<<x)}
                
    })
}
pub fn genetics_hypergroup(cardinality: &u64) -> impl Fn(usize, usize) -> u64 {
    let h = (1u64 << cardinality) - 1;
    move |a: usize, b: usize| {
        let lower_elements = (0..a.min(b)).fold(0u64, |acc, x| acc | (1 << x));
        h - lower_elements
    }
}
pub fn el_hypergroup<'a>(semigroup:&'a HyperGroupoid,pre_order:&'a Relation) -> impl Fn(usize, usize) -> u64 + use<'a> {
    match semigroup.is_associative()&&semigroup.collect_scalars().len()==semigroup.n as usize {
        true => match pre_order.is_pre_order() {
        true => {
            move |a: usize, b: usize| {
                let supp_ab= support(&&semigroup.mul_by_representation(&(1<<a), &(1<<b)), &semigroup.n);
                let ab:HashSet<u64> = (0..semigroup.n).into_iter().filter(|x|
                    supp_ab.iter()
                        .any(|i|pre_order.are_in_relations(&(*i as u64), x))
                    ).collect();
                subset_as_u64(&ab)
                
    }
        },
        false => panic!("The input relation is not a pre-order."),
    },
        false => panic!("The input structure must be a semigroup."),
    }
    
    
}

