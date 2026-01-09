use std::collections::HashSet;

use hyperstruc::{binary_relations::relations::Relation, hs::{hypergroupoids::{HyperGroupoid, QuotientHyperGroupoid}, ordered_semigroup::PreOrderedSemigroup}};
use itertools::Itertools;


fn main() {
let cardinality = 4u64;
 let semigroup_table = [
    0, 0, 0, 0,
    0, 1, 1, 1,
    0, 1, 2, 3,
    0, 1, 3, 2
].iter().map(|el|vec![*el]).collect_vec();
 let semigroup=HyperGroupoid::new_from_elements(&semigroup_table,&cardinality);
 let r = Relation { 
     a: (0..cardinality).collect::<HashSet<u64>>(), 
    b: (0..cardinality).collect::<HashSet<u64>>(), 
     rel: [
         (0, 0), (0, 1), (0, 2), (0, 3), 
         (1, 0), (1, 1), (1, 2), (1, 3), 
                         (2, 2), (2, 3), 
                         (3, 2), (3, 3)].to_vec() };
let transp = r.antiisomorphic();

let preordered_semigroup = PreOrderedSemigroup::new(&semigroup,&r).unwrap();
let el_hg = HyperGroupoid::new_el_hypergroup(preordered_semigroup);
let equivalence = Relation::new_from_preorder(&r);
let quotient = QuotientHyperGroupoid::new_from_regular_relation(&el_hg.1.semigrp, &equivalence.unwrap());
match quotient {
    Ok(q) => println!("quotient {}",q),
    Err(e) => println!("{}",e),
}; 
}
