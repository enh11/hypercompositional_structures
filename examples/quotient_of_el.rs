use std::collections::HashSet;

use hyperstruc::{binary_relations::relations::{Relation, pre_order_enumeration}, hg_5::semigroups_5::SEMIGROUP_5, hs::{hypergroupoids::{HyperGroupoid, QuotientHyperGroupoid}, hypergroups::HyperGroup, ordered_semigroup::PreOrderedSemigroup}};
use itertools::Itertools;
use rand::Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};


fn main() {
    let cardinality =4;
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

let preordered_semigroup = PreOrderedSemigroup::new(&semigroup,&r).unwrap();
let el_hg = HyperGroupoid::new_el_hypergroup(preordered_semigroup);
let equivalence = Relation::new_from_preorder(&r);
let quotient = QuotientHyperGroupoid::new_from_regular_relation(&el_hg.1.semigrp, &equivalence.unwrap());
match quotient {
    Ok(q) => println!("quotient {}",q),
    Err(e) => println!("{}",e),
}; 
let cardinality = 5u64;
let pre5= pre_order_enumeration(&cardinality);
let semigroup_table = SEMIGROUP_5[800].iter().map(|el|vec![*el]).collect_vec();
let semigroup=HyperGroupoid::new_from_elements(&semigroup_table,&cardinality);
semigroup.show();
let compatible_preorders = pre5.into_par_iter().filter(|r| semigroup.is_relation_compatible(&r.into_relation())).map(|r|r.into_relation()).collect::<Vec<Relation>>();
let n  =compatible_preorders.len();
println!("n is {n}");
let mut rng = rand::rng();
let k = rng.random_range(1..n);
let preorder  = &compatible_preorders[k];
println!("preorder random compatible is {}",preorder.zero_one_matrix().0);
let semigrp_preorder = PreOrderedSemigroup::new(&semigroup, &preorder).unwrap();
let el = HyperGroupoid::new_el_hypergroup(semigrp_preorder);
println!("el is {}",el.0);
let mu_preorder = Relation::new_from_preorder(&preorder).unwrap();
let quotient_set = mu_preorder.quotient_set();
quotient_set.iter().for_each(|s|println!("{:?}",s));
println!("mu is {}",mu_preorder.zero_one_matrix().0);
let qu = QuotientHyperGroupoid::new_from_regular_relation(&el.0, &mu_preorder);
println!("quotient is {}",qu.unwrap());
}
