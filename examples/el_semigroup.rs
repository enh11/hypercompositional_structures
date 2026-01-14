use std::collections::HashSet;

use hyperstruc::{binary_relations::relations::Relation, hs::{hypergroupoids::{HyperGroupoid, QuotientHyperGroupoid}, ordered_semigroup::PreOrderedSemigroup}};
use itertools::Itertools;


fn main() {
let cardinality = 3u64;
 let semigroup_table = [
    0, 0, 0,
    0, 1, 0,
    0, 0, 2,
].iter().map(|el|vec![*el]).collect_vec();
 let semigroup=HyperGroupoid::new_from_elements(&semigroup_table,&cardinality);
 let r = Relation { 
     a: (0..cardinality).collect::<HashSet<u64>>(), 
    b: (0..cardinality).collect::<HashSet<u64>>(), 
     rel: [
        (0, 0),
        (1, 0), (1, 1), 
                (2, 2) 
                         ].to_vec() };

let preordered_semigroup = PreOrderedSemigroup::new(&semigroup,&r).unwrap();
let el_hg = HyperGroupoid::new_el_hypergroup(preordered_semigroup);
el_hg.0.show();
let t = Relation::new_from_preorder(&r);
let q = QuotientHyperGroupoid::new_from_regular_relation(&el_hg.0, &t.unwrap());
println!("{}",q.unwrap());
 }
