use std::collections::HashSet;

use hyperstruc::{binary_relations::relations::Relation, hs::{hypergroupoids::HyperGroupoid, ordered_semigroup::PreOrderedSemigroup}};
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
let preordered_semigroup = PreOrderedSemigroup::new(&semigroup,&r).unwrap();
println!("a preodered semigroup {}",preordered_semigroup);
let class = preordered_semigroup.collect_isomorphism_class();
let el_h_rep = HyperGroupoid::new_el_hypergroup(class.0);
println!("class of length {}",class.1.len());
let el_h = class.1.iter().map(|s: &PreOrderedSemigroup|HyperGroupoid::new_el_hypergroup(s.clone())).collect_vec();
el_h.iter().all(|s|HyperGroupoid::new_el_hypergroup(s.1.clone())==el_h_rep);


}