use hyperstruc::{binary_relations::relations::{Relation, pre_order_enumeration}, enumeration::enumeration_preordered_semigroup_from_list, hg_2::semigroup_2::{self, SEMIGROUP_2, TAGS_REPRESENTANTS_SEMIGROUPS_2}, hs::hypergroupoids::HyperGroupoid};
use itertools::Itertools;
use nalgebra::OMatrix;
use rayon::iter::ParallelIterator;


fn main() {
    let cardinality =2u64;
    let semigrp_2=SEMIGROUP_2.iter().map(|table|table.iter().map(|el|vec![*el]).collect_vec()).collect_vec();
    let semigrp_2 = semigrp_2.iter().map(|s|HyperGroupoid::new_from_elements(s, &cardinality)).collect_vec();

    //let semigrp_2:Vec<HyperGroupoid>= TAGS_REPRESENTANTS_SEMIGROUPS_2.iter().map(|s|HyperGroupoid::new_from_tag_u128(s, &cardinality)).collect();
    let rels:Vec<Relation> = pre_order_enumeration(&cardinality).map(|s|s.into_relation()).collect();
    let enumeration = enumeration_preordered_semigroup_from_list(&semigrp_2, &rels);
    let enumeration = match enumeration {
        Ok(enumeration) => {
           let total:usize= enumeration.iter().map(|s|s.1.len()).sum();
           println!("total number of preoder is {}",total);
           enumeration
        },
        Err(e) => panic!("{}",e),
    };
    for (s,compatible_relations) in enumeration {
        s.show();
        println!("is compatible with\n");
        compatible_relations.iter().for_each(|r|println!("{}",r.zero_one_matrix().0));
    }
    

}