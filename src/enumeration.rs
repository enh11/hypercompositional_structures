use std::fmt::format;
use std::ops::Index;
use itertools::Itertools;
use permutation::Permutation;
use rayon::iter::plumbing::{bridge, Consumer, ProducerCallback, UnindexedConsumer};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use crate::unital_magma::UnitalMagma;
use crate::utilities::{get_min_max, write, U1024};
use crate::hs::HyperGroupoidMat;
use crate::utilities::representing_hypergroupoid;


pub fn collect_beta_equivalence(cardinality:&u64)->Vec<u128>{
    let size = cardinality.pow(3);
    (2u128.pow((size-cardinality) as u32)..2u128.pow(size as u32)).
        into_par_iter().filter(|i|
            representing_hypergroupoid(i,&cardinality)
            &&
            HyperGroupoidMat::new_from_tag(i, cardinality).beta_relation()
                .is_transitive()
        )
        .collect()
}
pub fn collect_hypergroupoid(cardinality:&u64)->Vec<u128>{
    let size = cardinality.pow(3);
    (2u128.pow((size-cardinality) as u32)..2u128.pow(size as u32))
        .into_par_iter()
        .filter(|i|
                representing_hypergroupoid(i,&cardinality)
        )
        .collect()
}
pub fn collect_hypergroupoid_with_scalar_identity(cardinality:&u64)->Vec<u128>{
    let size = cardinality.pow(3);
    (2u128.pow((size-cardinality) as u32)..2u128.pow(size as u32))
        .into_par_iter()
        .filter(|i|
            representing_hypergroupoid(&i,&cardinality)
            &&
            !(HyperGroupoidMat::new_from_tag(i, &cardinality).collect_scalar_identity().is_empty())
        )
        .collect()

}
pub fn collect_invertible_magmata(cardinality:&u64)->Vec<u128>{
    //find a way to improve this my looking at binary representation of tag
    let size = cardinality.pow(3);
    (2u128.pow((size-cardinality) as u32)..2u128.pow(size as u32))
        .into_par_iter()
        .filter(|i|
            representing_hypergroupoid(&mut i.clone(),&cardinality)
            &&
            HyperGroupoidMat::new_from_tag(i, &cardinality).collect_scalar_identity().len()==1
            &&
            UnitalMagma::new_from_tag(i, &cardinality).is_invertible_unital_magma()
        )
        .collect()

}
pub fn collect_hypergroups(cardinality:&u64)->Vec<u128>{
    let (min,max)= get_min_max(cardinality);
    (min..=max).into_par_iter()
        .filter(|i|
            representing_hypergroupoid(&mut i.clone(),&cardinality)&&HyperGroupoidMat::new_from_tag(i, cardinality).is_hypergroup()
        )
        .collect()
}
pub fn enumeration_hyperstructure(structure:&str,cardinality:&u64)->Vec<usize>{
    let tags= match structure {
        "hypergroups"=> collect_hypergroups(&cardinality),
        "unital magmata"=>collect_hypergroupoid_with_scalar_identity(&*cardinality),
        "invertible magmata"=> collect_invertible_magmata(&cardinality),
        _=>panic!("unknown structure! Works with 'hypergroups, unital magmata,invertible magmata, L_mosaics'. ")
    };
    //let tags = collect_hypergroups(&cardinality);
    let _= write(format!("{:?}",tags.clone()),&format!("tag_{structure}_{cardinality}"));
    let permut_vec:Vec<Vec<usize>> = (0..*cardinality as usize).permutations(*cardinality as usize).collect();
    let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
    let mut classes:Vec<(u64,Vec<u64>)>=Vec::new();

    for tag in tags {
    let mut isomorphism_classes:Vec<u64>=Vec::new();

    for sigma in &permutation {        
        let isomorphic_image_tag = HyperGroupoidMat::new_from_tag(&tag, &cardinality).isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag();
        isomorphism_classes.push(isomorphic_image_tag.try_into().unwrap());

    }
    isomorphism_classes=isomorphism_classes.iter().sorted().dedup().map(|x|*x).collect();
    let representant_of_class=isomorphism_classes.iter().min().unwrap();
    
        classes.push((*representant_of_class,isomorphism_classes));

    
}
    let classes:Vec<&(u64, Vec<u64>)>=classes.iter().sorted_by(|x,y|x.0.cmp(&y.0)).dedup().collect();
    let mut c:Vec<usize>=Vec::new();
    let mut c_k:Vec<&(u64,Vec<u64>)>;
    let mut s = String::new();
    for k in 1..=permut_vec.len(){
        c_k=classes.iter().filter(|y|(*y.1).len()==k).into_iter().map(|x|*x).collect_vec();
        c.push(c_k.len());
        let add_str=format!("{:?}\n",c_k);
        s.push_str(&add_str);
        let _ = write(s.clone(),&format!("enumeration_{structure}_{cardinality}"));
        
    }
    c
}
