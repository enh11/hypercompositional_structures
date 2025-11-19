
use itertools::Itertools;
use permutation::Permutation;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelBridge, ParallelIterator};
use rayon::vec;
use crate::binary_relations::relations::{Relation, pre_order_enumeration};
use crate::hg_2::semigroup_2::SEMIGROUP_2;
use crate::hg_3::semigroup_3::SEMIGROUP_3;
use crate::hg_4::semigroups_4::SEMIGROUP_4;
use crate::hg_5::semigroups_5::SEMIGROUP_5;
use crate::hs::HyperStructureError;
use crate::hs::hypergroupoids::HyperGroupoid;
use crate::hs::ordered_semigroup::PreOrderedSemigroup;
use crate::utilities::{get_min_max, get_min_max_u1024, parallel_tuples, representing_hypergroupoid_u1024, write, U1024RangeExt, U1024};
use crate::utilities::representing_hypergroupoid;


pub fn collect_beta_equivalence(cardinality:&u64)->Vec<u128>{
    let size = cardinality.pow(3);
    (2u128.pow((size-cardinality) as u32)..2u128.pow(size as u32)).
        into_par_iter().filter(|i|
            (representing_hypergroupoid(i,&cardinality))
            &
            (HyperGroupoid::new_from_tag_u128(i, cardinality).beta_relation().is_transitive())
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
            &
            !(HyperGroupoid::new_from_tag_u128(i, &cardinality).collect_scalar_identities().is_empty())
        )
        .collect()

}
pub fn collect_hypergroups(cardinality:&u64)->Vec<u128>{
    let (min,max)= get_min_max(cardinality);
    (min..=max).into_par_iter()
        .filter(|i|
            (representing_hypergroupoid(&mut i.clone(),&cardinality))
            &&
            (HyperGroupoid::new_from_tag_u128(i, cardinality).is_hypergroup())
        )
        .collect()
}
/// Collect all semigroup of cardinalty `n`. 
/// The procedure runs a brute-force in parallel. 
pub fn collect_semigroup(cardinality:&u64)->Vec<u128>{
    let p = parallel_tuples(*cardinality);
    p.map(|x|
        HyperGroupoid::new_from_elements(&x, cardinality))
        .filter(|hs|hs.is_associative())
        .map(|x|x.get_integer_tag()).collect()

}
pub fn enumeration_preordered_semigroup_from_list(semigroups:&Vec<HyperGroupoid>,relations:&Vec<Relation>)->Result<Vec<(HyperGroupoid,Vec<Relation>)>,HyperStructureError>{
    match semigroups.iter().all(|s|s.is_semigroup()) {
        true => match relations.iter().all(|r|r.is_pre_order()) {
            true => {
                let out=
                semigroups.par_iter().map(|s|{
                let compatible:Vec<Relation> = relations
                    .par_iter()
                    .filter(|r|
                        s.is_relation_compatible(&r)
                    ).map(|r|r.clone()
                    ).collect();
                (s.clone(),compatible)
                }
                ).collect();
                Ok(out)
            },
            false => Err(HyperStructureError::NotPreOrder),
        },
        false => Err(HyperStructureError::NotSemigroup),
    }

}
/*This works, but still to slow with respect to (u128..u128).into_par_iter() */
pub fn collect_hypergroups_u1024(cardinality:&u64)->Vec<U1024>{
    let (min,max)= get_min_max_u1024(cardinality);
    let hgs:Vec<U1024> = min.to(max+1).into_iter().par_bridge().filter(|i|
        (representing_hypergroupoid_u1024(i,&cardinality))
        &&
        (HyperGroupoid::new_from_tag_u1024(i, cardinality).is_hypergroup())
    )
    .collect();
let hgs:Vec<U1024>=hgs.iter().sorted().map(|x|*x).collect();
hgs
    
}
pub fn enumeration_hyperstructure(structure:&str,cardinality:&u64)->Vec<usize>{
    println!("Collecting all {} with cardinality {}...",structure,cardinality);
    let tags= match structure {
        "semigroups"=> collect_semigroup(&cardinality),
        "hypergroups"=> collect_hypergroups(&cardinality),
        "unital magmata"=>collect_hypergroupoid_with_scalar_identity(&*cardinality),
        _=>panic!("unknown structure! Works with 'hypergroups, unital magmata,invertible magmata, L_mosaics'. ")
    };
    let _= write(format!("{:?}",tags.clone()),&format!("tag_{structure}_{cardinality}"));
    let permut_vec:Vec<Vec<usize>> = (0..*cardinality as usize).permutations(*cardinality as usize).collect();
    println!("Completed.\nThere are {} {} of order {}",tags.len(),structure,cardinality);
    println!("Collecting classes of equivalence.");
    let classes:Vec<(U1024,Vec<U1024>)> = tags.iter()
        .map(|x|
            HyperGroupoid::new_from_tag_u128(x, cardinality).collect_isomorphism_class()
        )
        .unique()
        .collect();

    let c_k =
        (1..=permut_vec.len()).into_iter()
            .map(|k|
                classes.iter().filter(|x|x.1.len()==k).sorted_by(|x,y|x.0.cmp(&y.0)).collect_vec()
        ).collect_vec();

    let c:Vec<usize>= c_k.iter().map(|x|x.len()).collect();
    let s:String = c_k.iter().map(|x|format!("{:?}\n",x)).collect();
    let _ = write(s.clone(),&format!("enumeration_{}_{}",structure,cardinality));
    println!("Done!\nThe final enumeration is available in the file enumeration_{}_{}.txt",structure,cardinality);
    c
}
pub fn enumeration_hyperstructure_u1024(structure:&str,cardinality:&u64)->Vec<usize>{
    let tags= match structure {
        "hypergroups"=> collect_hypergroups_u1024(&cardinality),
        _=>panic!("unknown structure! Works with 'hypergroups, unital magmata,invertible magmata, L_mosaics'. ")
    };
    //let tags = collect_hypergroups(&cardinality);
    let _= write(format!("{:?}",tags.clone()),&format!("tag_{structure}_{cardinality}_{}","u1024"));
    let permut_vec:Vec<Vec<usize>> = (0..*cardinality as usize).permutations(*cardinality as usize).collect();
    let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
    let mut classes:Vec<(u64,Vec<u64>)>=Vec::new();

    for tag in tags {
    let mut isomorphism_classes:Vec<u64>=Vec::new();

    for sigma in &permutation {        
        let isomorphic_image_tag = HyperGroupoid::new_from_tag_u1024(&tag, &cardinality).isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag();
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
        let _ = write(s.clone(),&format!("enumeration_{structure}_{cardinality}_{}","u1024"));
        
    }
    c
}

// Seams to not work properly. It is ok up to the computation of preordered semigroup. 
// It seams to fail when computing classes up to isomorphism. 
// For example, wiith n = 2 we get 14 classes of equivalence of preordered semigroups, while 
// they are supposed to be 10.
//   

pub fn enumerate_el_hyperstructures(cardinality:&u64)->usize{
    let semigroups= match cardinality {
        2u64 => {SEMIGROUP_2.iter().map(|s|s.iter().map(|el|vec![*el]).collect_vec()).collect_vec()},
        3u64 => {SEMIGROUP_3.iter().map(|s|s.iter().map(|el|vec![*el]).collect_vec()).collect_vec()}
        _=>panic!("{}",HyperStructureError::ListOfSemigroupsNotAvailable)
        
    }; 
    let relations:Vec<Relation> = pre_order_enumeration(cardinality).map(|w|w.into_relation()).collect();
   
/*      /* Collect all preorders */
let relations :Vec<Relation>= pre_order_enumeration(&cardinality).map(|s|s.into_relation()).collect();
println!("there are {} preorder relations ", relations.len());
/* Collect all semigroups (if the list is not available) up to isomorphism */
let semigrps  = collect_semigroup(&cardinality)
    .iter()
    .map(|w|
        HyperGroupoid::new_from_tag_u128(w, &cardinality)
        )
    .map(|s|s.collect_isomorphism_class()
    ).sorted().dedup().collect_vec();
/* Store representants as Hypergroupoid  */
let semigroups = semigrps.iter().map(|s|HyperGroupoid::new_from_tag_u1024(&s.0, &cardinality)).collect_vec();
 */
    println!("there are {} possible semigroups " ,semigroups.len());
    let preordered_semigroup = semigroups.iter().cartesian_product(relations)
        .filter_map(|(s,r)|
/*         PreOrderedSemigroup::new(s, &r).ok()
 */            PreOrderedSemigroup::new(
        &HyperGroupoid::new_from_elements(s, cardinality), 
        &r).ok()
        
        )
        .sorted()
        .dedup()
        .collect_vec();
    println!("there are {} preodered semigroups of order {}",preordered_semigroup.len(),cardinality);
    let mut classes: Vec<(PreOrderedSemigroup, Vec<PreOrderedSemigroup>)>  = preordered_semigroup
        .par_iter()
        .map(|s|
            s.collect_isomorphism_class()
        ).filter_map(|s|
            s.ok()
        ).collect();
    classes.sort();
    classes.dedup();
println!("there are {} classes of preorders",classes.len());
/* println!("we print all classes of preorders");
for item in &classes {
    item.1.iter().for_each(|s|println!("{:?}",s.get_integer_tag_u1024()));
    println!("_______________________")
} */

/*
// This check that isomorphic preorders generate the same EL-hyperstrcuture
// This is a Theorem actually. But it is not an if and only if, i.e.,
// if two EL-hyperstructure coinced, then the generating preorders are not equal in general.

    classes.iter().all(|(rep,class)|


class.iter().all(|r|
    HyperGroupoid::new_el_hypergroup(r.clone())
    ==
    HyperGroupoid::new_el_hypergroup(rep.clone())
    )
); */

/* Fin qui tutto ok pare... 
 Da qui in poi quacosa non va ---
 */

    let el= classes
    .iter()
    .map(|s|
        HyperGroupoid::new_el_hypergroup(s.0.clone())
    )
    .sorted()
    .dedup()
    .collect_vec();

println!("The total number of EL-hyperstructures, considering the generating order is (not up to isomorphism) is {}",el.len());
let el_no_ordering = el.iter().map(|s|s.0.clone()).sorted().dedup().collect_vec();
println!("The total number of EL-hyperstructures, without considering the generating order is (not up to isomorphism) is {}",el_no_ordering.len());

/* Up to isomorphism */
let iso_el = el_no_ordering.iter().map(|s|s.collect_isomorphism_class()).sorted().dedup().collect_vec();
println!("The total number of EL-hyperstructures of order {} up to isomorphism is {}",cardinality,iso_el.len());
let out = iso_el.iter().filter(|s|HyperGroupoid::new_from_tag_u1024(&s.0,cardinality).is_hypergroup()).count();

    println!("there are {} EL-Hypergroups",out);
    /*for el in el{
        el.show();
    } */
out
    }
