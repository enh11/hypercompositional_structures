#![allow(unused)]
#![allow(unused_imports)]
#![allow(dead_code)]
use std::fs::File;
use std::io::prelude::*;
use hyperstruc::hs::{collect_hypergroupoid, collect_hypergroups};
use hyperstruc::utilities::{binary_to_u32, cartesian_product, collect_n_digits, from_tag_to_vec, get_subset, n_to_binary_vec, ones_positions, permutaton_matrix_from_permutation, power_set, to_set};
use hyperstruc::{hs::{get_random_hypercomposition_matrix, HyperGroupoidMat}, hyper_structure::{representation_random_hypercomposition_table, HyperStruct}};
use itertools::interleave;
use nalgebra::DMatrix;
use itertools::Itertools;
use rand::{seq::index, Rng};
use std::collections:: HashSet;
use permutation::Permutation;

fn main(){
    
    /*COLLECT ORDER THREE HYPERGROUPOIDS */
    let cardinality = 4u32;
    
    let mut hg_tags:Vec<u128>=collect_hypergroups(&cardinality);
    
    println!("Funded {} hypergroups", hg_tags.len());
    let s = format!("{:?}",hg_tags);
    write(s);
    
    
/* 
/*GETTING A HYPERSTRUCTUR FROM INTEGER TAG */
let mut rng = rand::thread_rng();
let n: u32=3;
let mut tag: u128 = rng.gen_range(2u128.pow(26)+1..=2u128.pow(27));
println!("{:b}",tag);
let subsets = from_tag_to_vec(&mut tag, &n);
println!("{:?}",subsets);
println!("tag is {}",tag);
let hyperstructure = HyperGroupoidMat::new_from_tag(tag,&n);
println!("{}",hyperstructure);
let tag = hyperstructure.get_integer_tag();
println!("tag is {}",tag);
*/
/*
let h_groupoid=  HyperGroupoidMat::new_random_from_cardinality(&n);
println!(" A new Hyper Groupoid : {}",h_groupoid);
println!("H is reproductive: {}",h_groupoid.is_reproductive());
let new_hg=h_groupoid.fix_reproductivity();
println!("A reproductive Hypergroupoid: {}",new_hg);
println!("H is repdocutive: {}",new_hg.is_reproductive());
println!("H is associativity: {}",new_hg.is_associative());

 */

/*GET HYPERSTRUCTURE FROM MATRIX */
/*
let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
let hypergroup=HyperGroupoidMat::new_from_matrix(&matrix);
println!("{}",hypergroup);
println!("H is hypergroup: {}",hypergroup.is_hypergroup());
println!("H is commutative: {}",hypergroup.is_commutative());
/*This is the set of all permutation of n */
let permut_vec:Vec<Vec<usize>> = (0..hypergroup.n as usize).permutations(hypergroup.n as usize).collect();
let alpha = Permutation::oneline(permut_vec[3].clone());
println!("alpha is {:?}",alpha);
/*This is the matrix corresponding to the permutation aplha */
let perm_mat = permutaton_matrix_from_permutation(&hypergroup.n, &alpha);
println!("permutation prova {}",perm_mat);

let prova_permut = hypergroup.permutation_of_table(&alpha);
println!("prova permut {}",prova_permut);
 */

/* 
(0..hypergroup.n as usize).into_iter().map(|i| permutation_matrix_from_sigma.row(i)=identity.row(sigma.apply_idx(i))).collect();
 */

/* 
Look for a random associative hypergroupoid of order 5


let n = 5u32;
let mut semihypergroup=   HyperGroupoidMat::new_random_from_cardinality(&n);

semihypergroup=loop {
    let x=HyperGroupoidMat::new_random_from_cardinality(&n);

    if semihypergroup.is_associative(){break x;}


}; */

/* 
let h=HyperGroupoid::new_random_from_cardinality(&n);
println!("A new hyperstructure: H = {}",h);
println!("H is reproductive: {}",h.is_reproductive());
println!("H is associative: {}",h.is_associative());

let new_h =h.fix_reproductivity();
println!("A modified H:\n{}",new_h);
println!("H is reproductive: {}",new_h.is_reproductive());

println!("H is associative: {}",new_h.is_associative());
let new_new = new_h.fix_associativity();
println!("A modified H:\n{}",new_new);
println!("H is reproductive: {}",new_new.is_reproductive());
println!("H is associative: {}",new_new.is_associative());
new_new.is_associative();
new_new.check_associativity();

 */
println!("THE END\n");


}
fn write(s:String)-> std::io::Result<()> {
    let mut file = File::create("foo.txt")?;

    file.write(&s.as_bytes())?;
    Ok(())
}
fn singleton(v:&Vec<u32>)->Vec<Vec<u32>>{
    v.iter().map(|x| vec![*x]).collect()
}

fn random_hypercomposition_table(h:&Vec<u32>)->Vec<((u32,u32),Vec<u32>)>{
    /*
    This builds a random Cayley table containing the products
    ab for any a,b in H  */
    let cartesian = cartesian_product(&h);
    let mut hyper_operation_table: Vec<((u32,u32),Vec<u32>)>=Vec::new();
    let mut rng = rand::thread_rng();
    let n=h.len() as u32;
    let k  =2u32.pow(n.try_into().unwrap());

    for item in cartesian {
        let n_subset: u32 =rng.gen_range(0u32..k);
        let ab=get_subset(&n_subset, &n);
        hyper_operation_table.push((item,ab) );
    }
    hyper_operation_table
}
fn hypercomposition(ht:Vec<((u32,u32),Vec<u32>)>)->Vec<Vec<u32>>{
    /*
    Elements here are elements of the hyperproduct table. 
    To access the element ab use ht[3a+b]*/
    ht.iter().map(|x|x.1.clone()).collect()
}
fn set_hypercomposition(h:&Vec<u32>,a:&Vec<u32>,b:&Vec<u32>)->HashSet<u32>{
    /*
    Let A,B be subset of H. Then, the product AB is define as the union of 
    ab with a in A and b in B */
    let h_set=to_set(&h);
    let a_set=to_set(&a);
    let b_set=to_set(&b);

    if !a_set.is_subset(&h_set)|!b_set.is_subset(&h_set){panic!("A and B must be subsets of H")}
    
    let ht=random_hypercomposition_table(h);
    let hyper_prod= hypercomposition(ht);
    let mut prod: Vec<_>=Vec::new();
    for x in a {
        for y in b{
            let c =hyper_prod.get((3*x+y) as usize).expect("No elements here");
    
            prod.push(c.to_vec());
            prod.concat();
        }
    }
    to_set(&prod.concat())
}