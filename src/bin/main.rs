use core::num;
use std::env;
use std::vec;
use std::time::Instant;
use hyperstruc::enumeration::collect_hypergroupoid_with_scalar_identity;
use hyperstruc::enumeration::enumeration_hyperstructure;
use hyperstruc::hg3::TAG_3;
use hyperstruc::hs::HyperGroupoidMat;
use hyperstruc::unital_magma::UnitalMagma;
use hyperstruc::unital_magma_3::TAG_UNITAL_MAGMATA_3;
use hyperstruc::utilities::binary_to_u32;
use hyperstruc::utilities::collect_n_digits;
use hyperstruc::utilities::from_tag_to_vec;
use hyperstruc::utilities::write;
use hyperstruc::utilities::{cartesian_product, get_subset, n_to_binary_vec, ones_positions, permutaton_matrix_from_permutation, power_set, representation_permutation_subset, subset_as_u32, vec_to_set};
use itertools::Itertools;
use nalgebra::coordinates::X;
use nalgebra::DMatrix;
use nalgebra::Matrix;
use rand::Rng;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use std::collections:: HashSet;
use permutation::Permutation;

fn main(){
/* 
let cardinality = 2;
let t=185;
println!("{:b}",t);
let new_hyperstructure_from_tag = HyperGroupoidMat::new_from_tag(&t,&cardinality);
let new_hyperstructure_from_matrix = HyperGroupoidMat::new_from_matrix(&DMatrix::from_row_slice(2usize,2usize,&[2,3,2,1]));
println!("{}", new_hyperstructure_from_tag);
let t= from_tag_to_vec(&t, &cardinality);

println!("{:?}",t);
println!("{}",new_hyperstructure_from_matrix);
println!("tag1 {}, tag2 {}", new_hyperstructure_from_matrix.get_integer_tag(),185); */
/*     let mat=DMatrix::from_row_slice(3, 3, &[1,2,4,2,5,7,4,2,1]);
    let h=HyperGroupoidMat::new_from_matrix(&mat);
    let magma=UnitalMagma{
        h:h,
        identity:1
    };

println!("unital magma: {}",magma.is_unital_magma());
let t  =24368401u128;
 let cardinality=3u32;
//let magma =UnitalMagma::new_from_tag(&t, &cardinality);
println!("magma {}",magma);
 println!("magma is invertible : {}",magma.is_invertible_unital_magma());
 let left_invertible= magma.collect_left_invertible();
 let left_inverses:Vec<(u32,Vec<u32>)>=left_invertible.iter().map(|x|(*x,magma.collect_left_inverses(x))).collect();

 let right_invertible=magma.collect_right_invertible();
 let right_inverses:Vec<(u32,Vec<u32>)>=right_invertible.iter().map(|x|(*x,magma.collect_right_inverses(x))).collect();

 println!("left inverses are {:?}",left_inverses);
 println!("right inverses are {:?}",right_inverses); */
/*
    let t=71663230u128;
    let cardinality=3u32;
    let m=UnitalMagma::new_from_tag(&t, &cardinality);
    let l_invertible=m.collect_left_invertible();
    let r_invertible=m.collect_right_invertible();


    println!("{m}");
    let l_inv_x=m.collect_left_inverses(&1u32);
    println!("left_ inverses of {} are {:?}",1u32.ilog2(),l_inv_x);
    let r_inv_x=m.collect_left_inverses(&1u32);
    println!("right_ inverses of {} are {:?}",1u32.ilog2(),r_inv_x);
    println!("left invertible are {:?}",l_invertible);
    println!("right invertible are {:?}",r_invertible);
    println!("m is invertible {}",m.is_invertible_unital_magma()); */
    

/*        /* COLLECT INVERTIBLE UNITAL MAGMATA (L-MOSAICS)*/ 
 let cardinality=2u32;
 let c= enumeration_hyperstructure("L_mosaics", &cardinality);
println!("c : {:?}",c);

 
 /* */
 

let cardinality=3u32;
let inv_magmata:Vec<u128> = TAG_UNITAL_MAGMATA_3.par_iter().filter(|m|UnitalMagma::new_from_tag(m, &cardinality).is_invertible_unital_magma()).map(|x|*x).collect::<Vec<u128>>();
    write(format!("{:?}",inv_magmata), &"invertible_unital_magmata_3");
println!("invertible magmata: {:?}",inv_magmata.len()); */

/* let cardinality=3u32;

for m in TAG_UNITAL_MAGMATA_3{
    let magma=UnitalMagma::new_from_tag(&m, &cardinality);
    if magma.is_invertible_unital_magma() {
        println!("{:b} with identity {:b}",m,magma.identity)
    }
  }
    
 */
    

         /*COLLECTING UNITAL MAGMATA */
/*     
let cardinality=2u32;
let now = Instant::now();
collect_hypergroupoid_with_scalar_identity(&cardinality); 
let end = now.elapsed();
    println!("Computation time:{:?}",end);
 */
      /*ISOMORPHIC HYPERGROUPS */
/*   let cardinality=3u32;
  let now = Instant::now();
let c= enumeration_hyperstructure("hypergroups", &cardinality);
let end = now.elapsed();
println!("Elapsed:{:?}",end);

println!("c : {:?}",c); */
/* 
let now = Instant::now();
    let e= enumeration_hypergroups(&4u32);
    println!("{:?}",e);
let end = now.elapsed();
println!("Elapsed:{:?}",end);
 */
      /*ISOMORPHIC MAGMATA */
/*   let cardinality=3u32;

let c= enumeration_hyperstructure("unital magmata", &cardinality);
println!("c : {:?}",c); */
/*
let now = Instant::now();
    let e= enumeration_hypergroups(&4u32);
    println!("{:?}",e);
let end = now.elapsed();
println!("Elapsed:{:?}",end);
 */ 
/*    let tag =25830028u128;
    let cardinality=3u32;
    println!("starting tag {}",tag);
    let hypergroup=HyperGroupoidMat::new_from_tag(&tag,&cardinality);
    println!("new from tag {}",hypergroup);
    println!("tag is hypergroup: {}",hypergroup.is_hypergroup());
    println!("tag {}",hypergroup.get_integer_tag());


for i in 0..hypergroup.n{
    if hypergroup.is_left_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u32.pow(i), &hypergroup.n));
        println!("Left identities {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_right_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u32.pow(i), &hypergroup.n));
        println!("Right identities {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u32.pow(i), &hypergroup.n));
        println!("Identity {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_left_scalar(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u32.pow(i), &hypergroup.n));
        println!("Left Scalar {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_right_scalar(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u32.pow(i), &hypergroup.n));
        println!("Right scalar {:?}",i_singleton)
    }
}

 */
/*     let args: Vec<String> = env::args().collect();  
    let number: u32 = match args[1].parse() {
        Ok(n) => {
            n
        },
        Err(_) => {
            eprintln!("error: Argument not an u32");
            return;
        },
    };    
    let now = Instant::now();
        let e= enumeration_hypergroups(&number);
        println!("{:?}",e);
    let end = now.elapsed();
    println!("Elapsed:{:?}",end);
 */
/*     let n=6u32;
    for i in 0..8 {

    let set =get_subset(&i, &3u32);
    println!("{} binary {:?}",i,n_to_binary_vec(&(i as u128), &3u32));
    println!("{} as vec is {:?}",i,set);
    let toset=to_set(&set);
    println!("{} to set is {:?}",i,toset);
    } */
/*     
    /*COLLECT ORDER THREE HYPERGROUPOIDS */
    
    let cardinality = 2u32;
    
    let mut hg_tags:Vec<u128>=collect_hypergroups(&cardinality);
    
    println!("Funded {} hypergroups", hg_tags.len());
    let s = format!("{:?}",hg_tags);
    write(s);
*/
    
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
/* 
GET HYPERSTRUCTURE FROM MATRIX */

let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
let hypergroup=HyperGroupoidMat::new_from_matrix(&matrix);
println!("{}",hypergroup);
println!("H is hypergroup: {}",hypergroup.is_hypergroup());



/*TEST NUMBER OF ISOMORPHISM IN TERMS OF PERMUTATIONS */
let mut count_isomorphism:u32=0;
let permut_vec:Vec<Vec<usize>> = (0..hypergroup.n as usize).permutations(hypergroup.n as usize).collect();
let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
for sigma in permutation {
    let alpha = sigma;
    let perm_mat = permutaton_matrix_from_permutation(&hypergroup.n, &alpha);
    let prova_permut = hypergroup.permutation_of_table(&alpha);
    let isomorphic_mat=perm_mat.clone()*prova_permut.hyper_composition*perm_mat.transpose();
    let isomorphic_hypergroup=HyperGroupoidMat::new_from_matrix(&isomorphic_mat);
    println!("sigma = {:?}",alpha);
    println!("isomorphic image {}",isomorphic_hypergroup);
    println!("isomorphic image is associative {}",isomorphic_hypergroup.is_associative());

    isomorphic_hypergroup.assert_associativity();
    if isomorphic_hypergroup.is_hypergroup() {
        count_isomorphism+=1;
    }

}
println!("number of isomorphism {}",count_isomorphism);

/* 
let now = Instant::now();
//let cardinality = 3;
let tag_2= collect_hypergroups(&CARDINALITY);
let permut_vec:Vec<Vec<usize>> = (0..CARDINALITY as usize).permutations(CARDINALITY as usize).collect();
let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
let mut classes:Vec<(u32,Vec<u32>)>=Vec::new();

for tag in tag_2 {
    let mut isomorphism_classes:Vec<u32>=Vec::new();

    for sigma in &permutation {        
        let isomorphic_image_tag = HyperGroupoidMat::new_from_tag(tag, &CARDINALITY).isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag();
        isomorphism_classes.push(isomorphic_image_tag);

    }
    let isomorphism_classes:Vec<u32>=isomorphism_classes.iter().sorted().dedup().map(|x|*x).collect();
    let representant_of_class=isomorphism_classes.iter().min().unwrap();
    
        classes.push((*representant_of_class,isomorphism_classes));

    
}
let classes:Vec<&(u32, Vec<u32>)>=classes.iter().sorted_by(|x,y|x.0.cmp(&y.0)).dedup().collect();
let mut c:Vec<usize>=Vec::new();
    let mut c_k:Vec<&(u32,Vec<u32>)>;
    let mut s = String::new();
    for k in 1..=permut_vec.len(){
        c_k=classes.iter().filter(|y|(*y.1).len()==k).into_iter().map(|x|x.clone()).collect_vec();
        c.push(c_k.len());
        let add_str=format!("{:?}\n",c_k);
        s.push_str(&add_str);
        write(s.clone());
        

    }
/* for tag in classes.iter().map(|x|x.0) {
    let hg  = HyperGroupoidMat::new_from_tag(tag as u128, &CARDINALITY);
    print!("{}",hg);
}  */

println!("isomorphism classes order 3 {:?}",c);
let end =now.elapsed();
println!("time= :{:?}",end);


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
    let h_set=vec_to_set(&h);
    let a_set=vec_to_set(&a);
    let b_set=vec_to_set(&b);

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
    vec_to_set(&prod.concat())
}