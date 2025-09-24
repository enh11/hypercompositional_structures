use core::time;
use std::time::{Duration, Instant};

use hyperstruc::{hs::{HyperGroupoid, QuotientHyperGroupoid}, hypergroups::HyperGroup, quotient_hg::{self, QuotientHyperGroup}};
fn main(){
//Example Hv-group of order 10
let cardinality = 10u64;
let input_array1 = vec![
        vec![0],vec![1],vec![2],vec![3],vec![4],vec![0,5],vec![6],vec![7],vec![8],vec![9],
        vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],
        vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],
        vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],
        vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],
        vec![0,5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],
        vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],
        vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],
        vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],
        vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8]
    ];

let input_array2 = vec![
        vec![0],vec![0,1],vec![0,2],vec![0,3],
        vec![0],vec![0,1],vec![0,2],vec![0,3],
        vec![0],vec![0,1],vec![0,2],vec![0,3],
        vec![0],vec![0,1],vec![0,2],vec![0,3]
    ];/* let input_array = vec![
        vec![0],vec![1],vec![0,2],vec![3],
        vec![1],vec![0,1],vec![3],vec![2,3],
        vec![0,2],vec![3],vec![0,1,2],vec![1,3],
        vec![3],vec![2,3],vec![1,3],vec![0,1,2,3]
    ]; */
    let hs1 = HyperGroupoid::new_from_elements(&input_array1, &cardinality);
    println!("{}",hs1);
    assert!(hs1.is_hv_group());
    let x  =1u64<<0;
    let y = 1u64<<1;
    let x_right_y = hs1.right_division(&x, &y);
    println!("x right y = {}",x_right_y);
    let pil = hs1.collect_partial_left_identities();
    let pir  = hs1.collect_partial_right_identities();
    println!("H1 is associative {}",hs1.is_associative());
    println!("H1 is reproductive {}", hs1.is_reproductive());
    println!("is weak associative {}",hs1.is_weak_associative());
    println!("H1 partial left identities are {:?}",pil);
    println!("H1 partial right identities are {:?}",pir);
    println!("identities are {:?}",hs1.collect_identities());
    let beta = hs1.beta_relation();
    let beta_star =beta.transitive_closure();
    assert!(beta_star.is_transitive());
    let quotient_hg = QuotientHyperGroupoid::new_from_equivalence_relation(&hs1, &beta_star);
    println!("quotient {}",quotient_hg);
    let cardinality = 4u64;
    let hs2 = HyperGroupoid::new_from_elements(&input_array2, &cardinality);
    println!("{}",hs2);
    println!("H2 is associative {}",hs2.is_associative());
    println!("H2 is reproductive {}", hs2.is_reproductive());
    let pil = hs2.collect_partial_left_identities();
    let pir  = hs2.collect_partial_right_identities();
    println!("H2 partial left identities are {:?}",pil);
    println!("H2 partial right identities are {:?}",pir);


    //hypergroup of order 4
    let cardinality = 4u64;
    let input_array = vec![
        vec![0,1,2,3],vec![1,2],vec![0,3],vec![0,1,2],
        vec![1,2],vec![0,3],vec![0,1,2],vec![0,1,2,3],
        vec![0,3],vec![0,1,2],vec![0,1,2,3],vec![1,2],
        vec![0,1,2],vec![0,1,2,3],vec![1,2],vec![0,3]
    ];
    let hs = HyperGroupoid::new_from_elements(&input_array,& cardinality);
    println!("hypergroup {}",hs);
    println!("is hg {}",hs.is_hypergroup());
    let class = hs.beta_relation().get_class(&2u64);
    println!("CLASS {:?}",class);
    println!("identities {:?}",hs.collect_identities());

    
    let cardinality = 7u64;
    let input_array = vec![
        vec![0],vec![1],vec![2],vec![3],vec![4],vec![5,6],vec![5,6],
        vec![1],vec![0],vec![4],vec![5,6],vec![2],vec![3],vec![3],
        vec![2],vec![5,6],vec![0],vec![4],vec![3],vec![1],vec![1],
        vec![3],vec![4],vec![5,6],vec![0],vec![1],vec![2],vec![2],
        vec![4],vec![3],vec![1],vec![2],vec![5,6],vec![0],vec![0],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
    ];
    let hypergroupoid = HyperGroupoid::new_from_elements(&input_array, &cardinality);
    let hg = match hypergroupoid.is_hypergroup() {
        true => hyperstruc::hypergroups::HyperGroup::new_from_hypergroupiod(&hypergroupoid),
        false => panic!()
    };
    let now = Instant::now();
    let heart1 = hg.heart_fast();
    let end1 = now.elapsed();
    println!("time {:?}",now.elapsed());
        let now = Instant::now();
    let heart2 = hg.heart();
    let end2= now.elapsed();
    println!("time {:?}",now.elapsed());
    assert_eq!(heart2,heart1);
    println!("ratio {}",end2.as_secs_f64()/end1.as_secs_f64());
    let beta = hg.collect_beta_classes();
    println!("beta {:?}",beta);
    let fundamental = hg.get_isomorphic_fundamental_group();
    println!("f {}",fundamental);
    let heart = hg.heart();
    println!("heart is {:?}",heart);
    let left_identities = hg.collect_left_identities();
    println!("left id are {:?}",left_identities);
    let right_identities  = hg.collect_right_identities();
    println!("right id are {:?}",right_identities);
    let identities= hg.collect_identities();
    println!("id are {:?}",identities);
        let right_scalars= hg.collect_right_scalars();
    println!("right scalars are {:?}",right_scalars);
    let left_scalars = hg.collect_left_scalars();
    println!("left scalars are {:?}",left_scalars);
    let scalars = hg.collect_scalars();
    println!("scalars are {:?}",scalars);
    assert!(hg.0.is_left_partial_identity_of_x(&1u64, &4u64));
    let partial_identities = hg.collect_partial_identities();
    println!("partial identities {:?}",partial_identities);
    let sub_hg = hg.collect_proper_subhypergroups();
    println!("sub are {:?}",sub_hg);
    let closed = hg.collect_proper_closed_subhypergroups();
    println!("closed are {:?}",closed);
    let normals = hg.collect_proper_normal_subhypergroups();
    println!("normals {:?}",normals);
    let reflexive = hg.collect_proper_reflexive_subhypergroups();
    println!("reflexive are {:?}",reflexive);
        let invertible = hg.collect_proper_invertible_subhypergroups();
    println!("reflexive are {:?}",invertible);
    

 
  /*  let cardinality = 3u64;
    let h =(1u64<<cardinality)-1; 

    let x = 5u64;
    let y = 6u64;
    

    println!("x meet y is {}. In binary is {:b}",x&y, x&y);
    println!("x join y is {}. In binary is {:b}",x|y,x|y);

    for k in 0..=h{
    println!(
        "complement of {:0width$b} is {:0width$b}",
        k,
        k^h,
        width = cardinality as usize
    );    }

  
     
    /*EXAMPLE USING MATRIX */
    let cardinality =3u64;
    let input_values  = vec![vec![0,1,2],  vec![0],  vec![2],
                                            vec![0,2],vec![1,2],vec![2],
                                            vec![1],vec![0,1,2],vec![2]];
    let hs = HyperGroupoid::new_from_elements(&input_values, cardinality);
    println!("hs is {}",hs);

    let cardinality =4u64;
    let matrix=
        DMatrix::from_row_slice(cardinality as usize, cardinality as usize, 
            &[2,2,3,14,5,14,3,3,1,11,12,7,7,3,8,8]);
    let hs = HyperGroupoid::new_from_matrix(&matrix);
    println!("{}",hs);

    let par_id_right :Vec<HashSet<u64>> = hs.collect_partial_right_identities().iter().map(|x|u64_to_set(x, &cardinality)).collect();
    println!("partial right identities are {:?}",par_id_right);


    /*EXAMPLE USING GENERATING FUNCTION*/
    let cardinality =7u64;
    let hs = HyperGroupoid::new_from_function(b_hypercomposition(), &cardinality);
    println!("{}",hs);

    /*EXAMPLE USING GENERATING FUNCTION*/
    let cardinality = 5;
    let hs = HyperGroupoid::new_from_function(tropical_hypergroup(), &cardinality);
    println!("{}",hs);

    let tag = hs.get_integer_tag_u1024();
    println!("hs is identified by {}",tag);

    println!("is hypergroup {}",hs.is_hypergroup());

    let hg = HyperGroup::new_from_tag_u1024(&tag, &cardinality);
    
    println!("is commutative {}",hg.is_commutative());
    println!("is transposition {}",hg.is_transposition());
    
    let sub_hg = hg.collect_proper_subhypergroups();
        for item in sub_hg {
            let subset = u64_to_set(&item, &cardinality);
            println!("{:?} is a subhypergroup ",subset);
        }
let identities = hg.collect_identities();
for id in identities {
            let subset = u64_to_set(&id, &cardinality);
            println!("{:?} is an identity ",subset);
        }
let scalar_elements = hg.collect_scalars();
for sc in scalar_elements {
            let subset = u64_to_set(&sc, &cardinality);
            println!("{:?} is a scalar element ",subset);
        }


 */


/* 
    let cardinality =5u64;
        let function = {|a:u64,b:u64| 
            if a!=b {return 1<<a.max(b)}
            else {return (0..=a).into_iter().fold(0, |acc,x|acc|1<<x)}};
        let hg = HyperGroup::new_from_function(function, &cardinality).unwrap();
        assert!(hg.is_transposition());
        println!("hg is {}",hg);
    let cardinality = 4u64;
        let matrix = DMatrix::from_row_slice(cardinality as usize, cardinality as usize, &[1,1,7,11,1,1,7,11,7,7,7,12,11,11,12,11]);
        let hg = HyperGroup::new_from_matrix(&matrix);
        let subset_a = 1u64;
        let subset_b = 7u64;
        let subset_c = 11u64;
        println!("hg {}",hg);

        assert!(!hg.subhypergroup_is_closed(&subset_a));
        assert!(!hg.subhypergroup_is_closed(&subset_b));
        assert!(!hg.subhypergroup_is_closed(&subset_c));
    let cardinality = 3u64;
    let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,4,1,2,4,7,7,7]);   
    let hg = HyperGroup::new_from_matrix(&matrix);
    let id = hg.collect_identities();
    let id_set = id.iter().map(|x| u64_to_set(x, &cardinality)).collect_vec();
    println!("id is {:?}",id_set);

    assert!(id.is_empty());
    let tag =24476983u128;
    let hg = HyperGroup::new_from_tag_u128(&tag, &cardinality);
    let identities = hg.collect_identities();
    println!("hg {}",hg);
    println!("identities are {:?}",identities);
    for u in identities {
        for a in hg.get_singleton() {
            let left_inv = hg.left_inverses_of_x(&a, &u);
            println!("left inverses of {} with respect to {} are {:?}",a,u,left_inv);
            let right_inv = hg.right_inverses_of_x(&a, &u);
            println!("right inverses of {} with respect to {} are {:?}",a,u,right_inv)
    }
}
for a in hg.get_singleton(){
let inv = hg.collect_inverses_of_x(&a);
println!("inverses of {} are {:?}",a.trailing_zeros(),inv)
} */
/* let cardinality =5u64;
let function = |a:u64,b:u64| 1<<a|1<<b;
let hg = match HyperGroup::new_from_function(function, &cardinality) {
    Ok(hg)=>hg,
    Err(HyperStructureError::NotHypergroup)=> panic!("Not hg")
    
};
println!("is canonical {}",hg.is_transposition());
let normal = hg.collect_proper_normal_subhypergroups();
let closed = hg.collect_proper_closed_subhypergroups();
let normal_closed =normal.iter().filter(|x|closed.contains(&x)).collect_vec();

for item in normal {
    let to_set = vec_to_set(&get_subset(&item, &cardinality));
    println!("closed subhg are {:?}",item);
} */

/* let cardinality=3u64;
let hgs= TAG_3_REPRESENTANTS.into_iter().map(|x|HyperGroup::new_from_tag_u128(&x, &cardinality)).collect_vec();
let hg = hgs.iter().filter(|x|x.is_quasicanonical()).collect_vec();
for item in hg {
    println!("hs quasi canonical {}",item);
} */

/*     /*TRY TRANSPOSITION */
    let cardinality = 4u64;
    let matrix=DMatrix::from_row_slice(cardinality as usize,cardinality as usize,&[1,3,5,9,3,2,6,10,5,6,4,12,9,10,12,8]);
let hs = HyperGroupoid::new_from_matrix(&matrix);
let tag  = U1024::from(1970661744247085768u128);
let hg = HyperGroup::new_from_tag_u1024(&tag, &cardinality);
println!("hs {}",hs);
println!("hg {}",hg.is_transposition());    
let circ = circumference_radius_d_filtered(&tag, &1usize, &cardinality);
println!("circ {:?}",circ);
 */
/* let cardinality =2u64;
let tag =U1024::from(111);
 */
/* 
    let cardinality = 3u64;
    let cardinality = 3u64;
    let total = get_min_max_u1024(&cardinality).1;
let mut c =TAG_HG_3_CLASS_3.to_vec();
let mut  a = ALMOST_TAG_HG_3_CLASS_3.to_vec();
c.retain(|x|!a.contains(x));
println!("not found in a {:?}",c);
for item in c.iter().map(|x|x.0) {
    let dist  = distance_tags_u1024(&U1024::from(item), &total, &cardinality);
    println!("dist {} from {} is {}",item,total,dist)
}
 */
/* missing representants in classes of order 3 are
 [22102794, 22239147, 22369620, 22377812, 31021851, 31134601, 57784611, 57784615]*/
/*Representants 22102794, 22239147,31021851,31134601 have neighborhood hypergroups */
/*Representant 22369620 has neighborhood hypergroups [22369621, 22377812, 89478484] */
/*Representant 22377812 has neighborhood hypergroups [22377813, 22369620, 89486676] */
/*Representant 57784611 has neighborhood hypergroups [57784615, 57780515, 57776419, 41007395, 24230179] */
/*Representant 57784615 has neighborhood hypergroups [57784611, 57780519, 57776423, 41007399, 24230183] */
/*Find neighborhood from missing classes */

    
 /*   /*NEW ENUMERATION ALGORITHM
   
   almost working....
   There exist hypergroups which has no neighborhood at distance one.  */ 
let now = Instant::now();
let cardinality = 2u64;
let total = get_min_max_u1024(&cardinality).1;
println!("min hs = {}",total);

/* let (classes,enumeration)= exploring_tree(&total, &cardinality);
let leafs2=collect_leafs_from_nodes(&enumeration,&cardinality);
println!("leafs {:?}",leafs2); */
// OK let new_nodes = find_nodes_not_enumerated(&leafs2, &enumeration, &1usize,&cardinality);
let new_nodes = OTHERS_DIST_2_FROM_LEAFS_3.iter().map(|x|U1024::from(*x)).collect_vec();
let new_classes = exploring_tree(&total, &cardinality);
/*
let new_classes:Vec<_>= new_nodes.par_iter().map(|x|exploring_tree(x, &cardinality)).collect();
 */
println!("new nodes {:?}",new_classes);
println!("new nodes {:?}",new_classes.1.len());
let classes = new_classes.0;
    let permut_vec:Vec<Vec<usize>> = (0..cardinality as usize).permutations(cardinality as usize).collect();
    let mut c:Vec<usize>=Vec::new();
    let mut c_k:Vec<(U1024,Vec<U1024>)>;
    let mut s = String::new();

    for k in 1..=permut_vec.len(){
        c_k=classes.iter().filter(|y|(*y.1).len()==k).into_iter().map(|x|x.clone()).collect_vec();
        c.push(c_k.len());
       
    }
   
    println!("c {:?}",c);
    
 */
/* /* GET TAGS HG AT DIST " FROM LEAFS" */
let mut s = String::new();
let mut  add_str=String::new();
    let leaf = LEAF_TAGS_3.iter().map(|x|U1024::from(*x)).collect_vec();
    let enumeration =TRY_ENUMERATION_3.iter().map(|x|U1024::from(*x)).collect_vec();

    let others :Vec<_>=leaf.par_iter().filter(|tag|
    circumference_radius_d_filtered(&U1024::from(**tag), &2usize, &cardinality).iter().any(|x|
    !enumeration.contains(&U1024::from(x)))).collect();

    s.clear();
    add_str.clear();
    add_str=format!("{:?}\n",others);
    s.push_str(&add_str);

    let _ = write(s.clone(),&format!("try_others_3_dist_2_from_leafs"));

    println!("others {:?}",others);
 */
/* /*COLLECT CLASSES FROM OTHERS */
let cardinality = 3u64;
let enumeration = TRY_ENUMERATION_3.iter().map(|x|U1024::from(*x)).collect_vec();
let others = OTHERS_DIST_2_FROM_LEAFS_3.iter().map(|x|U1024::from(*x)).collect_vec();
let new_hgs = others.iter().map(|x|hg_in_circumference_radius_one(x, &cardinality)).collect_vec();
let mut new_classes:Vec<(U1024,Vec<U1024>)> = others.par_iter().map(|tag|collect_isomorphism_class(tag, &cardinality)).collect();
new_classes.sort_by(|x,y|x.0.cmp(&y.0));
new_classes.dedup();
println!("new classes  {:?}",new_classes.len());
println!("new hgs {:?}",new_hgs); */
   /*      
    classes_dist_1.sort_by(|a, b| a.0.cmp(&b.0));
    classes_dist_1.dedup();
  println!("classes dist 1 {:?}",classes_dist_1);
    let mut visited = circumf_radius_1;
    visited.push(total);
    let mut circumf_radius_2:Vec<Vec<U1024>>=Vec::new();
    let representants: Vec<U1024>= classes_dist_1.iter().map(|x|x.0).collect();
    for tags in representants {
        circumf_radius_1 = hg_in_circumference_radius_one(&tags, &cardinality);
        circumf_radius_1.retain(|x|!visited.contains(x));
        if circumf_radius_1.is_empty(){continue;}
        else{
        circumf_radius_2.push(circumf_radius_1);
        }

    }
    let mut classes_dist_2:Vec<Vec<(U1024, Vec<U1024>)>> =circumf_radius_2.par_iter().map(|x|
        x.iter()
        .map(|tag|
            HyperGroup::new_from_tag_u1024(tag, &cardinality).collect_isomorphism_class()).collect::<Vec<_>>()).collect();
        classes_dist_2.iter_mut().map(|x|x.sort_by(|a,b|a.0.cmp(&b.0))).collect_vec();
        classes_dist_2.iter_mut().map(|x|x.dedup()).collect_vec();
        
        let mut t =classes_dist_2.concat();
        t.sort_by(|a, b| a.0.cmp(&b.0));
        t.dedup();

       println!("classes dist 2 {:?}",t);

 */
    /*     let cardinality = 3u64;
    let max_dist = TAG_HG_3.iter().combinations(2).map(|x| distance_tags(x[0], x[1], &cardinality)).max();
    println!("max dist is {}",max_dist.unwrap());
    let tag_max_dist:Vec<_> = TAG_HG_3.iter().combinations(2).filter(|x| distance_tags(x[0],x[1],&cardinality)==max_dist.unwrap()).collect();
    println!("tag_max {:?}",tag_max_dist); */

   /*  
    let cardinality = 2u64;
let total_hg  =U1024::from(2).pow(U1024::from(cardinality.pow(3)))-1;
    println!("total {}",total_hg);
    let mut visited:Vec<U1024> = Vec::new();
    visited.push(total_hg);

    let mut circumf_radius_1 = circumference_radius_d_filtered(&total_hg, &1usize, &cardinality);
    
    visited.append(&mut circumf_radius_1.clone());
    visited.sort();
    visited.dedup();println!("hg dist 1 {:?}",circumf_radius_1);
    println!("hg dist 1 {}",circumf_radius_1.len());
    let dist =4usize;
    let mut i=1usize;
    loop {
        let x = circumf_radius_1.par_iter().find_first(|x|circumference_radius_d_filtered(&x, &1usize, &cardinality).iter().any(|y|!visited.contains(y)));
        visited.append(&mut circumf_radius_1.clone());
        visited.sort();
        visited.dedup();
        println!("x is {}",x.unwrap());
        println!("dist from  total is {}",distance_tags_u1024(&total_hg, x.unwrap(), &cardinality));

        if i==dist{break println!("hypergroup {} dist from total {}",x.unwrap(),distance_tags_u1024(&total_hg, x.unwrap(), &cardinality));} 
        else {
            
        circumf_radius_1=circumference_radius_d_filtered(&x.unwrap(), &1usize, &cardinality);
        
        //println!("circunf {:?}",circumf_radius_1);
/*         
        circumf_radius_1.retain(|x|!current_circ.contains(x));
        println!("circunf  after retain{:?}",circumf_radius_1); */
        

        println!("visited {:?}",visited);
        }
        i+=1;
    } */ 
/* 
       /*A WAY TO GENERATE HPERGROUPS OF ORDER 6 */
let cardinality =5u64;
let total_hg3 = get_min_max_u1024(&cardinality).1;
let width = cardinality.pow(3u32);
let dist = 7;
let comb = (0..width).into_iter().rev().combinations(dist);/*Don't collect vector, to much memory. It will abort! */
for c in comb {
    let x = c.iter().fold(total_hg3, |total_hg3: U1024,x| total_hg3^(U1024::one()<<*x));
    if representing_hypergroupoid_u1024(&x, &cardinality){
    let hs = HyperGroupoid::new_from_tag_u1024(&x, &cardinality);
        if hs.is_hypergroup() {
            let hg  = HyperGroup::new_from_tag_u1024(&x, &cardinality);


            println!("{}",x);
        }
    }
}  */
/* /*TEST IN HG 6 */
/*This are tags of hypergroups of order 6 in the ball of radius 7 centered in the total_hg_7 */
let cardinality = 5u64;
let tag1 = U1024::from_dec_str("1495381495258030357016782942815387647").expect("error");
let tag2 = U1024::from_dec_str("1972423769135917645388022293062483967").expect("error");
let tag3 = U1024::from_dec_str("2315754257701731296489389389206519807").expect("error");
let tag4 = U1024::from_dec_str("2492135162217487135413640673202536447").expect("error");

let hg1= HyperGroup::new_from_tag_u1024(&tag1, &cardinality);
println!("hg1 {}",hg1);
let permut  = (0..cardinality as usize).permutations(cardinality as usize).into_iter().map(|x|Permutation::oneline(x));
let mut isomorphic_to_tag1:Vec<U1024>=permut.into_iter().map(|sigma| hg1.isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag_u1024()).collect();
isomorphic_to_tag1.dedup();
println!("iso_tag1 {:?}",isomorphic_to_tag1);

let hg2= HyperGroup::new_from_tag_u1024(&tag2, &cardinality);
println!("hg1 {}",hg2);
let permut  = (0..cardinality as usize).permutations(cardinality as usize).into_iter().map(|x|Permutation::oneline(x));
let mut isomorphic_to_tag2:Vec<U1024>=permut.into_iter().map(|sigma| hg2.isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag_u1024()).collect();
isomorphic_to_tag2.dedup();
println!("iso_tag1 {:?}",isomorphic_to_tag2.len());

let hg3= HyperGroup::new_from_tag_u1024(&tag2, &cardinality);
println!("hg3 {}",hg3);
let permut  = (0..cardinality as usize).permutations(cardinality as usize).into_iter().map(|x|Permutation::oneline(x));
let mut isomorphic_to_tag3:Vec<U1024>=permut.into_iter().map(|sigma| hg3.isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag_u1024()).collect();
isomorphic_to_tag3.dedup();
println!("iso_tag1 {:?}",isomorphic_to_tag3.len());


let hg4= HyperGroup::new_from_tag_u1024(&tag2, &cardinality);
println!("hg3 {}",hg4);
let permut  = (0..cardinality as usize).permutations(cardinality as usize).into_iter().map(|x|Permutation::oneline(x));
let mut isomorphic_to_tag4:Vec<U1024>=permut.into_iter().map(|sigma| hg4.isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag_u1024()).collect();
isomorphic_to_tag4.dedup();
println!("iso_tag1 {:?}",isomorphic_to_tag4.len());
 */


/* for tag in TAG_HG_3 {
    let one_from_tag:Vec<_> = TAG_HG_3.par_iter().filter(|x|distance_tags(&x, &tag, &cardinality)==1).collect();
println!("one_from_{} = {:?}",tag,one_from_tag);
} */
/*     let tag = get_min_max_u1024(&cardinality).1;
    let d = 20usize;

    let circ = circunference_radius_d_filtered(&tag, &d, &cardinality);
    println!("circ {:?}",circ);  */   



/*     let cardinality = 3u64;
 */  /*   let tag = 33025917u128;
    let hg = HyperGroup::new_from_tag_u128(&tag, &cardinality);
    let reflexive = vec_to_set(&get_subset(&hg.find_reflexive_subhypergroup().unwrap(),&cardinality));
    println!("hg {}",hg);
    println!("ref {:?}",reflexive);
 */
/*     let (min, max)=get_min_max_u1024(&cardinality);
 *//*     for tag in min.to(max).rev() {
        if representing_hypergroupoid_u1024(&tag, &cardinality) {
            let hs  = HyperGroupoid::new_from_tag_u1024(&tag, &cardinality);
            if hs.is_hypergroup() {
                let hg = HyperGroup::new_from_tag_u1024(&tag, &cardinality);
                if hg.is_transposition(){
                    println!("found transposition hg {}",tag);
                }
            }
        }
    } */
 /*   let ex :Vec<u128>= TAG_HG_3.into_par_iter().filter(|x|
    representing_hypergroupoid(x, &cardinality)&&
HyperGroupoid::new_from_tag_u128(x, &cardinality).is_hypergroup()&&
HyperGroup::new_from_tag_u128(x, &cardinality).is_transposition()&&
HyperGroup::new_from_tag_u128(x, &cardinality).find_reflexive_subhypergroup().is_some()).collect();
println!("{:?}",ex);
println!("{:?}",ex.len());
let dist:Vec<_> = ex.iter().cartesian_product(ex.clone()).map(|(x,y)|distance_tags(&x, &y, &cardinality)).collect();

println!("dist {:?}",dist); */

/*    for tag in min.to(max).rev() {
        if representing_hypergroupoid_u1024(&tag, &cardinality) {
            let hs  = HyperGroupoid::new_from_tag_u1024(&tag, &cardinality);
            if hs.is_hypergroup() {
                let hg = HyperGroup::new_from_tag_u1024(&tag, &cardinality);
            if hg.is_transposition()&&hg.find_reflexive_subhypergroup().is_some() {
                println!("hg is {} and {} is reflexive subhypergroup",hg,hg.find_reflexive_subhypergroup().unwrap())
            }
            }
        }
    } */


    /*   let cardinality =4u64;
  let min = 2305843009213693951;
  let some = 4575657221408423935;
  let max = get_min_max(&cardinality).1;
  let dist = distance_tags(&min, &max, &cardinality);
  println!("dist {}",dist);
  println!("log2dist = {}",dist.ilog2());
    let any  =(2305843009213693951 +1..some).into_par_iter().find_first(|x|(representing_hypergroupoid(x, &cardinality))&&(HyperGroupoid::new_from_tag_u128(&*x, &cardinality).is_hypergroup()));    /*Some Tests With Metric */
    let any_hg= HyperGroupoid::new_from_tag_u128(&any.unwrap(), &cardinality);
    println!("tag of any is {}",any_hg.get_integer_tag());
    println!("tag of max is {}",max);

    let dist = distance_tags(&any.unwrap(), &some, &cardinality);

    println!("dist some max is  ={}",dist);
    println!("log2dist = {}",dist.ilog2());
    println!("{}",any_hg.is_hypergroup());
    any_hg.assert_associativity();
    any_hg.is_reproductive();  */
  /*  
  let cardinality =2u64;
    let mut dist :Vec<(u64,Vec<u128>)> = Vec::new();
    for i in 0..cardinality.pow(3) {
        let tag_i_from_first= TAG_HG_2.iter().filter(|x|distance_tags(&TAG_HG_2[0], &x, &cardinality) ==i).map(|x|*x).collect_vec();
        dist.push((i,tag_i_from_first));
    }
    println!("{:?}",dist);
    let min = get_min_max(&cardinality);
    let dist_last = distance_tags(&min.1, &TAG_HG_2[0], &cardinality);
    println!("dist last 2  ={}",dist_last);
    let dist_first = distance_tags(&min.0, &TAG_HG_2[0], &cardinality);
    println!("dist first 2  ={}",dist_first);
    let second = 107u128;
    let second_dist_min=distance_tags(&second, &min.0, &cardinality);
    let second_dist_max=distance_tags(&second, &min.1, &cardinality);
    println!("second dist min = {} second dist max = {}",second_dist_min,second_dist_max);

    let cardinality =3u64;
let mut dist :Vec<(u64,usize)> = Vec::new();
for i in 0..cardinality.pow(3) {
    let tag_i_from_first= TAG_HG_3.iter().filter(|x|distance_tags(&TAG_HG_3[0], &x, &cardinality) ==i).map(|x|*x).collect_vec();
    dist.push((i,tag_i_from_first.len()));
}
println!("{:?}",dist);
let min = get_min_max(&cardinality);
    let dist_first = distance_tags(&min.0, &TAG_HG_3[0], &cardinality);
    let dist_last = distance_tags(&min.1, &TAG_HG_3[0], &cardinality);
    println!("dist last 2  ={}",dist_last);
    println!("dist first 2  ={}",dist_first);
let v1: [u128; 6]=[20749812, 25020900, 31960177, 33006930, 39137148, 74559292];
let v2: [u128;6]=[20758004, 25029092, 31960181, 33006934, 55914364, 91336508];
let dist:Vec<u64>= v1.iter().zip(v2).map(|(x,y)|distance_tags(x, &y, &cardinality)).collect();
println!("dist classes {:?}",dist);  */
    
/*  let cardinality=3u64;
    let nbeta = collect_beta_not_equivalence(&cardinality);
    println!("nbeta : {}",nbeta.len()); */
     
 /*    let tag=TAG_HG_2[3];
    let hg2=HyperGroupoid::new_from_tag_u128(&tag, &2u64);
    println!("core is {:?}",hg2.beta_relation().rel);
    let tag=TAG_HG_3[100];
    let hg2=HyperGroupoid::new_from_tag_u128(&tag, &3u64);
    println!("core is {:?}",hg2.beta_relation());
    println!("beta is equivalence: {}",hg2.beta_relation().is_equivalence()); */

 /*Example 1 Karim ABBASI, Reza AMERI, Yahya TALEBI-ROSTAMI 
 This is not a hypergroup, so we get a false result!   */
    /* let cardinality=  3u64;
    let hs = HyperGroupoid::new_from_matrix(&DMatrix::from_row_slice(cardinality as usize, cardinality as usize, &[1,6,1,6,1,6,1,6,1]));
    println!("Is hypergroup {}",hs.is_hypergroup());
    println!("Hypergroupoid : {}",hs);
let ph :Vec<HashSet<u64>> = hs.collect_ph().iter().map(|x|vec_to_set(&get_subset(x, &cardinality))).collect();
println!("ph is {:?}",ph);
let beta:Vec<(HashSet<u64>,HashSet<u64>)> = hs.beta_relation().rel.iter().map(|(x,y)| (vec_to_set(&get_subset(x, &cardinality)),vec_to_set(&get_subset(y, &cardinality)))).collect();

println!("beta {:?}",beta);
println!("beta is rif: {}",hs.beta_relation().is_reflexive());
println!("beta is symm: {}",hs.beta_relation().is_symmetric());

println!("beta is trans: {}",hs.beta_relation().is_transitive());
let eq_classes=hs.beta_relation().collect_classes();
println!("classes are {:?}",eq_classes);
 */ 
/*  
/*Example 4.2  Pourhaghani, Anvariyen, Davvaz (OK)
    NOTE:
    This is not working as in the paper, even though it is associative. Why?*/
println!("Example 4.2  Pourhaghani, Anvariyen, Davvaz");

let cardinality=4u64;
let hypergroupoid = HyperGroupoid::new_from_matrix(
    &DMatrix::from_row_slice(
        cardinality as usize,
        cardinality as usize,
        &[1,3,5,9,1,3,5,9,1,2,4,8,1,2,4,8]));
println!("is associative {}",hypergroupoid.is_associative());
println!("Hypergroupoid : {}",hypergroupoid);
let ph = hypergroupoid.collect_ph();
println!("ph is {:?}",ph); 
let beta:Vec<(HashSet<u64>,HashSet<u64>)> = hypergroupoid.beta_relation().rel.iter().map(|(x,y)| (vec_to_set(&get_subset(x, &cardinality)),vec_to_set(&get_subset(y, &cardinality)))).collect();
println!("beta {:?}",beta);
println!("beta is equivalence: {}",hypergroupoid.beta_relation().is_equivalence());
let classes = hypergroupoid.beta_relation().collect_classes();
println!("classes {:?}",classes); 
    
 */
/* 
let cardinality=  5u64;
let hs = HyperGroupoid::new_from_matrix(
    &DMatrix::from_row_slice(
        cardinality as usize, 
        cardinality as usize, 
        &[1,11,1,11,11,1,2,1,11,11,1,11,5,11,31,1,11,1,11,11,1,11,5,11,31]));
println!("Hypergroupoid : {}",hs);
assert!(hs.assert_associativity());
let ph = hs.collect_ph();
println!("ph is {:?}",ph);
let beta = hs.collect_beta_classes();
println!("beta {:?}",beta);
 */    /* 
let cardinality = 2;
let t=185;
println!("{:b}",t);
let new_hyperstructure_from_tag = HyperGroupoid::new_from_tag_u128(&t,&cardinality);
let new_hyperstructure_from_matrix = HyperGroupoid::new_from_matrix(&DMatrix::from_row_slice(2usize,2usize,&[2,3,2,1]));
println!("{}", new_hyperstructure_from_tag);
let t= from_tag_to_vec(&t, &cardinality);

println!("{:?}",t);
println!("{}",new_hyperstructure_from_matrix);
println!("tag1 {}, tag2 {}", new_hyperstructure_from_matrix.get_integer_tag(),185); */
/*     let mat=DMatrix::from_row_slice(3, 3, &[1,2,4,2,5,7,4,2,1]);
    let h=HyperGroupoid::new_from_matrix(&mat);
    let magma=UnitalMagma{
        h:h,
        identity:1
    };

println!("unital magma: {}",magma.is_unital_magma());
let t  =24368401u128;
 let cardinality=3u64;
//let magma =UnitalMagma::new_from_tag_u128(&t, &cardinality);
println!("magma {}",magma);
 println!("magma is invertible : {}",magma.is_invertible_unital_magma());
 let left_invertible= magma.collect_left_invertible();
 let left_inverses:Vec<(u64,Vec<u64>)>=left_invertible.iter().map(|x|(*x,magma.collect_left_inverses(x))).collect();

 let right_invertible=magma.collect_right_invertible();
 let right_inverses:Vec<(u64,Vec<u64>)>=right_invertible.iter().map(|x|(*x,magma.collect_right_inverses(x))).collect();

 println!("left inverses are {:?}",left_inverses);
 println!("right inverses are {:?}",right_inverses); */
/*
    let t=71663230u128;
    let cardinality=3u64;
    let m=UnitalMagma::new_from_tag:u128(&t, &cardinality);
    let l_invertible=m.collect_left_invertible();
    let r_invertible=m.collect_right_invertible();


    println!("{m}");
    let l_inv_x=m.collect_left_inverses(&1u64);
    println!("left_ inverses of {} are {:?}",1u64.ilog2(),l_inv_x);
    let r_inv_x=m.collect_left_inverses(&1u64);
    println!("right_ inverses of {} are {:?}",1u64.ilog2(),r_inv_x);
    println!("left invertible are {:?}",l_invertible);
    println!("right invertible are {:?}",r_invertible);
    println!("m is invertible {}",m.is_invertible_unital_magma()); */
    

/*         /* COLLECT INVERTIBLE UNITAL MAGMATA (L-MOSAICS)*/ 
        let now = Instant::now();

 let cardinality=3u64;
 let c= enumeration_hyperstructure("hypergroups", &cardinality);
println!("c : {:?}",c);
let end = now.elapsed();
    println!("Computation time:{:?}",end);

  */
 
/* let cardinality=3u64;

for m in TAG_UNITAL_MAGMATA_3{
    let magma=UnitalMagma::new_from_tag_u128(&m, &cardinality);
    if magma.is_invertible_unital_magma() {
        println!("{:b} with identity {:b}",m,magma.identity)
    }
  }
    
 */
    

         /*COLLECTING UNITAL MAGMATA */
/*     
let cardinality=2u64;
let now = Instant::now();
collect_hypergroupoid_with_scalar_identity(&cardinality); 
let end = now.elapsed();
    println!("Computation time:{:?}",end);
 */
      /*ISOMORPHIC HYPERGROUPS */
/*   let cardinality=3u64;
  let now = Instant::now();
let c= enumeration_hyperstructure("L_mosaics", &cardinality);
let end = now.elapsed();
println!("Elapsed:{:?}",end);

println!("c : {:?}",c); */
/* 
let now = Instant::now();
    let e= enumeration_hypergroups(&4u64);
    println!("{:?}",e);
let end = now.elapsed();
println!("Elapsed:{:?}",end);
 */
      /*ISOMORPHIC MAGMATA */
/*   let cardinality=3u64;

let c= enumeration_hyperstructure("unital magmata", &cardinality);
println!("c : {:?}",c); */
/*
let now = Instant::now();
    let e= enumeration_hypergroups(&4u64);
    println!("{:?}",e);
let end = now.elapsed();
println!("Elapsed:{:?}",end);
 */ 
/*    let tag =25830028u128;
    let cardinality=3u64;
    println!("starting tag {}",tag);
    let hypergroup=HyperGroupoid::new_from_tag_u128(&tag,&cardinality);
    println!("new from tag {}",hypergroup);
    println!("tag is hypergroup: {}",hypergroup.is_hypergroup());
    println!("tag {}",hypergroup.get_integer_tag());


for i in 0..hypergroup.n{
    if hypergroup.is_left_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u64.pow(i), &hypergroup.n));
        println!("Left identities {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_right_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u64.pow(i), &hypergroup.n));
        println!("Right identities {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_identity(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u64.pow(i), &hypergroup.n));
        println!("Identity {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_left_scalar(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u64.pow(i), &hypergroup.n));
        println!("Left Scalar {:?}",i_singleton)
    }
}
for i in 0..hypergroup.n{
    if hypergroup.is_right_scalar(&i) {

        let i_singleton=vec_to_set(&get_subset(&2u64.pow(i), &hypergroup.n));
        println!("Right scalar {:?}",i_singleton)
    }
}

 */
/*     let args: Vec<String> = env::args().collect();  
    let number: u64 = match args[1].parse() {
        Ok(n) => {
            n
        },
        Err(_) => {
            eprintln!("error: Argument not an u64");
            return;
        },
    };    
    let now = Instant::now();
        let e= enumeration_hypergroups(&number);
        println!("{:?}",e);
    let end = now.elapsed();
    println!("Elapsed:{:?}",end);
 */
/*     let n=6u64;
    for i in 0..8 {

    let set =get_subset(&i, &3u64);
    println!("{} binary {:?}",i,n_to_binary_vec(&(i as u128), &3u64));
    println!("{} as vec is {:?}",i,set);
    let toset=to_set(&set);
    println!("{} to set is {:?}",i,toset);
    } */
    
/*
let h_groupoid=  HyperGroupoid::new_random_from_cardinality(&n);
println!(" A new Hyper Groupoid : {}",h_groupoid);
println!("H is reproductive: {}",h_groupoid.is_reproductive());
let new_hg=h_groupoid.fix_reproductivity();
println!("A reproductive Hypergroupoid: {}",new_hg);
println!("H is repdocutive: {}",new_hg.is_reproductive());
println!("H is associativity: {}",new_hg.is_associative());

 */
/* 
/* GET HYPERSTRUCTURE FROM MATRIX */

let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
let hypergroup=HyperGroupoid::new_from_matrix(&matrix);
println!("{}",hypergroup);
println!("H is hypergroup: {}",hypergroup.is_hypergroup());
let a = 1u64;
let b =2u64;
let a_right_b=hypergroup.right_division(&a,&b);
println!("a / b = {}",a_right_b);
 */
/* 
/*COLLECT ISOMORPHISM CLASS OF A HYPERGROUPS */ 
let cardinality= 3u64;
let tag=U1024::from(20819966u128);
let hg=HyperGroup::new_from_tag_u1024(&tag, &cardinality);
let iso=hg.collect_isomorphism_class();
println!("iso {:?}",iso);
 */
/* /*TEST NUMBER OF ISOMORPHISM IN TERMS OF PERMUTATIONS */
let mut count_isomorphism:u64=0;
let cardinality =3u64;
let tag  =131071999u128;
let hypergroup = HyperGroup::new_from_tag_u128(&tag, &cardinality);
let permut_vec:Vec<Vec<usize>> = (0..hypergroup.cardinality() as usize).permutations(hypergroup.cardinality() as usize).collect();
let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
for sigma in permutation {
    let isomorphic_hg = hypergroup.isomorphic_hypergroup_from_permutation(&sigma);
    println!("isomorphic image {}",isomorphic_hg);
    count_isomorphism+=1;
    }
println!("number of isomorphism {}",count_isomorphism);
 */
/* 
let now = Instant::now();
//let cardinality = 3;
let tag_2= collect_hypergroups(&CARDINALITY);
let permut_vec:Vec<Vec<usize>> = (0..CARDINALITY as usize).permutations(CARDINALITY as usize).collect();
let permutation:Vec<Permutation> = permut_vec.iter().map(|sigma| Permutation::oneline(sigma.clone())).collect();
let mut classes:Vec<(u64,Vec<u64>)>=Vec::new();

for tag in tag_2 {
    let mut isomorphism_classes:Vec<u64>=Vec::new();

    for sigma in &permutation {        
        let isomorphic_image_tag = HyperGroupoid::new_from_tag_u128(tag, &CARDINALITY).isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag();
        isomorphism_classes.push(isomorphic_image_tag);

    }
    let isomorphism_classes:Vec<u64>=isomorphism_classes.iter().sorted().dedup().map(|x|*x).collect();
    let representant_of_class=isomorphism_classes.iter().min().unwrap();
    
        classes.push((*representant_of_class,isomorphism_classes));

    
}
let classes:Vec<&(u64, Vec<u64>)>=classes.iter().sorted_by(|x,y|x.0.cmp(&y.0)).dedup().collect();
let mut c:Vec<usize>=Vec::new();
    let mut c_k:Vec<&(u64,Vec<u64>)>;
    let mut s = String::new();
    for k in 1..=permut_vec.len(){
        c_k=classes.iter().filter(|y|(*y.1).len()==k).into_iter().map(|x|x.clone()).collect_vec();
        c.push(c_k.len());
        let add_str=format!("{:?}\n",c_k);
        s.push_str(&add_str);
        write(s.clone());
        

    }
/* for tag in classes.iter().map(|x|x.0) {
    let hg  = HyperGroupoid::new_from_tag_u128(tag as u128, &CARDINALITY);
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


let n = 5u64;
let mut semihypergroup=   HyperGroupoid::new_random_from_cardinality(&n);

semihypergroup=loop {
    let x=HyperGroupoid::new_random_from_cardinality(&n);

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
