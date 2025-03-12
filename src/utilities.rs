
use std::{collections::HashSet, fs::File, io::Write};
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
/// Converts a Vec<u32> into a HashSet<u32>
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// let v=vec![1,2,0];
/// let v_as_set=vec_to_set(&v);
/// println!("{:?}",v_as_set);
///
pub fn vec_to_set(v:&Vec<u32>)->HashSet<u32>{
    v.iter().map(|x|*x).collect()
}
pub fn power_set (n:&u32)->Vec<Vec<u32>>{
    (0..2u32.pow(*n)).into_iter().map(|i|  get_subset(&i, &n)).collect()
  }
/// Represents the integer $k$ as a subset of the set H={0,1,..,n-1}.
/// There are 2^n different integer representing subsets of H. It will panic if $k is greater then 2^n.
/// The subset S represented by k is given by the binary representation of k: 
/// i is in S if and only if the i-th digit of binary representation of k is a one.
/// Output is Vec<u32>. Use fn vec_to_set to convert it into a HashSet<u32>.
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// use hyperstruc::utilities::get_subset;
/// use std::collections::HashSet;
/// 
/// let cardinality = 4u32;
/// let k=6;
/// let subset_as_vec=get_subset(&k,&cardinality);
/// let subset= vec_to_set(&subset_as_vec);
/// 
/// println!("{:?}",subset);
/// let test_subset:HashSet<u32>= (1..=2).into_iter().collect();
/// assert_eq!(subset,test_subset);
///
/// let cardinality = 3u32;
/// let k=2;
/// let subset_as_vec=get_subset(&k,&cardinality);
/// let subset= vec_to_set(&subset_as_vec);
/// 
/// println!("{:?}",subset);
/// let test_subset:HashSet<u32>= [1].into_iter().collect();
/// assert_eq!(subset,test_subset);
pub fn get_subset(k:&u32,cardinality:&u32)->Vec<u32>{
    let mut subset: Vec<u32> = Vec::new();
    if k>=&2u32.pow(*cardinality){panic!("k can't be grater then 2^n");}
    for i in 0..*cardinality {
        if (k >> i)&1==1{
            subset.push(i);
            }
    }
    subset.iter().map(|x|*x).collect()
}
pub fn ones_positions(k:u32,n:&usize)->Vec<u32>{
    (0..*n as u32).into_iter().filter(|x|(k>>x)&1==1).collect()

}

pub fn cartesian_product(set: &Vec<u32>) -> Vec<(u32, u32)> {
    let mut product: Vec<(u32, u32)> = Vec::new();
    for &a in set {
        for &b in set {
            product.push((a, b));
        }
    }
    product
}
/// Represents a subset $S$ as an integer in H={0,1,..,n-1}.
/// There are 2^n different subsets of H. It will panic if $k is greater then 2^n.
/// Elements in S correspond to indices of occurrences of ones for the binary representation of the integer k which identifies S.
/// Therefore, k equals the sum of power two with exponents in S. 
/// # Example
/// ```
/// use hyperstruc::utilities::vec_to_set;
/// use hyperstruc::utilities::subset_as_u32;
/// use hyperstruc::utilities::get_subset;
/// use std::collections::HashSet;
/// 
/// let subset:HashSet<u32>= (1..=2).into_iter().collect();
/// let subset_as_integer= subset_as_u32(&subset);
/// let test_integer=6u32;
/// assert_eq!(subset_as_integer,test_integer);
/// 
/// let cardinality=4;
/// let k=8;
/// let subset=vec_to_set(&get_subset(&k,&cardinality));
/// let subset_as_integer= subset_as_u32(&subset);
/// assert_eq!(subset_as_integer,k);
/// 
///
pub fn subset_as_u32(k:&HashSet<u32>)->u32{
    k.iter().map(|x|2u32.pow(*x)).sum()
}

pub fn permutaton_matrix_from_permutation(n:&u32,sigma:&Permutation)->DMatrix<u32>{
    let identity: DMatrix<u32>=DMatrix::identity(*n as usize,*n as usize);
    let rows:Vec<Vec<u32>> = identity.row_iter().map(|x|x.iter().map(|z|*z).collect()).collect();

    let x: Vec<u32> =sigma.clone().inverse().apply_slice(rows).concat();
    DMatrix::from_row_slice(*n as usize, *n as usize, &x).transpose()
}
/// Gets binary representation of an integer k with width-numbers-of-bit.
/// # Example
/// ```
/// use hyperstruc::utilities::n_to_binary_vec;/// 
/// let cardinality=4;
/// let k=6;
/// let binary=n_to_binary_vec(&k,&cardinality);
/// let v=vec![0,1,1,0];
/// assert_eq!(v,binary);
///
pub fn n_to_binary_vec(k: &u128, width: &u32) -> Vec<u32> {
	format!("{:0>width$}", format!("{:b}", k), width = *width as usize)
		.chars()
		.map(|x| if x == '1' { 1u32 } else { 0u32 })
		.collect()
}
pub fn binary_to_u32(binary_vec:&Vec<u32>)->u32 {
    let s= binary_vec.iter().map(|x|x.to_string()).collect::<String>();
    u32::from_str_radix(&s, 2).unwrap()

}
pub fn representation_permutation_subset (k:&u128,sigma:&Permutation)->u32 {
    /*
    The input value is k in (0..2^n), therefore it represent a subset of H with |H|=n.
    Any occurrence of 1 in the binary representation of k correspond to an element in the subset S corresponding to k 
    Example: k=5="101"-> S={2,0}. 
    The input value sigma is a permutation of S_n. We build the corresponding permutation matrix and we make it act on the binary representation of k. 
    Then, we convert the result into u32.
    
    We prefer to inverse and normalize sigma
    */
    let binary_k=n_to_binary_vec(&k,&(sigma.len() as u32)).iter().rev().map(|x|*x).collect_vec();
    let x =sigma.apply_slice(&binary_k).iter().rev().map(|x|*x).collect_vec();
    
    binary_to_u32(&x)
}
pub fn representing_hypergroupoid(n:&mut u128,cardinality:&u32)->bool{
    let mut hypergroupoid = true;
    while *n>=2u128.pow(*cardinality) {
        if n.trailing_zeros()>=*cardinality {
            hypergroupoid=false;
            return hypergroupoid;
        }
        *n>>=cardinality;
    }
    hypergroupoid
}
pub fn collect_n_digits(width: &u32,m_representation_hypergroupoid:&u128)->Vec<u32>{
    let binary_n = n_to_binary_vec(m_representation_hypergroupoid, &width);
    let out :Vec<u32>=(0..*width)
        .rev()
        .into_iter()
        .map(|i|
            binary_n.iter().rev().nth(i as usize).unwrap()).
            into_iter()
            .map(
                |x| if *x == 1 { 1u32 } 
                    else { 0u32 })
            .collect();
    out
}
pub fn from_tag_to_vec(tag:&u128, n:&u32)->Vec<Vec<u32>>{
    let mut tag = tag.clone();
    let mut tag_vec:Vec<Vec<u32>>=Vec::new();
    let power_set_cardinality = 2u128.pow(*n);

    loop{
        if tag.trailing_zeros()>=*n{
            println!("{} zeroes found",tag.trailing_zeros());
            panic!("Can't be an hypergroupoid");
        } else {
            tag_vec.push(collect_n_digits(n, &tag));
            tag>>=*n as u128;
        }
        if tag<power_set_cardinality {
            tag_vec.push(collect_n_digits(n, &tag));

            break;}
    }
    tag_vec.iter().map(|x|x.clone()).rev().collect()
    

}
pub fn write(s:String,name:&str)-> std::io::Result<()> {
    let file_name=format!("{}.txt",name);

    let mut file = File::create(file_name)?;
    
    file.write(&s.as_bytes())?;
    Ok(())
}
pub fn get_min_max(cardinality:&u32)->(u128,u128){
    let size= cardinality.pow(3);
    let square= cardinality.pow(2);
let min: u128 = (0..square).into_iter().fold(0, |acc: u128,x: u32|acc+2u128.pow(*cardinality*x));

let max = 2u128.pow(size)-1;
(min,max)

 }

