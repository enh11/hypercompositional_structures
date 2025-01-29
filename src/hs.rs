use core::panic;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::{collections::HashSet, fmt::Display, vec};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use rand:: Rng;
use crate::utilities::{binary_to_u32, cartesian_product, from_tag_to_vec, get_subset, n_to_binary_vec, ones_positions, permutaton_matrix_from_permutation, representation_permutation_subset, representing_hypergroupoid, subset_as_u32, vec_to_set};
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroupoidMat{
    pub h:HashSet<u32>,
    pub hyper_composition:DMatrix<u32>,
    pub n:u32
}
impl HyperGroupoidMat {
/// Generate a new random hyperstructure with cardinality n.
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoidMat;
/// 
/// let n  =4u32;
/// let hyperstructure = HyperGroupoidMat::new_random_from_cardinality(&n);
/// println!("{hyperstructure}");
///  
   pub fn new_random_from_cardinality(n:&u32)->Self{
    let h_vec=(0..*n as u32).into_iter().collect();
        let ht = get_random_hypercomposition_matrix(&n);
        HyperGroupoidMat { 
            h: h_vec, 
            hyper_composition: ht, 
            n: *n}
   }
/// Generate a new hyperstructure given a square matrix. Every entry in the matrix are u32 and represent a subset of H={0,1,2,...,n},
/// where n is the size of the matrix, i.e., the cardinality of the new hyperstructure.
/// In particular, if x,y are elements of H, then x*y is the entries in position (x,y). 
/// For more detail about representation, see pub fn get_subset() in utilities.rs.
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoidMat;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoidMat::new_from_matrix(&matrix);
/// println!("{hyperstructure}");
///
   pub fn new_from_matrix(matrix:&DMatrix<u32>)->Self{
    if !matrix.is_square(){panic!("Matrix must be a square matrix!")}
    let a:Vec<&u32>= matrix.iter().filter(|a|**a==0).collect();
    if !a.is_empty(){panic!("In order to have a hypergroupoid, matrix can't contain zeroes!")}
    let h:HashSet<u32>= (0..matrix.ncols() as u32).into_iter().collect();
    let n=matrix.ncols();
    HyperGroupoidMat{
        h:h,
        hyper_composition:matrix.clone(),
        n:n as u32
   }
}
pub fn new_from_tag(mut tag:&u128,cardinality:&u32)->Self{
    let vector_of_subsets_as_integers=from_tag_to_vec(&mut tag, &cardinality);
    let vector_of_subsets_as_integers: Vec<u32>=vector_of_subsets_as_integers.iter().rev().map(|x|binary_to_u32(x)).collect();

    let hyper_composition_matrix = DMatrix::from_row_slice(*cardinality as usize, *cardinality as usize, &vector_of_subsets_as_integers);
        HyperGroupoidMat::new_from_matrix(&hyper_composition_matrix)
}
pub fn get_integer_tag(&self)->u32{

    binary_to_u32(&self.hyper_composition
        .transpose()//transpose is required because iteration in matrix is by column
        .iter()
        .map(|x|n_to_binary_vec(&(*x as u128), &self.n))
        .rev()
        .concat())

}
pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
    let permutation_hypergroupoids = &self.hyper_composition;
    let alfa =DMatrix::from_iterator(self.n as usize, self.n as usize, 
        permutation_hypergroupoids.iter()
            .map(|x| representation_permutation_subset(&(*x as u128),&sigma)));
    
    HyperGroupoidMat { 
        h: self.h.clone(), 
        hyper_composition:alfa, 
        n: self.n}
}
pub fn isomorphic_hypergroup_from_permutation(&self, sigma:&Permutation)->Self{
    let perm_mat = permutaton_matrix_from_permutation(&self.n, &sigma.clone());
    let isomorphic_matrix=perm_mat.clone()*self.permutation_of_table(sigma).hyper_composition*perm_mat.transpose();
    HyperGroupoidMat::new_from_matrix(&isomorphic_matrix)
}
pub fn is_hypergroup(&self)->bool{
    self.is_associative()&&self.is_reproductive()

}
pub fn is_commutative(&self)->bool{
    for a in self.get_singleton().iter(){
        for b in self.get_singleton().iter(){
            let ab=self.mul_by_representation(a, b);
            let ba=self.mul_by_representation(b, a);
            if ab!=ba {return false;}
        }
    }
    true
}
pub fn is_left_identity(&self,e:&u32)->bool{
    if e>=&self.n {panic!("Not an element in hypergroupoid!")}
    let row_e=self.hyper_composition.row(*e as usize);
    (0..self.n).into_iter().zip(row_e.iter()).all(|x|x.1>>x.0&1==1)
}
pub fn is_right_identity(&self,e:&u32)->bool{
    if e>=&self.n {panic!("Not an element in hypergroupoid!")}
    let col_e=self.hyper_composition.column(*e as usize);
    (0..self.n).into_iter().zip(col_e.iter()).all(|x|x.1>>x.0&1==1)}
pub fn is_identity(&self,e:&u32)->bool{
    self.is_left_identity(&e)&&self.is_right_identity(&e)
}
pub fn is_left_scalar(&self,s:&u32)->bool{
    if s>=&self.n {panic!("Not an element in hypergroupoid!")}
    self.hyper_composition.column(*s as usize).iter().all(|x|x.is_power_of_two())
}
pub fn is_right_scalar(&self,s:&u32)->bool{
    if s>=&self.n {panic!("Not an element in hypergroupoid!")}
    self.hyper_composition.row(*s as usize).iter().all(|x|x.is_power_of_two())
}
pub fn collect_left_identity(&self)->Vec<u32>{
    self.get_singleton().iter()
        .filter(|e|self.is_left_identity(e))
        .map(|e|*e)
        .collect_vec()

}
pub fn collect_scalars(&self)->Vec<u32>{
    self.get_singleton().iter()
        .filter(|s|self.is_left_scalar(&s)&&self.is_right_scalar(&s))
        .map(|x|*x)
        .collect::<Vec<u32>>()
}
pub fn collect_scalar_identity(&self)->Vec<u32>{
    self.collect_scalars()
        .par_iter()
        .filter(|s|self.is_identity(s))
        .map(|x|*x)
        .collect::<Vec<u32>>()
}
/// Represents the integer $k$ as a subset of the set H={0,1,..,n-1}.
/// There are 2^n different integer representing subsets of H. It will panic if $k is greater then 2^n.
/// The subset S represented by k is given by the binary representation of k: 
/// i is in S if and only if the i-th digit of binary representation of k is a one.
/// Output is Vec<u32>. Use fn vec_to_set to convert it into a HashSet<u32>.
/// Reverse function is subset_as_u32() in utilities.rs
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoidMat;
/// use std::collections::HashSet;
/// 
/// let cardinality = 4u32;
/// let hyperstructure = HyperGroupoidMat::new_random_from_cardinality(&cardinality);
/// let k=6;
/// let subset=hyperstructure.get_subset_from_k(&k);
/// println!("{:?}",subset);
/// let test_subset:HashSet<u32>= (1..=2).into_iter().collect();
/// assert_eq!(subset,test_subset);
/// 
/// let k=8;
/// let subset=hyperstructure.get_subset_from_k(&k);
/// println!("{:?}",subset);
/// let test_subset:HashSet<u32>= vec![3].into_iter().collect();
/// assert_eq!(subset,test_subset);
///
pub fn get_subset_from_k(&self,k:&u32)->HashSet<u32>{
    let cardinality = self.h.len() as u32;
    let subset: Vec<u32> = get_subset(&k, &cardinality);
    vec_to_set(&subset)
}
   pub fn mul_by_representation(&self,int_k:&u32,int_l:&u32)->u32{
    let ones_k=ones_positions(*int_k, &self.h.len());
    let ones_l= ones_positions(*int_l, &self.h.len());
    let mut indexes:Vec<(u32,u32)>=Vec::new();
    for a in &ones_k{
        for b in &ones_l{
                indexes.push((*a,*b));
        }
    }
    indexes.iter().fold(0u32, |acc, x| acc|self.hyper_composition[(x.0 as usize,x.1 as usize)])
}
pub fn mul(&self,subset_k:&HashSet<u32>,subset_l:&HashSet<u32>)->u32{
    if !subset_k.is_subset(&self.h)||!subset_l.is_subset(&self.h) { panic!("K and L must be a subsets of H!")};
    let int_k=subset_as_u32(&subset_k);
    let int_l=subset_as_u32(&subset_l);
self.mul_by_representation(&int_k, &int_l)   
}
pub fn left_division(&self,a:&u32,b:&u32)->u32{    
    let sub_a=vec_to_set(&get_subset(&2u32.pow(*a), &self.n));
    let sub_b=2u32.pow(*b);
    self.get_singleton().iter()
    .filter(
        |x| sub_a.is_subset(
            &vec_to_set(&get_subset(
                        &self.mul_by_representation(&sub_b, x), &self.n)
                    )
                )
            ).fold(0, |acc,t|acc|t)

   
}
pub fn right_division(&self,a:&u32,b:&u32)->u32{
        /*This function compute the value a/b={x in H s.t. a in xb} */
        let sub_a=vec_to_set(&get_subset(&2u32.pow(*a), &self.n));
    let sub_b=2u32.pow(*b);
    self.get_singleton().iter()
    .filter(
        |x| sub_a.is_subset(
            &vec_to_set(&get_subset(
                        &self.mul_by_representation(x,&sub_b), &self.n)
                    )
                )
            ).fold(0, |acc,t|acc|t)

   

}
   pub fn is_reproductive(&self)->bool{
    let h:Vec<u32>=Vec::from_iter(0..self.n).iter().map(|_|2u32.pow(self.n)-1).collect();
    /*xH is row_sum */
    let row_sum:Vec<u32> = self.hyper_composition.row_iter().map(|x|x.iter().fold(0u32, |acc,element|acc|element)).collect();
    /*Hx is column sum */
    let col_sum:Vec<u32> = self.hyper_composition.column_iter().map(|x|x.iter().fold(0u32, |acc,element|acc|element)).collect();
    if h==row_sum&&h==col_sum {
        true
    }
    else {
        false
    }
   }
pub fn assert_associativity(&self)->bool{
    for a in &self.get_singleton(){
        for b in &self.get_singleton(){
            for c in &self.get_singleton(){
                let ab_c=self.mul_by_representation(
                    &self.mul_by_representation(&a, &b),&c);
                let a_bc = self.mul_by_representation(&a, &self.mul_by_representation(&b, &c));
                assert_eq!(a_bc,ab_c)
            }
        }
    }
true
}
pub fn is_associative(&self)->bool{
    for a in &self.get_singleton(){
        for b in &self.get_singleton(){
            for c in &self.get_singleton(){
                let ab_c=self.mul_by_representation(
                    &self.mul_by_representation(&a, &b),&c);
                let a_bc = self.mul_by_representation(&a, &self.mul_by_representation(&b, &c));
                if a_bc==ab_c {continue;} else {
                        return false;
                    }
            }
        }
    }
true
}
pub fn get_singleton(&self)->DMatrix<u32>{
    DMatrix::from_row_iterator(1, self.n as usize, (0..self.n).into_iter().map(|i|2u32.pow(i)))
}
}
impl Display for HyperGroupoidMat{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table:DMatrix<String>=DMatrix::from_iterator(self.n as usize, self.n as usize, 
            self.hyper_composition.iter().map(|x|format!("{:?}",vec_to_set(&get_subset(x, &self.n)))));
        
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {}Size:{}\n", self.h, table, self.hyper_composition, self.n )
    }
}


pub fn get_random_hypercomposition_table(n:&u32)->Vec<Vec<u32>>{
    let vec: Vec<u32>=(0u32..*n as u32).into_iter().map(|x|x).collect();
    let index_cartesian=cartesian_product(&vec);
    let mut rng = rand::thread_rng();
    let mut hypercomposition_table=vec![vec![0u32;*n as usize];*n as usize];
    
    for item in index_cartesian {
        hypercomposition_table[item.0 as usize][item.1 as usize]=rng.gen_range(1..2u32.pow(*n as u32))
} 
hypercomposition_table
}
pub fn get_random_hypercomposition_matrix(n:&u32)->DMatrix<u32>{
    let mut rng = rand::thread_rng();
    let m  =DMatrix::from_iterator(*n as usize, *n as usize, (0..n.pow(2)).into_iter().map(|_|rng.gen_range(1..2u32.pow(*n as u32))));
    m
} 

