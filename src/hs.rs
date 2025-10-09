use core::panic;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{collections::HashSet, fmt::Display, vec};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use rand::Rng;
use crate::{binary_relations::relations::Relation, fuzzy::FuzzySubset, utilities::{all_triplets, binary_to_n, from_tag_to_vec, from_tag_u1024_to_vec, get_min_max_u1024, get_subset, n_to_binary_vec, permutaton_matrix_from_permutation, representation_permutation_subset, representing_hypergroupoid_u1024, subset_as_u64, support, u64_to_set, vec_to_set, U1024}};
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroupoid{
    pub h:HashSet<u64>,
    pub hyper_composition:DMatrix<u64>,
    pub n:u64
}
impl HyperGroupoid {
/// Generate a new random hyperstructure with cardinality n.
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// 
/// let n  =4u64;
/// let hyperstructure = HyperGroupoid::new_random_from_cardinality(&n);
/// println!("{hyperstructure}");
///  
   pub fn new_random_from_cardinality(n:&u64)->Self{
    let h_vec=(0..*n as u64).into_iter().collect();
        let ht = get_random_hypercomposition_matrix(&n);
        HyperGroupoid { 
            h: h_vec, 
            hyper_composition: ht, 
            n: *n}
   }
/// 
/// 
/// Generates a new hyperstructure from a binary function Fn(u64, u64) -> u64.
///
/// The function is applied to all pairs `(a, b)` in the Cartesian product `H × H`.
///
/// For each pair:
/// - The result of `f(a, b)` is converted to a binary string.
/// - All binary strings are concatenated into one large binary string.
/// - A tag is computed from this final string to define the new hyperstructure.
///
/// ### Note
/// If you want to treat elements 'a' and 'b' as **singleton sets** (i.e., `{a}` and `{b}`),
/// you need to encode them using their power-of-two representation: `2^a` and `2^b`.
/// This can be computed efficiently with bit shifting:
/// 
/// ```rust
/// let a =0u64;
/// let b = 1u64;
/// let singleton_a = 1 << a;
/// let singleton_b = 1 << b;
/// ```
/// Use these in your function if it expects `a` and `b` as sets rather than elements of H.
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use hyperstruc::hypergroups::HyperGroup;
/// 
/// let cardinality =3u64;
/// let function = |a:u64,b:u64| 1<<a|1<<b;
/// let hs  = HyperGroupoid::new_from_function(function, &cardinality);
/// println!("hs {}",hs);
/// assert!(hs.is_hypergroup());
/// let hg = HyperGroup::new_from_tag_u1024(&hs.get_integer_tag_u1024(), &cardinality);
/// println!("is transpositional {}",hg.is_transposition());
/// assert!(hg.is_transposition());
/// 
///
pub fn new_from_function<F>(function:F,cardinality:&u64)->Self
    where F: Fn(u64,u64) -> u64{
        let tag_vec = (0..*cardinality)
            .into_iter()
            .cartesian_product(0..*cardinality)
                .into_iter()
                .map(|(a,b)|
                    n_to_binary_vec(&(function(a,b) as u128),cardinality)
                ).concat();
        let tag = U1024::from_binary_vec(&tag_vec);
        HyperGroupoid::new_from_tag_u1024(&tag, cardinality)
}

/// 
/// 
/// Generate a new hyperstructure given a square matrix. Every entry in the matrix is u64 and represent a subset of `H = {0,1,2,...,n}`,
/// where `n` is the size of the matrix, i.e., the cardinality of the new hyperstructure.
/// In particular, if `a`,`b` are elements of `H`, then `ab` is the entry in position `(a,b)`. 
/// For more detail about representation, see pub fn get_subset() in utilities.rs.
/// 
/// # Example
/// ```
/// 
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///     &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// println!("{hyperstructure}");
/// 
/// 
   pub fn new_from_matrix(matrix:&DMatrix<u64>)->Self{
    if !matrix.is_square(){panic!("Matrix must be a square matrix!")}
    if matrix.iter().any(|a|*a==0){panic!("In order to have a hypergroupoid, matrix can't contain zeroes!")}
    let h:HashSet<u64>= (0..matrix.ncols() as u64).into_iter().collect();
    let n=matrix.ncols();
    HyperGroupoid{
        h,
        hyper_composition:matrix.clone(),
        n:n as u64
   }
}
/// Constructs a `HyperGroupoid` from a flattened Cayley table.
///
/// This function takes a vector of vectors representing the hyperoperation table, 
/// where each inner vector corresponds to the image of a pair `(a, b) ∈ H × H` under 
/// the hyperoperation. The index in the flat vector is computed as `cardinality * a + b`, 
/// assuming row-major order. Each subset is encoded as a bitmask via the sum of `1 << x` 
/// for all `x ∈ H`.
///
/// # Arguments
/// - `input_array`: A `Vec<Vec<u64>>` of length `n^2` (where `n` is the cardinality),
///   representing the full Cayley table of the hyperoperation.
/// - `cardinality`: The size `n` of the underlying set `H`.
///
/// # Returns
/// A `HyperGroupoid` instance built from the encoded table.
///
/// # Panics
/// Panics if `input_array.len() != cardinality^2`, i.e., if the table is not square.
///
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// 
/// let cardinality = 3u64;
/// let input_cayley = vec![
///     vec![0, 1, 2], vec![0],     vec![2],
///     vec![0, 2],    vec![1, 2],  vec![2],
///     vec![1],       vec![0, 1, 2], vec![2]];
/// let hypergroupoid = HyperGroupoid::new_from_elements(&input_cayley, &cardinality);
/// 
/// let tag = 120776892u128;
/// let hypergroupoid = HyperGroupoid::new_from_tag_u128(&tag, &cardinality);
/// 
pub fn new_from_elements(input_array: &Vec<Vec<u64>>, cardinality:&u64)->Self{
    if input_array.len() as u64!=cardinality.pow(2u32) {panic!("Cardinality is {} and the input vector length is {}. It should be a square-length vector!",cardinality,input_array.len())}
    let function = |a:u64,b:u64| 
        input_array[(cardinality*a+b) as usize]
            .iter()
            .fold(
                0u64, 
                |acc,x|acc|1<<x);
    HyperGroupoid::new_from_function(function, &cardinality)
    
}
/// Generate a new hyperstructure given a tag and the cardinality of the set H. If n is the cardinality, then tag is a u128 less than or equal to n^3. 
/// Its binary representation, divided in groups of n-bits, provide the table of hyperoperation; each group of n-bits corresponds to a subset of H. 
/// For example, if n=2, then a tag must be less or equal to 2^8. The tag t=185 has binary representation 10111000, divided in groups of 2-bits it is
/// 10-11-10-01. The bits 10 represent the subset {1}, the bits 11 represents {0,1}, the bits 01 represents {0}.
/// With this example it follows that 0*0={1}, 0*1={0,1}, 1*0 = {1} and 1*1 = emptyset.
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let cardinality = 2;
/// let t=185;
/// let new_hyperstructure_from_tag = HyperGroupoid::new_from_tag_u128(&t,&cardinality);
/// let new_hyperstructure_from_matrix = HyperGroupoid::new_from_matrix(
///         &DMatrix::from_row_slice(
///         2usize,
///         2usize,
///         &[2,3,2,1]));
/// assert_eq!(new_hyperstructure_from_tag, new_hyperstructure_from_matrix)
/// 
/// 
pub fn new_from_tag_u128(tag:&u128,cardinality:&u64)->Self{
    let vector_of_subsets_as_integers=from_tag_to_vec(&tag, &cardinality);
    let vector_of_subsets_as_integers: Vec<u64>=vector_of_subsets_as_integers.iter().map(|x|binary_to_n(x)).collect();

    let hyper_composition_matrix = DMatrix::from_row_slice(*cardinality as usize, *cardinality as usize, &vector_of_subsets_as_integers);
        HyperGroupoid::new_from_matrix(&hyper_composition_matrix)
}
/// Generate a new hyperstructure given a tag and the cardinality of the set H. If n is the cardinality, then tag is a u128 less than or equal to n^3. 
/// Its binary representation, divided in groups of n-bits, provide the table of hyperoperation; each group of n-bits corresponds to a subset of H. 
/// For example, if n=2, then a tag must be less or equal to 2^8. The tag t=185 has binary representation 10111000, divided in groups of 2-bits it is
/// 10-11-10-01. The bits 10 represent the subset {1}, the bits 11 represents {0,1}, the bits 01 represents {0}.
/// With this example it follows that 0*0={1}, 0*1={0,1}, 1*0 = {1} and 1*1 = emptyset.
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use hyperstruc::utilities::U1024;
/// 
/// let cardinality = 6u64;
/// let hs = HyperGroupoid::new_random_from_cardinality(&cardinality);
/// let tag_hs = hs.get_integer_tag_u1024();
/// let check_hs = HyperGroupoid::new_from_tag_u1024(&tag_hs, &cardinality);
/// assert_eq!(check_hs,hs)
/// 
/// 
pub fn new_from_tag_u1024(mut tag:&U1024,cardinality:&u64)->Self{
    let vector_of_subsets_as_integers=from_tag_u1024_to_vec(&mut tag, &cardinality).iter().map(|x|binary_to_n(x)).collect_vec();
    assert_eq!(vector_of_subsets_as_integers.len(),cardinality.pow(2u32) as usize,"Can't represent an hypergroupoid!");
    let hyper_composition_matrix = DMatrix::from_row_slice(*cardinality as usize, *cardinality as usize, &vector_of_subsets_as_integers);
        
    HyperGroupoid::new_from_matrix(&hyper_composition_matrix)
}
 /// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let cardinality = 2;
/// let t=185;
/// let new_hyperstructure_from_tag = HyperGroupoid::new_from_tag_u128(&t,&cardinality);
/// let new_hyperstructure_from_matrix = HyperGroupoid::new_from_matrix(
///         &DMatrix::from_row_slice(
///             2usize,
///             2usize,
///             &[2,3,2,1]));
/// assert_eq!(new_hyperstructure_from_tag.get_integer_tag(), new_hyperstructure_from_matrix.get_integer_tag())
/// 
/// 
pub fn get_integer_tag(&self)->u128{

    binary_to_n(&self.hyper_composition
        .transpose()//transpose is required because iteration in matrix is by column
        .iter()
        .map(|x|n_to_binary_vec(&(*x as u128), &(self.n as u64)))
        .concat())

}
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use hyperstruc::utilities::U1024;
/// 
/// let cardinality = 2;
/// let t=185;
/// let new_hyperstructure_from_tag = HyperGroupoid::new_from_tag_u128(&t,&cardinality);
/// let new_hyperstructure_from_matrix = HyperGroupoid::new_from_matrix(
///         &DMatrix::from_row_slice(
///             2usize,
///             2usize,
///             &[2,3,2,1]));
/// let tag1 = new_hyperstructure_from_tag.get_integer_tag_u1024();
/// let tag2 = new_hyperstructure_from_matrix.get_integer_tag_u1024();
/// 
/// assert_eq!(tag1,tag2);
/// assert_eq!(tag1,U1024::from(t));
///
/// 
pub fn get_integer_tag_u1024(&self)->U1024 {
    let binary_vector:Vec<u64> = self.hyper_composition
        .transpose()//transpose is required because iteration in matrix is by column
        .iter()
        .map(|x|n_to_binary_vec(&(*x as u128), &(self.n as u64)))
        .concat()
        .into_iter()
        .collect();
    U1024::from_binary_vec(&binary_vector)    //binary_vector.iter().fold(U1024::from(0u64), |acc,x|acc+U1024::from(2).pow(U1024::from(*x)))
}
pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
    let permutation_hypergroupoids = &self.hyper_composition;
    let alfa =DMatrix::from_iterator(self.n as usize, self.n as usize, 
        permutation_hypergroupoids.iter()
            .map(|x| representation_permutation_subset(&(*x as u128),&sigma).try_into().unwrap()));
    
    HyperGroupoid { 
        h: self.h.clone(), 
        hyper_composition:alfa, 
        n: self.n}
}
pub fn isomorphic_hypergroup_from_permutation(&self, sigma:&Permutation)->Self{
    let perm_mat = permutaton_matrix_from_permutation(&(self.n as u64), &sigma.clone());
    let isomorphic_matrix=perm_mat.clone()*self.permutation_of_table(sigma).hyper_composition*perm_mat.transpose();
    HyperGroupoid::new_from_matrix(&isomorphic_matrix)
}
///
/// Collect isomorphism class of an hypergroup. Elements in the class are obtained by permutation of the set defining the hyperstructure.
/// The representant of the class is chosen to be the smaller among the tags in the class.
/// It returns a tuple (representant, class), where class is a vector of tags.
/// 
pub fn collect_isomorphism_class(&self)->(U1024,Vec<U1024>){
    let cardinality = self.n as usize;
    let permutation_vec:Vec<Vec<usize>> = (0..cardinality).permutations(cardinality ).collect();
    let permutation:Vec<Permutation> = permutation_vec.par_iter()
        .map(|sigma| 
                Permutation::oneline(sigma.clone())
            ).collect();
    let mut isomorphism_classes:Vec<U1024>=permutation.par_iter()
        .map(|sigma|
                self.isomorphic_hypergroup_from_permutation(&sigma).get_integer_tag_u1024()
            ).collect();
   isomorphism_classes.sort();
   isomorphism_classes.dedup();
   let representant_of_class= isomorphism_classes.clone()
        .into_iter()
        .sorted().dedup()
            .into_iter().min().unwrap();
        
       (representant_of_class,isomorphism_classes)

}
///
/// Return true if hyperstructure is both associative and reproductive.
///  
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// 
/// let cardinality = 3u64;
/// let tag = 22150143u128;
/// let hs = HyperGroupoid::new_from_tag_u128(&tag,&cardinality);
/// 
/// assert!(hs.is_hypergroup());
/// 
/// use hyperstruc::utilities::U1024;
/// 
/// let cardinality  =2u64;
/// use hyperstruc::hg_2::tag_hypergroups_2::TAG_HG_2;
/// 
/// let tag_2:Vec<U1024> = TAG_HG_2.iter().map(|x| U1024::from(*x)).collect();
/// for tags in tag_2 {
///     let hg =  HyperGroupoid::new_from_tag_u1024(&tags,&cardinality);
///     assert!(hg.is_hypergroup())
/// }
/// 
/// 
pub fn is_hypergroup(&self)->bool{
    self.is_associative()&self.is_reproductive()
}
pub fn is_weak_associative(&self)->bool{
    self.get_singleton().iter().all(|a|
        self.get_singleton().iter().all(|b|
            self.get_singleton().iter()
                .all(|c| {
                    self.mul_by_representation(
                    &self.mul_by_representation(a, b),
                    c)
                    &
                    self.mul_by_representation(
                    a,
                    &self.mul_by_representation(b, c))!=0
                })  
        )
    )

}
pub fn is_hv_group(&self)->bool {
    self.is_weak_associative()&&self.is_reproductive()
}
pub fn is_isomorphic_to(&self,other: &Self)->bool{
    let total_tag = get_min_max_u1024(&self.n).1;
    let total_hg = &HyperGroupoid::new_from_tag_u1024(&total_tag, &self.n);
    if self.hamming_distance_u1024(total_hg)==other.hamming_distance_u1024(total_hg){
        self.collect_isomorphism_class().1.contains(&other.get_integer_tag_u1024())
    }
    else{

        false
    }
}
/// Checks whether the structure is commutative.
/// 
/// # Returns
///
/// `true` if commutative, `false` otherwise.
///
/// # Examples
/// 
/// Commutative case:
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use hyperstruc::generating_functions::b_hypercomposition;
/// 
/// let cardinality = 7u64;
/// let hs = HyperGroupoid::new_from_function(b_hypercomposition(), &cardinality);
/// assert!(hs.is_commutative());
/// ```
/// 
/// Non-commutative case:
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// 
/// let cardinality = 4u64;
/// let matrix = DMatrix::from_row_slice(
///     4, 4,
///     &[2, 2, 3, 14,
///       5, 14, 3, 3, 
///       1, 11, 12, 7,
///       7, 3, 8, 8]);
/// let hs = HyperGroupoid::new_from_matrix(&matrix);
/// assert!(!hs.is_commutative());
/// ```
/// 
pub fn is_commutative(&self)->bool{
    self.get_singleton().into_iter().combinations(2).into_iter()
        .all(|v|
            self.mul_by_representation(&v[0], &v[1])
            ==
            self.mul_by_representation(&v[1], &v[0]))
}
pub fn is_left_partial_identity_of_x(&self,e:&u64,x:&u64) -> bool{
    if !e.is_power_of_two()|!x.is_power_of_two() {panic!("Some input value is not an element in hypergroupoid!")}
    *x&self.mul_by_representation(e, x)==*x
}
pub fn is_right_partial_identity_of_x(&self,e:&u64,x:&u64) -> bool{
    if !e.is_power_of_two()|!x.is_power_of_two() {panic!("Some input value is not an element in hypergroupoid!")}
    *x&self.mul_by_representation(x, e)==*x
}
pub fn collect_partial_left_identities_of_x(&self,x:&u64)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|e|
            self.is_left_partial_identity_of_x(&e, x))
        .collect()
}
pub fn collect_partial_right_identities_of_x(&self,x:&u64)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|e|
            self.is_right_partial_identity_of_x(&e, x))
        .collect()
}
pub fn collect_partial_left_identities(&self)->Vec<u64>{
    let mut i_pl = self.get_singleton()
        .into_iter()
        .map(|x|
            self.collect_partial_left_identities_of_x(&x))
        .concat();
    i_pl.sort();
    i_pl.dedup();
    i_pl
}
pub fn collect_partial_right_identities(&self)->Vec<u64>{
    let mut i_pr = self.get_singleton()
        .into_iter()
        .map(|x|
            self.collect_partial_right_identities_of_x(&x))
        .concat();
    i_pr.sort();
    i_pr.dedup();
    i_pr
}
pub fn collect_partial_identities(&self)->Vec<u64>{
    let mut i_pl = self.collect_partial_left_identities();
    let mut i_pr = self.collect_partial_right_identities();
    i_pl.append(&mut i_pr);
    i_pl.sort();
    i_pl.dedup();
    i_pl
}
pub fn is_left_identity(&self,e:&u64)->bool{
    if !e.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
    let e=e.trailing_zeros();//the number of trailing_zeros in a power of two integer is equal to the exponent of that power.
    let row_e=self.hyper_composition.row(e as usize);
    (0..self.n).into_iter().all(|x| (row_e.index(x as usize)>>x)&1==1)

}
pub fn is_right_identity(&self,e:&u64)->bool{
    if !e.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
    let e=e.trailing_zeros();//the number of trailing_zeros in a power of two integer is equal to the exponent of that power.
    let col_e=self.hyper_composition.column(e as usize);
    (0..self.n).into_iter().all(|x| (col_e.index(x as usize)>>x)&1==1)
    }
pub fn is_identity(&self,e:&u64)->bool{
    self.is_left_identity(&e)&&self.is_right_identity(&e)
}
pub fn collect_left_identities(&self)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|e|self.is_left_identity(e))
        .collect()
}
pub fn show_left_identities(&self){
    let li = self.collect_left_identities();
    if li.is_empty() {println!("No left identity found.")}
    else {
        let left_id:Vec<String> = li
            .iter()
            .map(|e| 
                format!("{}",e.trailing_zeros())
            ).collect();
        println!("Left identities: {:?}",left_id.join(", "));
        }
    
}
pub fn collect_right_identities(&self)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|e|self.is_right_identity(e))
        .collect()
}
pub fn show_right_identities(&self){
    let ri = self.collect_right_identities();
    if ri.is_empty() {println!("No right identity found.")}
    else {
        let right_id:Vec<String> = ri
            .iter()
            .map(|e| 
                format!("{}",e.trailing_zeros())
            ).collect();
        println!("Right identities: {:?}",right_id.join(", "));
        }
    
}
pub fn collect_identities(&self)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|e|self.is_identity(&e))
        .collect()
}
pub fn show_identities(&self) {
    let identities = self.collect_identities();
    if identities.is_empty() {println!("No identity found.")}
    else {
        let id:Vec<String> = identities
            .iter()
            .map(|e| 
                format!("{}",e.trailing_zeros())
            ).collect();
        println!("Left identities: {:?}",id.join(", "));
        
    }
}
pub fn is_left_scalar(&self,s:&u64)->bool{
    if !s.is_power_of_two() {panic!("Not an element in hypergroupoid!")}
    let s=s.trailing_zeros();
    self.hyper_composition.row(s as usize).iter().all(|x|x.is_power_of_two())
}
pub fn is_right_scalar(&self,s:&u64)->bool{
    if !s.is_power_of_two()  {panic!("Not an element in hypergroupoid!")}
    let s=s.trailing_zeros();
    self.hyper_composition.column(s as usize).iter().all(|x|x.is_power_of_two())
}
pub fn collect_left_scalars(&self)->Vec<u64> {
    self.get_singleton().into_iter().filter(|x|self.is_left_scalar(x)).collect()
}
pub fn collect_right_scalars(&self)-> Vec<u64> {
        self.get_singleton().into_iter().filter(|x|self.is_right_scalar(x)).collect()
}
pub fn is_scalar(&self,s:&u64)->bool {
    self.is_left_scalar(s)&&self.is_right_scalar(s)
}
pub fn collect_scalars(&self)->Vec<u64>{
    self.get_singleton().into_iter()
        .filter(|s|self.is_scalar(&s))
        .collect()
}
pub fn collect_scalar_identities(&self)->Vec<u64>{
    self.collect_scalars()
        .into_iter()
        .filter(|s|self.is_identity(s))
        .collect()
}
/// 
/// Return the subset of right inverses of x with respect to the identity e.
/// Return a u64 representing the subset `x\e-{e}`. 
/// To convert it into a HashSet use `u64_to_set()`.
/// 
pub fn right_inverses_of_x(&self, x:&u64,e:&u64)->u64{
    let h = (1<<self.n)-1;
    assert!(x.is_power_of_two()&&x<&h,"{} is not an element in H",x);
    assert!(self.is_identity(e),"{} is not an identity in H",e);
    let right_inverses=self.left_division(e, x);
    let e_as_element = e.trailing_zeros();
    right_inverses & !(1 << e_as_element)//remove the occurrence of 1 in e^th position (if there is) in order to remove e from the set of inverses.

}
/// 
/// Return the subset of left inverses of x with respect to the identity e.
/// Return a u64 representing the subset `e/x-{e}`. 
/// To convert it into a HashSet use `u64_to_set()`.
/// 
pub fn left_inverses_of_x(&self, x:&u64,e:&u64)->u64{
    let h = (1<<self.n)-1;
    assert!(x.is_power_of_two()&&x<&h,"{} is not an element in H",x);
    assert!(self.is_identity(e),"{} is not an identity in H",e);    
    
    let left_inverses=self.right_division(e, x);
    let e_as_element = e.trailing_zeros();
    left_inverses & !(1 << e_as_element)//remove the occurrence of 1 (if there is) in order to remove u from the set of inverses. 
}
/// 
/// Returns left inverses of `x` in the hyperstructure `H`. Returns a `Vec<(u64,u64)>`, where the first element in the tuple is 
/// the identity `u` with respect to which the inverses are computed, and the second `u64` in the tuple is the integer representation 
/// of the left identities with respect to `u`. To convert them as HashSet use `u64_to_set()`.
/// 
pub fn collect_left_inverses_of_x(&self,x:&u64)->Vec<(u64,u64)>{
    let left_identities_of_x = Vec::new();
    let units = self.collect_identities();
    if units.is_empty() {return left_identities_of_x}
    units.iter().map(|u|(*u,self.left_inverses_of_x(x,u))).collect_vec()
}
pub fn show_letf_inverses_of_x(&self,singleton_x:&HashSet<u64>)->Vec<(HashSet<u64>,HashSet<u64>)> {
    let x = subset_as_u64(singleton_x);
    assert!(x.is_power_of_two(),"Input subset must be a singleton!");
    self.collect_left_inverses_of_x(&x)
        .iter()
        .map(|(identity,inverses)|
            (u64_to_set(identity, &self.n),u64_to_set(inverses, &self.n))).collect()
}
///
/// 
/// Returns right inverses of `x` in the hyperstructure `H`. Returns a `Vec<(u64,u64)>`, where the first element in the tuple is 
/// the unit `u` with respect to which the identities are computed, and the second `u64` in the tuple is the integer representation 
/// of the right identities with respect to `u`. To convert them as HashSet `use u64_to_set()`.
/// 
/// 
pub fn collect_right_inverses_of_x(&self,x:&u64)->Vec<(u64,u64)>{
    let right_identities_of_x = Vec::new();
    let units = self.collect_identities();
    if units.is_empty() {return right_identities_of_x}
    units.iter().map(|u|(*u,self.right_inverses_of_x(x,u))).collect_vec()
}
pub fn show_right_inverses_of_x(&self,singleton_x:&HashSet<u64>)->Vec<(HashSet<u64>,HashSet<u64>)> {
    let x = subset_as_u64(singleton_x);
    assert!(x.is_power_of_two(),"Input subset must be a singleton!");
    self.collect_right_inverses_of_x(&x)
        .iter()
        .map(|(identity,inverses)|
            (u64_to_set(identity, &self.n),u64_to_set(inverses, &self.n))).collect()
}
pub fn collect_inverses_of_x(&self,x:&u64)->Vec<(u64,u64)>{
    self.collect_left_inverses_of_x(x).iter()
    .filter(|inverses|
        self.collect_right_inverses_of_x(x)
        .contains(&inverses))
    .map(|x|*x)
    .collect_vec()
}
pub fn show_inverses_of_x(&self,singleton_x:&HashSet<u64>)->Vec<(HashSet<u64>,HashSet<u64>)> {
    let x = subset_as_u64(singleton_x);
    assert!(x.is_power_of_two(),"Input subset must be a singleton!");
    self.collect_inverses_of_x(&x)
        .iter()
        .map(|(identity,inverses)|
            (u64_to_set(identity, &self.n),u64_to_set(inverses, &self.n))).collect()
}
pub fn collect_all_finite_hyperproducts(&self)->(Vec<u64>,u64){
    let mut a_current = Vec::new();
    let  mut a_next=self.get_singleton();
    let mut stabilization_index = 0u64;
    while a_next!=a_current {
        a_current = a_next;
        a_next=a_current.iter().cartesian_product(a_current.iter()).map(|(x,y)|self.mul_by_representation(x, y)).unique().collect();
        a_next.extend(a_current.iter());
        a_next.sort();
        a_next.dedup();
        stabilization_index+=1;
    }
    (a_current,stabilization_index)
}
/// # Computing β Relation in Hypergroupoid
///
/// This module implements functionality for computing the elements of the
/// **β relation** in a hypergroupoid.
///
/// The β relation is defined over finite hyperproducts of a hypergroupoid \( H \).
///
/// ## Reference
/// This implementation is based on the algorithm proposed by Pourhaghani, Anvariyeh, and Davvaz:
///
/// [An algorithm to calculate the members of the relation β in hypergroupoids](https://www.researchgate.net/publication/386742470_An_algorithm_to_calculate_the_members_of_the_relation_beta_in_hypergroupoids)
///
/// where they proved that 
///
/// β = ⋃_{q ∈ P(H)} q × q
/// 
/// 
/// where:
/// - `P(H)` is the set of all finite hyperproducts of the hypergroupoid `H`,
/// - `q × q` is the Cartesian product of each such set `q`.
///
/// ## Usage
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// 
/// let cardinality=4u64;
/// let hs = HyperGroupoid::new_from_matrix(
///     &DMatrix::from_row_slice(
///     cardinality as usize,
///     cardinality as usize,
///      &[6,10,10,10,
///       10,10,10,10,
///       10,10,10,10,
///       10,10,10,10]));
/// let beta = hs.beta_relation();
/// let expected_beta: Vec<(u64, u64)>  = vec![(0, 0), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (3, 1), (3, 3)];
/// assert_eq!(beta.rel,expected_beta);
///
pub fn beta_relation(&self)->Relation{
    let ph=self.collect_all_finite_hyperproducts();

    let mut beta:Vec<(u64, u64)>=ph.0.iter()
        .map(
            |q|
                support(q, &self.n).into_iter()
                .cartesian_product(
                    support(q, &self.n).into_iter())
                .map(|(x,y)|(x as u64,y as u64))
                .collect()
            ).concat();
    beta.sort();
    beta.dedup();

        Relation{
            a:self.h.clone(),
            b:self.h.clone(),
            rel:beta
        }

}
pub fn collect_beta_classes(&self)->Vec<(u64,Vec<u64>)>{
    self.beta_relation().quotient_set()
}
/// Represents the integer $k$ as a subset of the set H={0,1,..,n-1}.
/// There are 2^n different integer representing subsets of H. It will panic if $k is greater then 2^n.
/// The subset S represented by k is given by the binary representation of k: 
/// i is in S if and only if the i-th digit of binary representation of k is a one.
/// Output is Vec<u64>. Use fn vec_to_set to convert it into a HashSet<u64>.
/// Reverse function is subset_as_u64() in utilities.rs
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use std::collections::HashSet;
/// 
/// let cardinality = 4u64;
/// let hyperstructure = HyperGroupoid::new_random_from_cardinality(&cardinality);
/// let k=6;
/// let subset=hyperstructure.get_subset_from_k(&k);
/// println!("{:?}",subset);
/// let test_subset:HashSet<u64>= (1..=2).into_iter().collect();
/// assert_eq!(subset,test_subset);
/// 
/// let k=8;
/// let subset=hyperstructure.get_subset_from_k(&k);
/// println!("{:?}",subset);
/// let test_subset:HashSet<u64>= vec![3].into_iter().collect();
/// assert_eq!(subset,test_subset);
///
pub fn get_subset_from_k(&self,k:&u64)->HashSet<u64>{
    let cardinality = self.h.len() as u64;
    let subset: Vec<u64> = get_subset(&k, &cardinality);
    vec_to_set(&subset)
}
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// 
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let a =2u64;
/// let b=1u64;
/// let ab=2u64;
/// let mul=hyperstructure.mul_by_representation(&a,&b);
/// assert_eq!(ab,mul);
pub fn mul_by_representation(&self,subset_a:&u64,subset_b:&u64)->u64{
    support(subset_a, &(self.n)).iter()
        .cartesian_product(support(subset_b, &(self.n)))
            .into_iter()
            .fold(
            0u64, 
            |acc,(x,y)|
                acc|self.hyper_composition[(*x ,y)])
}
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let a = [1].into();
/// let b = [0].into();
/// let ab = 2u64;
/// let mul=hyperstructure.mul(&a,&b);
/// assert_eq!(ab,mul);
/// 
pub fn mul(&self,subset_k:&HashSet<u64>,subset_l:&HashSet<u64>)->u64{
    if !subset_k.is_subset(&self.h)||!subset_l.is_subset(&self.h) { panic!("K and L must be a subsets of H!")};
    let int_k=subset_as_u64(&subset_k);
    let int_l=subset_as_u64(&subset_l);
self.mul_by_representation(&int_k, &int_l)   
}
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let a = [1].into();
/// let b = [0].into();
/// let expected_ab = [1].into();
/// let ab=hyperstructure.subsets_multiplication(&a,&b);
/// assert_eq!(ab,expected_ab);
/// 
pub fn subsets_multiplication(&self,set_a:&HashSet<u64>,set_b:&HashSet<u64>)->HashSet<u64> {
    let a = subset_as_u64(&set_a);
    let b = subset_as_u64(&set_b);
    let ab = self.mul_by_representation(&a, &b);
    
    u64_to_set(&ab, &self.n)

}
/// 
/// 
/// Returns `b\a = {x in H : a in bx}`.
/// Input `subset_a` and `subset_b` are `HashSet<u64>`, representing subsets of `H`. 
/// 
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// 
/// let cardinality = 3u64;
/// let input_array=vec![
///     vec![0],     vec![1],     vec![0,1,2],
///     vec![1],     vec![0,1,2], vec![0,1,2],
///     vec![0,1,2], vec![0,1,2], vec![0,2]
///     ];
/// let hyperstructure=HyperGroupoid::new_from_elements(&input_array,&cardinality);
/// let subset_a:HashSet<u64> = [1].into();
/// let subset_b:HashSet<u64> = [2].into();
/// let expected_division:HashSet<u64> = [0,1].into();
/// let a_left_divided_by_b = hyperstructure.subset_left_division(&subset_a,&subset_b);
/// assert_eq!(expected_division,a_left_divided_by_b);
/// 
pub fn subset_left_division(&self,subset_a:&HashSet<u64>,subset_b:&HashSet<u64>)->HashSet<u64> {
    let a = subset_as_u64(subset_a);
    let b = subset_as_u64(subset_b);
    let a_left_divided_by_b = self.left_division(&a, &b);
    
    u64_to_set(&a_left_divided_by_b, &self.n)
}
/// 
/// 
/// Returns `a/b={x in H : a in xb}`.
/// Input `subset_a` and `subset_b` are `HashSet<u64>`, representing subsets of `H`. 
/// 
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// 
/// let cardinality = 3u64;
/// let input_array=vec![
///     vec![0],     vec![1],     vec![0,1,2],
///     vec![1],     vec![0,1,2], vec![0,1,2],
///     vec![0,1,2], vec![0,1,2], vec![0,2]
///     ];
/// let hyperstructure=HyperGroupoid::new_from_elements(&input_array,&cardinality);
/// let a:HashSet<u64> = [0].into();
/// let b:HashSet<u64> = [1].into();
/// let expected_division:HashSet<u64> = [1,2].into();
/// let a_right_divided_by_b=hyperstructure.subset_right_division(&a,&b);
/// assert_eq!(expected_division,a_right_divided_by_b);
/// 
pub fn subset_right_division(&self,subset_a:&HashSet<u64>,subset_b:&HashSet<u64>)->HashSet<u64> {
    let a = subset_as_u64(subset_a);
    let b = subset_as_u64(subset_b);
    let a_right_divided_by_b = self.right_division(&a, &b);
    
    u64_to_set(&a_right_divided_by_b, &self.n)
}
/// 
/// 
/// Compute the integer representation of `b\a = {x in H : a in bx}`.
/// Input `a` and `b` must be type `u64`, representing subsets of `H`. 
/// You can use the function `subset_as_u64()` to convert an `HashSet<u64>` into its unique integer representation.
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let a = 2u64;
/// let b = 4u64;
/// let expected_division = 3u64;
/// let a_left_divided_by_b = hyperstructure.left_division(&a,&b);
/// assert_eq!(expected_division,a_left_divided_by_b);
pub fn left_division(&self,a:&u64,b:&u64)->u64{    
    self.get_singleton().iter()
        .filter(|x| 
            a&(&self.mul_by_representation(&b,&x))!=0 //This is equivalent to {a} meets {b}{x}, i.e., {a} intersection {b}{x} not empty. Therefore this algorithm also compute A\B
        )
        .fold(0, |acc, x| acc|x)  
}
/// 
/// 
/// Returns the integer representation of `a/b={x in H : a in xb}`.
/// Input `a` and `b` must be type `u64`, representing non empty subset of `H`.
/// You can use the function `subset_as_u64()` to convert an `HashSet<u64>` into its unique integer representation.
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// use std::collections::HashSet;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let a = 1u64;
/// let b = 2u64;
/// let expected_division = 6u64;
/// let a_right_divided_by_b=hyperstructure.right_division(&a,&b);
/// assert_eq!(expected_division,a_right_divided_by_b);
/// 
pub fn right_division(&self,a:&u64,b:&u64)->u64{    
    self.get_singleton().iter()
        .filter(|x| 
            a&(&self.mul_by_representation(&x,&b))!=0
        )            
        .fold(0, |acc, x| acc|x)
}

///
/// Return true if hyperstructure is `reproductive`, i.e., `xH = H = Hx` holds for all `x` in `H`.
///
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// println!("{hyperstructure}");
/// assert!(hyperstructure.is_reproductive())
///
pub fn is_reproductive(&self)->bool{
    let h = (1<<self.n)-1;
    self.hyper_composition
        .row_iter()
        .all(|x|
            x.iter().fold(0u64, |acc,element|
                acc|element)==h)
    &&
    self.hyper_composition
        .column_iter()
        .all(|x|
            x.iter().fold(0u64, |acc,element|
                acc|element)==h)
}
/// 
/// 
/// Return true if hyperstructure is associative, i.e., `(xy)z = x(zy)` holds for all `x` in `H`.
///
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// println!("{hyperstructure}");
/// assert!(hyperstructure.assert_associativity())
///
pub fn assert_associativity(&self)->bool{
    all_triplets(self.n)
        .all(
        |(a,b,c)|
                {
                self.mul_by_representation(
                    &self.mul_by_representation(&(1<<a), &(1<<b)),&(1<<c))
                == 
                self.mul_by_representation(&(1<<a), &self.mul_by_representation(&(1<<b), &(1<<c)))
            }
    )
}
/// Return true if hyperstructure is associative, i.e., (xy)z = x(zy) holds for all x in H.
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(
///         3usize,
///         3usize,
///         &[1,2,7,2,7,7,7,7,5]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// println!("{hyperstructure}");
/// assert!(hyperstructure.is_associative())
///
pub fn is_associative(&self)->bool{
    self.get_singleton().iter().all(|a|
        self.get_singleton().iter().all(|b|
            self.get_singleton().iter()
                .all(|c| {
                    self.mul_by_representation(
                    &self.mul_by_representation(a, b),
                    c)
                    ==
                    self.mul_by_representation(
                    a,
                    &self.mul_by_representation(b, c))
                })  
        )
    )
}
pub fn get_singleton(&self)->Vec<u64>{
    (0..self.n).into_iter().map(|i|1<<i).collect()
}
pub fn get_sets_singleton(&self)->Vec<HashSet<u64>>{
    (0..self.n).into_iter().map(|i|[i].into()).collect()
}
/// Calculate the distance between two hyperstructures. The distance is defined as the the 
/// Hamming distance between binary representations of hyperstructure's tags, i.e., 
/// the number of positions at which the corresponding binary tags are different.
/// 
///
pub fn hamming_distance(&self,other:&HyperGroupoid)->usize {
    assert_eq!(self.n,other.n);
    (self.get_integer_tag()^other.get_integer_tag()).count_ones() as usize
}
pub fn hamming_distance_u1024(&self, other:&HyperGroupoid)->usize{
    assert_eq!(self.n, other.n);
    let dist:u32 = (self.get_integer_tag_u1024()^other.get_integer_tag_u1024()).to_little_endian().iter().map(|x|x.count_ones()).sum();
    dist as usize
}
pub fn get_corsini_fuzzysubset(&self)->FuzzySubset{
        FuzzySubset::new(self.clone(), self.get_corsini_membership_function())
    }
}
impl Display for HyperGroupoid{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table:DMatrix<String>=DMatrix::from_iterator(self.n as usize, self.n as usize, 
            self.hyper_composition.iter().map(|x|format!("{:?}",vec_to_set(&get_subset(x, &self.n)))));
        
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {}Size:{}\n", self.h, table, self.hyper_composition, self.n )
    }
}
pub fn get_random_hypercomposition_table(n:&u64)->Vec<Vec<u64>>{
    let hypercomposition_table=vec![vec![0u64;*n as usize];*n as usize];
    let mut rng = rand::thread_rng();
    hypercomposition_table.iter()
        .map(|a|
            a.iter()
            .map(
                |_| rng.gen_range(1u64..1<<n)
            ).collect_vec()
        ).collect_vec()
}
pub fn get_random_hypercomposition_matrix(n:&u64)->DMatrix<u64>{
    let mut rng = rand::thread_rng();
    DMatrix::from_iterator(
        *n as usize,
        *n as usize,
        (0..n.pow(2))
            .into_iter()
            .map(|_| rng.gen_range(1..(1<<n))))
    
} 
pub fn distance_tags(tag1:&u128,tag2:&u128,cardinality:&u64)->u64 {
    let width = cardinality.pow(3);
    let binary_tag1 = n_to_binary_vec(tag1, &width);
    let binary_tag2 = n_to_binary_vec(tag2, &width);

    binary_tag1.iter().zip(binary_tag2).into_iter().filter(|(x,y)|*x!=y).count() as u64
}
pub fn distance_tags_u1024(tag1:&U1024,tag2:&U1024,cardinality:&u64)->usize{
    let width = cardinality.pow(3);
    let binary_tag1 = from_tag_u1024_to_vec(tag1, &width).concat();
    let binary_tag2 = from_tag_u1024_to_vec(tag2, &width).concat();

    binary_tag1.iter().zip(binary_tag2).into_iter().filter(|(x,y)|*x!=y).count()
}
///
/// Collects all binary strings that differ by d bits from the tag binary string.
/// The distance is intended to be the Hamming's distance. 
/// 
/// 
pub fn circumference_radius_d(tag:&U1024,d:&usize,cardinality:&u64)->Vec<U1024>{
    let width = cardinality.pow(3u32);
    let circunference:Vec<_> = (0..width).into_iter().combinations(*d)
    .into_iter()
    .map(|pos|
        pos
            .iter()
            .fold(*tag,
                |acc,x| acc^(U1024::from(U1024::one()<<*x))
            )).collect();

    circunference
}
///
/// Collects all binary strings that differ by d bits from the tag binary string and 
/// filter them to take only those representing hypergroups.
/// The distance is intended to be the Hamming's distance. 
/// 
/// 
pub fn circumference_radius_d_filtered(tag:&U1024,d:&usize,cardinality:&u64)->Vec<U1024>{
    circumference_radius_d(tag, d, cardinality).par_iter().filter(|x|representing_hypergroupoid_u1024(&x, cardinality)&&HyperGroupoid::new_from_tag_u1024(x, cardinality).is_hypergroup())
    .map(|x|*x).collect::<Vec<U1024>>()
}
///
/// Collects all binary strings that differ by 1 bits from the tag binary string and 
/// filter them to take only those representing hypergroups.
/// The distance is intended to be the Hamming's distance. 
/// 
/// 
pub fn hg_in_circumference_radius_one(tag:&U1024,cardinality:&u64)->Vec<U1024>{
    circumference_radius_d_filtered(tag, &1usize, cardinality)
}

/*QUOTIENT HYPERGROUPOIDS */
#[derive(Debug,PartialEq,Clone)]
pub struct QuotientHyperGroupoid{
    pub base_hypergroup:HyperGroupoid,
    pub equivalence_relation:Relation,
    pub hyper_composition:DMatrix<Vec<Vec<u64>>>,
    pub n:u64
}
impl QuotientHyperGroupoid {
    pub fn new_from_equivalence_relation(base_hypergroupoid:&HyperGroupoid,equivalence:&Relation)->Self{

        assert!(equivalence.is_equivalence(),"The input relation is not an equivalence! The quotinet is not defined!");
        let classes = equivalence.quotient_set();
        let n = classes.len() as u64;
        let representants:Vec<u64> = classes.iter().map(|x|x.0).collect();
        let function  = |a:usize,b:usize| 
            support(
                &base_hypergroupoid.mul_by_representation(
                    &(1<<representants[a]), &(1<<representants[b])),&base_hypergroupoid.n).iter()
                    .map(|x|equivalence.get_class(&(*x as u64)).1)
                    .sorted()
                    .unique()
                    .collect_vec();
        let hyper_composition: nalgebra::Matrix<Vec<Vec<u64>>, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<Vec<Vec<u64>>, nalgebra::Dyn, nalgebra::Dyn>> = DMatrix::from_fn(
            n as usize, 
            n as usize, 
            function
            );
        QuotientHyperGroupoid { 
            base_hypergroup:base_hypergroupoid.clone(), 
            equivalence_relation:equivalence.clone(), 
            hyper_composition, 
            n,
            }
       
    }
    
}
impl Display for QuotientHyperGroupoid{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let representants = self.equivalence_relation.quotient_set().iter().map(|x|x.0).collect_vec();
        let table:DMatrix<String>=DMatrix::from_iterator(self.n as usize, self.n as usize, 
            self.hyper_composition.iter().map(|x| format!("{:?}", x.concat())));
        
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} Size:{}\n", representants, table, self.n )
    }
}


