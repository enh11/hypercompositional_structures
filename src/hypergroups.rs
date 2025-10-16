use core::panic;
use std::{collections::{HashMap, HashSet}, fmt::{self, Display}};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use num_rational::Rational64;
use permutation::Permutation;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use crate::{binary_relations::relations::Relation, fuzzy::FuzzySubset, hs::{circumference_radius_d_filtered, hg_in_circumference_radius_one, HyperGroupoid}, quotient_hg::QuotientHyperGroup, utilities::{get_complement_subset, support, u64_to_set, vec_to_set, U1024}};
#[derive(Debug, Clone)]
pub enum HyperStructureError {
    NotHypergroup,
}
impl fmt::Display for HyperStructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HyperStructureError::NotHypergroup => write!(f, "The structure is not a valid hypergroup."),
        }
    }
}
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroup(pub HyperGroupoid);

impl HyperGroup {
    pub fn new_from_hypergroupoid(h : HyperGroupoid)->Self{
        match h.try_into() {
            Ok(hg) => HyperGroup(hg),
            _ => panic!("The input hypergroupoid is not a hypergroup."),
        }
    }
    pub fn new_from_matrix(matrix:&DMatrix<u64>)->Self{
        let hg= HyperGroupoid::new_from_matrix(matrix);
        HyperGroup::new_from_hypergroupoid(hg)

    }
    pub fn new_from_tag_u128(tag:&u128,cardinality:&u64)->Self{
        let hg= HyperGroupoid::new_from_tag_u128(tag,cardinality);
        HyperGroup::new_from_hypergroupoid(hg)
    }
    pub fn new_from_tag_u1024(tag:&U1024,cardinality:&u64)->Self{
        let hg= HyperGroupoid::new_from_tag_u1024(tag,cardinality);
        HyperGroup::new_from_hypergroupoid(hg)
    }
    pub fn new_from_elements(input_array: &Vec<Vec<u64>>, cardinality:&u64)->Self{
        let hg = HyperGroupoid::new_from_elements(input_array, &cardinality);
        HyperGroup::new_from_hypergroupoid(hg)

    }
/// Generates a new hypergroup from a binary function Fn(u64, u64) -> u64.
///
/// The function is sent to HyperGroupoids::new_from_function`, which compute the new hyperstructure and check if it is a hypergroup.
/// If it is, it returns OK(HyperGroup), otherwise an Error;
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
/// use hyperstruc::hypergroups::HyperGroup;
/// 
/// let cardinality =3u64;
/// let function = |a:u64,b:u64| 1<<a|1<<b;
/// let hg = HyperGroup::new_from_function(function, &cardinality);
/// assert!(hg.unwrap().is_transposition());
///
pub fn new_from_function<F>(function:F,cardinality:&u64)->Result<Self,HyperStructureError>
    where F: Fn(u64,u64) -> u64{
        let hs = HyperGroupoid::new_from_function(function, cardinality);
        match hs.is_hypergroup() {
            true => Ok(HyperGroup::new_from_tag_u1024(&hs.get_integer_tag_u1024(), cardinality)),
            false => Err(HyperStructureError::NotHypergroup),            
        }
}
pub fn cardinality(&self)->u64{
        self.0.n
}
pub fn get_singleton(&self)->Vec<u64>{
    self.0.get_singleton()
}
pub fn mul_by_representation(&self, subset_a:&u64,subset_b:&u64)->u64{
    self.0.mul_by_representation(subset_a, subset_b)
}
pub fn left_division(&self, subset_a:&u64,subset_b:&u64)->u64{
    self.0.left_division(subset_a, subset_b)
}
pub fn right_division(&self, subset_a:&u64,subset_b:&u64)->u64{
    self.0.right_division(subset_a, subset_b)
}
pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
    let alpha = self.0.permutation_of_table(sigma).get_integer_tag_u1024();
    HyperGroup::new_from_tag_u1024(&alpha, &self.cardinality())
}
pub fn collect_isomorphism_class(&self)->(U1024,Vec<U1024>){
    self.0.collect_isomorphism_class()
}
pub fn isomorphic_hypergroup_from_permutation(&self, sigma:&Permutation)->Self{
    HyperGroup::new_from_tag_u1024(&self.0.isomorphic_hypergroup_from_permutation(sigma).get_integer_tag_u1024(),&self.cardinality())
}
///
/// Return true if the two hypergroups are isomorphic. 
/// 
/// # Example
/// ```
/// use hyperstruc::utilities::U1024;
/// use hyperstruc::hypergroups::HyperGroup;
/// 
/// let cardinality = 3u64;
/// let tag_1 = U1024::from(22097724u128);
/// 
/// let hg_1=HyperGroup::new_from_tag_u1024(&tag_1,&cardinality);
/// let tag_2 = U1024::from(31958100u128);
/// let hg_2=HyperGroup::new_from_tag_u1024(&tag_2,&cardinality);
/// 
/// assert!(hg_1.is_isomorphic_to(&hg_2))
/// 
pub fn is_isomorphic_to(&self,other:&Self)->bool{
    self.0.is_isomorphic_to(&other.0)
}
pub fn get_integer_tag_u1024(&self)->U1024{
    self.0.get_integer_tag_u1024()
}
///
/// 
/// Return true if the hypergroup is transposition, i.e., if `a\b meets c/d implies ad meets bc`.
///
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hypergroups::HyperStructureError;
/// use hyperstruc::hypergroups::HyperGroup;
/// 
/// let cardinality =5u64;
/// let function = |a:u64,b:u64| 1<<a|1<<b;
/// let hg = match HyperGroup::new_from_function(function, &cardinality) {
/// Ok(hg)=>hg,
/// Err(HyperStructureError::NotHypergroup)=> panic!("Not hg")
/// };
/// assert!(hg.is_transposition())
/// 
/// 
pub fn is_transposition(&self)->bool {
    for a in self.get_singleton(){
        for b in self.get_singleton(){
            for c in self.get_singleton(){
                for d in self.get_singleton(){
                    if self.left_division(&a, &b)&self.right_division(&c, &d)!=0{
                        if self.mul_by_representation(&a, &d)&self.mul_by_representation(&b, &c)==0{
                            return  false;
                        }
                    }
                }
            }
        }
    }
    true
}
pub fn is_commutative(&self)->bool{
    self.0.is_commutative()
}
pub fn is_quasicanonical(&self)->bool{
    if !self.is_transposition(){return false;}
    !self.0.collect_scalar_identities().is_empty() //if a scalar identity exists, it is unique.
}
pub fn is_join_space(&self)->bool{
    self.is_commutative()&&self.is_transposition()
}
pub fn is_canonical(&self)->bool{
    self.is_commutative()&&self.is_quasicanonical()
}
pub fn is_sub_hypergroup(&self,k:&u64)->bool{
    let power_set_cardinality = 1<<self.cardinality();
    assert!(*k<power_set_cardinality,"K is not a subset of H!");
    let ones:Vec<usize> = support(k, &self.0.n).iter().map(|x|*x as usize).collect();
    let sub_matrix = self.0.hyper_composition.select_columns(&ones).select_rows(&ones);
        sub_matrix.row_iter().all(
            |row|
                row
                    .iter()
                    .fold(0,|acc,x|acc|x)==*k)
        &&
        sub_matrix.row_iter().all(
            |row|
                row
                    .iter()
                    .fold(0,|acc,x|acc|x)==*k)
    }
/// 
/// 
/// Returns the vector of all identities of the hypergroup. They are elements in H and are represented as integers in [0,2^n-1].
/// Therefore, they are powers of two. To identify them as elements use u.trailing_zeros() ore 
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use nalgebra::DMatrix;
/// let matrix=DMatrix::from_row_slice(3usize,3usize,&[1,2,4,1,2,4,7,7,7]);
/// let hyperstructure=HyperGroupoid::new_from_matrix(&matrix);
/// let identities = hyperstructure.collect_identities();
/// assert!(identities.is_none())
/// 
///
pub fn collect_identities(&self)->Option<Vec<u64>>{
    self.0.collect_identities()
}
pub fn collect_left_identities(&self)-> Option<Vec<u64>> {
    self.0.collect_left_identities()
}
pub fn collect_right_identities(&self)-> Option<Vec<u64>> {
    self.0.collect_right_identities()
}
pub fn collect_scalars(&self)->Vec<u64>{
    self.0.collect_scalars()
}
pub fn collect_scalars_identities(&self)->Vec<u64>{
    self.0.collect_scalar_identities()
}
pub fn collect_left_scalars(&self)->Vec<u64> {
    self.0.collect_left_scalars()
}
pub fn collect_right_scalars(&self)-> Vec<u64> {
    self.0.collect_right_scalars()
}
pub fn collect_partial_identities(&self)-> Vec<u64> {
    self.0.collect_partial_identities()
}
pub fn left_inverses_of_x(&self, x:&u64,u:&u64)->u64{
    self.0.left_inverses_of_x(x, u)
}
pub fn right_inverses_of_x(&self,x:&u64,u:&u64)->u64{
    self.0.right_inverses_of_x(x, u)
}
pub fn collect_right_inverses_of_x(&self,x:&u64)->Option<Vec<(u64,u64)>>{
    self.0.collect_right_inverses_of_x(x)
}
pub fn collect_left_inverses_of_x(&self,x:&u64)->Option<Vec<(u64,u64)>>{
    self.0.collect_left_inverses_of_x(x)
}
pub fn collect_inverses_of_x(&self,x:&u64)->Option<Vec<(u64,u64)>>{
    self.0.collect_inverses_of_x(x)
}
pub fn get_corsini_fuzzysubset(&self)->FuzzySubset{
        self.0.get_corsini_fuzzysubset()
    }
pub fn collect_proper_subhypergroups(&self)->Option<Vec<u64>> {
    let power_set_cardinality: u64 =(1<<self.cardinality())-1;
    let sub_hg:Vec<u64>=(1..power_set_cardinality)
        .into_iter()
        .filter(|x|self.is_sub_hypergroup(x))
        .collect();
    match sub_hg.is_empty() {
        true => None,
        false => Some(sub_hg)
    }
}
pub fn show_proper_subhypergroups(&self) {
match  self.collect_proper_subhypergroups(){
    Some(sub_hgs) => {
        println!("Proper subhypergroups of H are:\n{}",
        (0..sub_hgs.len()).into_iter().map(|i|
            format!("K_{} = {:?}",
            i,
            u64_to_set(&sub_hgs[i],&self.cardinality()))
        ).join("\n")
    )
    },
    None => println!("H has no non-trivial subhypergroups."),
}
}
pub fn subhypergroup_is_closed(&self,subset_k:&u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    let kc=get_complement_subset(subset_k, &self.cardinality());
    let kc_k= self.mul_by_representation(&kc, subset_k);
        if kc!=kc_k {return false;}
    let k_kc= self.mul_by_representation(subset_k, &kc);
        if kc!=k_kc {return false;}
    true
    }
pub fn subhypergroup_is_reflexive(&self,subset_k:&u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    self.get_singleton()
        .iter()
        .all(|x|
            self.right_division(&subset_k, x)
            ==
            self.left_division(&subset_k, x))
}

///
/// 
/// Return true if the sub-hypergroup is `normal`, i.e., if `aN = Na` hold for any `a` in `H`.
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hypergroups::HyperStructureError;
/// use hyperstruc::hypergroups::HyperGroup;
/// 
/// let cardinality =5u64;
/// let function = |a:u64,b:u64| 1<<a|1<<b;
/// let hg = match HyperGroup::new_from_function(function, &cardinality) {
/// Ok(hg)=>hg,
/// Err(HyperStructureError::NotHypergroup)=> panic!("Not hg")
/// };
/// let subhypergroup = 10u64; //This is the subset K = {1,3} of H = {0,1,2,3,4}
/// assert!(hg.is_sub_hypergroup(&subhypergroup));
/// assert!(hg.subhypergroup_is_normal(&subhypergroup));
/// 
/// 
pub fn subhypergroup_is_normal(&self,subset_k:&u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    self.get_singleton()
        .iter()
        .all(|x|
            self.mul_by_representation(&subset_k, x)
            ==
            self.mul_by_representation(&subset_k, x))
}
pub fn subhypergroup_is_right_invertible(&self,subset_k:&u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    self.get_singleton()
        .iter()
        .cartesian_product(self.get_singleton())
        .into_iter()
        .all(|(x,y)|
            ((x&self.mul_by_representation(subset_k, &y))!=*x) 
            || 
            ((y&self.mul_by_representation(subset_k, x))==y)
            )
}
pub fn subhypergroup_is_left_invertible(&self,subset_k:&u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    self.get_singleton()
        .iter()
        .cartesian_product(self.get_singleton())
        .into_iter()
        .all(|(x,y)|
            ((x&self.mul_by_representation(&y, subset_k))!=*x) 
            || 
            ((y&self.mul_by_representation( x,subset_k))==y)
            )
}
pub fn subhypergroup_is_invertible(&self,subset_k:&u64)->bool {
    self.subhypergroup_is_left_invertible(subset_k)&&self.subhypergroup_is_right_invertible(subset_k)
}
pub fn collect_proper_invertible_subhypergroups(&self)->Option<Vec<u64>>{
    match self.collect_proper_subhypergroups() {
        Some(sub_hgs) => {
            let sub_invertible:Vec<u64> = sub_hgs.into_iter()
                .filter(|x|
                self.subhypergroup_is_invertible(&x)
                )
            .collect();
            match sub_invertible.is_empty() {
                true => None,
                false => Some(sub_invertible)
            }
        },
        None => None,
    }
}
pub fn show_proper_invertible_subhypergroups(&self) {
match  self.collect_proper_invertible_subhypergroups(){
    Some(sub_hg_invertible) => {
        println!("Proper invertible subhypergroups of H are:\n{}",
        (0..sub_hg_invertible.len()).into_iter().map(|i|
            format!("I_{} = {:?}",
            i,
            u64_to_set(&sub_hg_invertible[i],&self.cardinality()))
        ).join("\n")
    )
    },
    None => println!("H has no non-trivial invertible subhypergroup."),
}
}
pub fn collect_proper_closed_subhypergroups(&self)->Option<Vec<u64>>{
    match self.collect_proper_subhypergroups() {
        Some(sub_hgs) => {
            let sub_closed:Vec<u64> = sub_hgs.into_iter()
                .filter(|x|
                self.subhypergroup_is_closed(&x)
                )
            .collect();
            match sub_closed.is_empty() {
                true => None,
                false => Some(sub_closed)
            }
        },
        None => None,
    }
}
pub fn show_proper_closed_subhypergroups(&self) {
match  self.collect_proper_closed_subhypergroups(){
    Some(sub_hg_closed) => {
        println!("Proper closed subhypergroups of H are:\n{}",
        (0..sub_hg_closed.len()).into_iter().map(|i|
            format!("C_{} = {:?}",
            i,
            u64_to_set(&sub_hg_closed[i],&self.cardinality()))
        ).join("\n")
    )
    },
    None => println!("H has no non-trivial closed subhypergroup."),
}
}
pub fn collect_proper_reflexive_subhypergroups(&self)->Option<Vec<u64>>{
    match self.collect_proper_subhypergroups() {
        Some(sub_hgs) => {
            let sub_reflexive:Vec<u64> = sub_hgs.into_iter()
                .filter(|x|
                self.subhypergroup_is_reflexive(&x)
                )
            .collect();
            match sub_reflexive.is_empty() {
                true => None,
                false => Some(sub_reflexive)
            }
        },
        None => None,
    }
}
pub fn show_proper_reflexive_subhypergroups(&self) {
match  self.collect_proper_reflexive_subhypergroups(){
    Some(sub_hg_reflexive) => {
        println!("Proper reflexive subhypergroups of H are:\n{}",
        (0..sub_hg_reflexive.len()).into_iter().map(|i|
            format!("R_{} = {:?}",
            i,
            u64_to_set(&sub_hg_reflexive[i],&self.cardinality()))
        ).join("\n")
    )
    },
    None => println!("H has no non-trivial reflexive subhypergroup."),
}
}
pub fn collect_proper_normal_subhypergroups(&self)->Option<Vec<u64>>{
    match self.collect_proper_subhypergroups() {
        Some(sub_hgs) => {
            let sub_normal:Vec<u64> = sub_hgs.into_iter()
                .filter(|x|
                self.subhypergroup_is_normal(&x)
                )
            .collect();
            match sub_normal.is_empty() {
                true => None,
                false => Some(sub_normal)
            }
        },
        None => None,
    }
}
pub fn show_proper_normal_subhypergroups(&self) {
match  self.collect_proper_normal_subhypergroups(){
    Some(sub_hg_normal) => {
        println!("Proper normal subhypergroups of H are:\n{}",
        (0..sub_hg_normal.len()).into_iter().map(|i|
            format!("N_{} = {:?}",
            i,
            u64_to_set(&sub_hg_normal[i],&self.cardinality()))
        ).join("\n")
    )
    },
    None => println!("H has no non-trivial normal subhypergroup."),
}
}

///
/// Explore the graph of hypergroups starting from it. It returns a tuple ((representant,classes),enumeration of hypergroup in the tree)).
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::hypergroups::HyperGroup;
/// use hyperstruc::utilities::get_min_max_u1024;
/// use hyperstruc::hypergroups::collect_classes_with_respect_to_cardinality;
/// 
/// let cardinality = 2u64;
/// let total_tag = get_min_max_u1024(&cardinality).1;
/// let total_hg = HyperGroup::new_from_tag_u1024(&total_tag,&cardinality);
/// let (classes, enumeration)=total_hg.exploring_tree();
/// assert_eq!(enumeration.len(),14);
/// 
/// let enumeration_classes = collect_classes_with_respect_to_cardinality(&classes,&cardinality).1;
/// let expected_enumeration = [2,6].to_vec();
/// assert_eq!(enumeration_classes,expected_enumeration);
/// 
/// 
pub fn exploring_tree(&self)->(Vec<(U1024,Vec<U1024>)>,Vec<U1024>) {
    exploring_tree(&self.get_integer_tag_u1024(), &self.cardinality())
    
}
pub fn get_qq_u(&self,u:&u64)->Vec<(u64,u64)>{
    self.0.get_qq_u(u)
}
pub fn get_alpha_u(&self,u:&u64)-> Rational64{
    self.0.get_alpha_u(u)
}
pub fn get_mu_u(&self,u:&u64)->Rational64{
    self.0.get_mu_u(u)
}
pub fn get_fuzzy_grade(&self)->usize{
    self.0.get_fuzzy_grade()
}
pub fn get_strong_fuzzy_grade(&self)->usize{
    self.0.get_strong_fuzzy_grade()
}
pub fn beta_relation(&self)->Relation {
    self.0.beta_relation()
}
pub fn collect_beta_classes(&self)->Vec<(u64,Vec<u64>)>{
    self.0.beta_relation().quotient_set()
}
pub fn get_fundamental_group(&self)-> QuotientHyperGroup {
    let beta = self.beta_relation();
    QuotientHyperGroup::new_from_equivalence_relation(self.clone(), beta)
}
/// Returns the isomorphic image of the fundamental group as a standard hypergroup, 
/// represented over the set `H = {0, 1, ..., n-1}`, where `n` is the cardinality of 
/// the quotient set.
///
/// This function constructs a new hypergroup by ordering the equivalence classes of the 
/// fundamental group and assigning each class a unique integer identifier from `0` to `n - 1`. 
/// The result is a standard representation of the quotient structure, isomorphic to the 
/// original fundamental group.
///
/// # Example
/// ```
/// use hyperstruc::hypergroups::HyperGroup;
/// use nalgebra::DMatrix;
///
/// let cardinality = 4u64;
/// let hs = HyperGroup::new_from_matrix(
///     &DMatrix::from_row_slice(
///         cardinality as usize,
///         cardinality as usize,
///         &[1, 6, 6, 8,
///           6, 8, 8, 1,
///           6, 8, 8, 1,
///           8, 1, 1, 6]
///     )
/// );
///
/// let fundamental_group = hs.get_fundamental_group();
/// println!("fundamental {}", fundamental_group);
///
/// let fundamental_group = hs.get_isomorphic_fundamental_group();
/// let cardinality_fg = 3u64;
///
/// let expected_fundamental_group = HyperGroup::new_from_elements(
///     &vec![
///         vec![0], vec![1], vec![2],
///         vec![1], vec![2], vec![0],
///         vec![2], vec![0], vec![1],
///     ],
///     &cardinality_fg,
/// );
///
/// assert_eq!(fundamental_group, expected_fundamental_group);
/// ```
///
/// # Panics
/// This function may panic if the internal composition structure is invalid, 
/// such as missing equivalence class mappings or inconsistent dimensions.
///
/// # See Also
/// - [`get_fundamental_group`](Self::get_fundamental_group): Returns the untransformed fundamental group.
/// - [`HyperGroup::new_from_elements`](HyperGroup::new_from_elements): Constructs a hypergroup from raw data.
pub fn get_isomorphic_fundamental_group(&self)->HyperGroup{
    let fg = self.get_fundamental_group();
    let mut classes:Vec<Vec<u64>>= self.beta_relation().quotient_set().into_iter().map(|(_,y)|y).collect();
    classes.sort();
    let cardinality = classes.len();

    let hash:HashMap<Vec<u64>,usize>= classes.into_iter().enumerate().map(|(i, x)| (x, i)).collect();    
    let generating_function = |a:u64,b:u64| 
           {
           let get =  match hash.get(
                &fg.0.hyper_composition[(a as usize,b as usize)][0]) {
                    Some(x)=>*x,
                    None=> panic!("not found!")               
           };

           1<<get as u64
           };

    HyperGroup::new_from_function(generating_function, &(cardinality as u64)).unwrap()
}
/// Compute the Heart of a hypergroup. This is the β*‐class-class of the identity of the fundamental group  H / β*.
/// This method returns the heart as a set of elements, rather than its integer representation.
///
/// # Returns
/// A `Option<HashSet<usize>>` representing the heart of the hypergroup.
///
/// # Example
/// ```
/// use hyperstruc::hypergroups::HyperGroup;
/// use std::collections::HashSet;
///
/// let cardinality = 3u64;
/// let tag = 57784611u128;
/// let hg = HyperGroup::new_from_tag_u128(&tag, &cardinality);
///
/// let heart = hg.heart();
/// let expected = Some(HashSet::from([0, 1]));
///
/// assert_eq!(heart, expected);
/// ```
pub fn heart(&self)->Option<HashSet<u64>>{
    let fundamental_group = self.get_isomorphic_fundamental_group();
    let identity = fundamental_group.collect_identities();
    //assert!(identity.len()==1);
    match identity {
        Some(e) => {
            let id = e[0] as u64;
            let beta_identity = self.beta_relation().get_class(&id);
                Some(vec_to_set(&beta_identity.1))}

        None => None,
    }
}
pub fn heart_fast(&self)->HashSet<u64> {
    let pid = self.collect_partial_identities();
    let pid = pid[0].trailing_zeros() as u64; //convert from singleton to element in H
     self.beta_relation().get_class(&pid).1.into_iter().collect()
}
}
impl Display for HyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
    
}
pub fn collect_classes_from_circumference(center_tag:&U1024,visited_tags:&Vec<U1024>,cardinality:&u64)->(Vec<U1024>,Vec<(U1024,Vec<U1024>)>){
    let mut circumferences: Vec<Vec<U1024>>=Vec::new();
    let mut circumference_radius_1=hg_in_circumference_radius_one(center_tag, &cardinality);
    circumference_radius_1.retain(|x|!visited_tags.contains(x));
    circumferences.push(circumference_radius_1.clone());
   
    let classes_from_circumferences:Vec<_>=circumferences.par_iter().map(
        |x|x.iter().map(|tag|
                HyperGroup::new_from_tag_u1024(tag, &cardinality)
                .collect_isomorphism_class())
            .collect_vec())
            .collect();
    let mut classes_from_circumferences=classes_from_circumferences.concat();
    classes_from_circumferences.retain(|x|!x.1.is_empty());
    classes_from_circumferences.sort_by(|x,y|x.0.cmp(&y.0));
    classes_from_circumferences.dedup();
    (circumference_radius_1,classes_from_circumferences)
}
pub fn collect_isomorphism_class(tag:&U1024,cardinality:&u64)->(U1024,Vec<U1024>){
    let hg  = HyperGroup::new_from_tag_u1024(tag, cardinality);
    hg.0.collect_isomorphism_class()
}
pub fn tag_has_sons(tag:&U1024,cardinality:&u64)->bool{
    hg_in_circumference_radius_one(tag, cardinality).iter().any(|x|x<tag)
}
pub fn tag_is_leaf(tag:&U1024,cardinality:&u64)->bool {
    !tag_has_sons(tag, cardinality)
}
pub fn collect_leafs_from_nodes(nodes: &Vec<U1024>,cardinality:&u64)->Vec<U1024>{
    nodes.par_iter()
        .filter(|x|
            tag_is_leaf(*x, cardinality)
                ).map(|x|*x)
                .collect()
}
///
/// Explore the graph of hypergroups starting from its tag. It returns a tuple ((representant,classes),enumeration of hypergroup in the tree)).
/// 
/// # Example
/// 
/// ```
/// use hyperstruc::utilities::get_min_max_u1024;
/// use hyperstruc::hypergroups::exploring_tree;
/// 
/// let cardinality = 2u64;
/// let total_tag = get_min_max_u1024(&cardinality).1;
/// let (classes, enumeration)=exploring_tree(&total_tag,&cardinality);
/// assert_eq!(enumeration.len(),14);
/// 
/// 
pub fn exploring_tree(starting_tag:&U1024,cardinality:&u64)->(Vec<(U1024,Vec<U1024>)>,Vec<U1024>){
    let mut visited_tags: Vec<U1024> = Vec::new();
    let mut center_tags:Vec<U1024> = Vec::new();
    center_tags.push(*starting_tag);
    visited_tags.push(*starting_tag);
    let mut centers_and_classes:Vec<(Vec<U1024>,Vec<(U1024,Vec<U1024>)>)>;
    let mut classes_from_circumferences :Vec<Vec<(U1024,Vec<U1024>)>>;
    let mut classes:Vec<(U1024,Vec<U1024>)> = [collect_isomorphism_class(starting_tag, cardinality)].to_vec();
    let mut cl1 = collect_classes_from_circumference(starting_tag, &visited_tags, &cardinality);
    classes.append(&mut cl1.1.clone());
    center_tags=cl1.1.iter().map(|x|x.0).collect();
    visited_tags.append(&mut cl1.0);

    while !center_tags.is_empty() {
        centers_and_classes=center_tags.par_iter().map(|x|
            collect_classes_from_circumference(x, &visited_tags, &cardinality)
            ).collect();
        classes_from_circumferences = centers_and_classes.iter().map(|y|y.1.clone()).collect();
        classes_from_circumferences.concat();
        classes_from_circumferences.retain(|x|!x.is_empty());
        classes_from_circumferences.sort();
        classes_from_circumferences.dedup();

        classes.append(&mut classes_from_circumferences.clone().concat());
        classes.sort_by(|x,y|x.0.cmp(&y.0));
        classes.dedup();
        let visited_tags_vec=centers_and_classes.iter().map(|y|y.0.clone()).collect_vec();
        center_tags=visited_tags_vec.concat();
        center_tags.sort();
        center_tags.dedup();
        visited_tags.append(&mut center_tags.clone());
        visited_tags.sort();
        visited_tags.dedup();
        
    }
    let enumeration = classes.iter().map(|x|x.1.clone()).collect_vec();
    let enumeration = enumeration.concat();
    
    (classes,enumeration)
}
pub fn find_nodes_not_enumerated(leafs:&Vec<U1024>,enumeration:&Vec<U1024>,radius:&usize,cardinality:&u64)->Vec<U1024>{
    leafs.par_iter().filter(|tag|
        circumference_radius_d_filtered(&U1024::from(**tag), radius, &cardinality).iter().any(|x|
        !enumeration.contains(&U1024::from(x)))
    ).map(|x|*x).collect()
    
}
pub fn collect_classes_with_respect_to_cardinality(classes:&Vec<(U1024,Vec<U1024>)>,cardinality:&u64)->(Vec<Vec<(U1024,Vec<U1024>)>>,Vec<usize>){
    let permut_vec:Vec<Vec<usize>> = (0..*cardinality as usize).permutations(*cardinality as usize).collect();
    let mut enumeration_classes:Vec<usize>=Vec::new();
    let mut c_k:Vec<(U1024,Vec<U1024>)>;
    let mut collected_classes:Vec<Vec<(U1024,Vec<U1024>)>>=Vec::new();

    for k in 1..=permut_vec.len(){
        c_k=classes.iter().filter(|y|(*y.1).len()==k).into_iter().map(|x|x.clone()).collect_vec();
        enumeration_classes.push(c_k.len());
        collected_classes.push(c_k);
       
    }
    (collected_classes,enumeration_classes)
}

