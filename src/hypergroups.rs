use core::panic;
use std::fmt::Display;
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use crate::{hs::{circumference_radius_d_filtered, hg_in_circumference_radius_one, HyperGroupoidMat}, utilities::{get_complement_subset, get_min_max_u1024, get_subset, ones_positions, representation_permutation_subset, representing_hypergroupoid, vec_to_set, U1024}};
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroup(HyperGroupoidMat);

impl HyperGroup {
    pub fn new_from_matrix(matrix:&DMatrix<u64>)->Self{
        let hg= HyperGroupoidMat::new_from_matrix(matrix);
        assert!(hg.is_hypergroup(),"Not an hypergroup!");
        HyperGroup(hg)

    }
    pub fn new_from_tag_u128(tag:&u128,cardinality:&u64)->Self{
        let hg= HyperGroupoidMat::new_from_tag(tag,cardinality);
        assert!(hg.is_hypergroup(),"Not an hypergroup! Tag is {}",tag);
        HyperGroup(hg)
    }
    pub fn new_from_tag_u1024(tag:&U1024,cardinality:&u64)->Self{
        let hg= HyperGroupoidMat::new_from_tag_u1024(tag,cardinality);
        assert!(hg.is_hypergroup(),"Not an hypergroup!");
        HyperGroup(hg)
    }
    pub fn cardinality(&self)->u64{
        self.0.n
    }
    pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
        let alpha = self.0.permutation_of_table(sigma).get_integer_tag_u1024();
        HyperGroup::new_from_tag_u1024(&alpha, &self.cardinality())

    }
    pub fn collect_isomorphism_class(&self)->(U1024,Vec<U1024>){
        self.0.collect_isomorphism_class()
        
    }
    pub fn isomorphic_hypergroup_from_permutation(&self, sigma:&Permutation)->Self{
        HyperGroup::new_from_tag_u1024(&self.0.isomorphic_hypergroup_from_permutation(sigma).get_integer_tag_u1024(),&self.0.n)
    }
    pub fn get_integer_tag_u1024(&self)->U1024{
        self.0.get_integer_tag_u1024()
    }
    pub fn is_transposition(&self)->bool {
        for a in self.0.get_singleton(){
            for b in self.0.get_singleton(){
                for c in self.0.get_singleton(){
                    for d in self.0.get_singleton(){
                        if self.0.left_division(&a, &b)&self.0.right_division(&c, &d)==1{
                            if self.0.mul_by_representation(&a, &d)&self.0.mul_by_representation(&b, &c)!=1{
                                return  false;
                            }
                        
                        }
                    }
                }
            }
        }
        true
    }
    pub fn is_sub_hypergroup(&self,k:&u64)->bool{
        let ones:Vec<usize> = ones_positions(k, &self.0.n).iter().map(|x|*x as usize).collect();
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
    pub fn collect_proper_subhypergroups(&self)->Vec<u64> {
        let power_set =2u64.pow(self.cardinality() as u32);
        (1..power_set-1).into_iter().filter(|x|!x.is_power_of_two()&&self.is_sub_hypergroup(x)).collect_vec()
    }
    pub fn subhypergroup_is_closed(&self,subset_k:&u64)->bool {
        if !self.is_sub_hypergroup(&subset_k) {return false;}
        let kc=get_complement_subset(subset_k, &self.0.n);
        let kc_k= self.0.mul_by_representation(&kc, subset_k);
        if kc!=kc_k {return false;}
        let k_kc= self.0.mul_by_representation(subset_k, &kc);
        if kc!=k_kc {return false;}
        true
    }
pub fn subhypergroup_is_reflexive(&self,subset_k:u64)->bool {
    if !self.is_sub_hypergroup(&subset_k) {return false;}
    self.0.get_singleton().iter().all(|x|self.0.right_division(&subset_k, x)==self.0.left_division(&subset_k, x))
}
pub fn collect_reflexive_subhypergroup(&self)->Vec<u64>{
    self.collect_proper_subhypergroups()
        .into_iter()
        .filter(|x| self.subhypergroup_is_reflexive(*x))
        .collect()
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
//println!("cl1 {:?}",cl1);
    classes.append(&mut cl1.1.clone());
//println!("classes {:?}",classes);
    center_tags=cl1.1.iter().map(|x|x.0).collect();
    visited_tags.append(&mut cl1.0);
// println!("new centers  {:?}",center_tags);

    while !center_tags.is_empty() {
        centers_and_classes=center_tags.par_iter().map(|x|
            collect_classes_from_circumference(x, &visited_tags, &cardinality)
            ).collect();
        classes_from_circumferences = centers_and_classes.iter().map(|y|y.1.clone()).collect();
        classes_from_circumferences.concat();
        classes_from_circumferences.retain(|x|!x.is_empty());
        classes_from_circumferences.sort();
        classes_from_circumferences.dedup();
        //println!("classes_from tags in center {:?}",classes_from_circumferences);

        classes.append(&mut classes_from_circumferences.clone().concat());
        classes.sort_by(|x,y|x.0.cmp(&y.0));
        classes.dedup();
        //println!("update classes {:?}",classes);
        let visited_tags_vec=centers_and_classes.iter().map(|y|y.0.clone()).collect_vec();
        center_tags=visited_tags_vec.concat();
        center_tags.sort();
        center_tags.dedup();
        //println!("new centers {:?}",center_tags);
        visited_tags.append(&mut center_tags.clone());
        visited_tags.sort();
        visited_tags.dedup();
        //println!("visited tags are {:?}",visited_tags);
        
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
