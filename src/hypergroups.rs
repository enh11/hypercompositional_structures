use core::panic;
use std::fmt::Display;
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use crate::{hs::HyperGroupoidMat, utilities::{get_complement_subset, get_subset, ones_positions, representation_permutation_subset, representing_hypergroupoid, vec_to_set, U1024}};
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
        let power_set =2u64.pow(self.0.n as u32);
        (1..power_set-1).into_iter().filter(|x|!x.is_power_of_two()&&self.is_sub_hypergroup(x)).collect_vec()
    }
    pub fn subhypergroup_is_closed(&self,k:&u64)->bool {
        assert!(self.is_sub_hypergroup(k));
        let kc=get_complement_subset(k, &self.0.n);
        let kc_k= self.0.mul_by_representation(&kc, k);
        if kc!=kc_k {return false;}
        let k_kc= self.0.mul_by_representation(k, &kc);
        if kc!=k_kc {return false;}
        true
    }
    /*ADD SUBHYPERGROUP CHECK!!! */
pub fn subhypergroup_is_reflexive(&self,subset_k:u64)->bool {
    self.0.get_singleton().iter().all(|x|self.0.right_division(&subset_k, x)==self.0.left_division(&subset_k, x))
}
pub fn find_reflexive_subhypergroup(&self)->Option<u64>{
    self.collect_proper_subhypergroups()
        .into_iter()
        .find(|x| self.subhypergroup_is_reflexive(*x))
}

}
impl Display for HyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
    
}
