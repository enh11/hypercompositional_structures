use core::panic;
use std::fmt::Display;
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use crate::{hs::HyperGroupoidMat, utilities::{get_complement_subset, get_subset, ones_positions, representing_hypergroupoid, vec_to_set, U1024}};
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
        assert!(hg.is_hypergroup(),"Not an hypergroup!");
        HyperGroup(hg)
    }
    pub fn new_from_tag_u1024(tag:&U1024,cardinality:&u64)->Self{
        let hg= HyperGroupoidMat::new_from_tag_u1024(tag,cardinality);
        assert!(hg.is_hypergroup(),"Not an hypergroup!");
        HyperGroup(hg)
    }
    pub fn is_sub_hypergroup(&self,k:&u64)->bool{
        let ones:Vec<usize> = ones_positions(k, &self.0.n).iter().map(|x|*x as usize).collect();
        let sub_matrix = self.0.hyper_composition.select_columns(&ones).select_rows(&ones);
        let sub_hyperstructure=HyperGroupoidMat::new_from_matrix(&sub_matrix);
/*         println!("hyperstructure {}",sub_hyperstructure);
 */        sub_matrix.row_iter().all(
            |row|
                row.iter()
                    .fold(0,|acc,x|acc|x)==*k)
                    &&
        sub_matrix.row_iter().all(
            |row|
                row.iter()
                    .fold(0,|acc,x|acc|x)==*k)

    }
    pub fn collect_proper_subhypergroups(&self)->Vec<u64> {
        let power_set =2u64.pow(self.0.n as u32);
        (1..power_set-1).into_iter().filter(|x|!x.is_power_of_two()&&self.is_sub_hypergroup(x)).collect_vec()
    }
    pub fn is_closed_subhypergroup(&self,k:&u64)->bool {
        assert!(self.is_sub_hypergroup(k));
        let kc=get_complement_subset(k, &self.0.n);
        let kc_k= self.0.mul_by_representation(&kc, k);
        if kc!=kc_k {return false;}
        let k_kc= self.0.mul_by_representation(k, &kc);
        if kc!=k_kc {return false;}
        true



    }

}
impl Display for HyperGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.0)
    }
    
}
