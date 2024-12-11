use std::{collections::HashSet, vec};
use itertools::Itertools;
use rand::Rng;
use crate::utilities::{cartesian_product, get_subset, ones_positions, representing_hypergroupoid, subset_as_u32, to_set};

#[derive(Debug, Clone,PartialEq)]
pub struct HyperStruct{
    pub h:HashSet<u32>,
    pub hyper_composition:Vec<((u32,u32),HashSet<u32>)>
}
impl HyperStruct {
    pub fn new_hypergroup_from_cardinality(n:u32)->Self{
        loop {
            let x=HyperStruct::new_random_from_cardinality(n);
            if x.is_associative()&&x.is_reproductive() {break x;}
        }
    }
    pub fn new_random_from_cardinality(n:u32)->Self{
        let h_vec:Vec<u32>=(0..n).into_iter().collect();
        let h_set:HashSet<u32>=to_set(&h_vec);
        let ht: Vec<((u32, u32), Vec<u32>)> = random_hypercomposition_table(&h_vec);

        let mut ht_set:Vec<((u32, u32), HashSet<u32>)>=Vec::new();
            for item in ht {
                let set=to_set(&item.1);
                ht_set.push((item.0,set));
            }
            HyperStruct { 
                h: h_set, 
                hyper_composition: ht_set
            }
    }
    pub fn elements_multiplication(&self,a:&u32,b:&u32)->HashSet<u32>{
        assert!(self.h.contains(&a),"{} is not an element in H",a);
        assert!(self.h.contains(&b),"{} is not an element in H",b);
        self.hyper_composition[(self.h.len() as u32 *a+b)as usize].1.clone()
    }
    pub fn subsets_multiplication(&self, a:&HashSet<u32>,b:&HashSet<u32>)->HashSet<u32>{
        assert!(a.is_subset(&self.h),"{:?} is not a subset of H",a);
        assert!(b.is_subset(&self.h),"{:?} is not a subset of H",b);
        let mut vec_product:Vec<u32>= Vec::new();
        for x in a {
            for y in b {
                for item in self.elements_multiplication(x,y){
                    vec_product.push(item);
                }
            }
        }
        to_set(&vec_product)
    }
    pub fn singleton(&self)->Vec<HashSet<u32>>{
        //self.h.iter().map(|x| vec![*x]).collect()
        let mut singletons = Vec::new();
        let sort_h:Vec<_>=self.h.clone().into_iter().sorted().collect();
        for item in sort_h {
            singletons.push(HashSet::from([item]));
        }
        singletons
        
    }
    pub fn get_subset_from_k(&self,k:&u32)->HashSet<u32>{
        /*
        k is a number in 0..2^n. We use its binary representation to build a set
        whose elements are the non-zero bits of n*/
        let n = self.h.len() as u32;
        let mut subset: Vec<u32> = Vec::new();
        if k>=&2u32.pow(n){panic!("k can't be grater then 2^n");}
        for i in 0..n {
            if (k >> i)&1==1{
                subset.push(i);
                }
        }
        to_set(&subset)
    }
    
    pub fn is_reproductive(&self)->bool{
        let mut reproductive = false;
        for element in self.singleton() {
            let x_h=self.subsets_multiplication(&element, &self.h);
            let h_x=self.subsets_multiplication(&self.h, &element);
            if(x_h==self.h)&&(h_x==self.h)
                {reproductive=true}
                else {return false;}
        }
        reproductive
    }
    pub fn is_associative(&self)->bool{
        let mut associativity:bool=false;
        for a in self.singleton() {
            for b in self.singleton(){
                for c in self.singleton(){
                    let ab_c=self.subsets_multiplication(&self.subsets_multiplication(&a, &b),&c);
                    let a_bc = self.subsets_multiplication(&a, &&self.subsets_multiplication(&b, &c));
                if a_bc!=ab_c{return false}else{associativity= true}
                }
            }
        }
        associativity
    }
    pub fn subsets_multiplication_represetation(&self, k_subset:&HashSet<u32>,l_subset:&HashSet<u32>)->Vec<(u32,u32)>{
        if !k_subset.is_subset(&self.h)||!k_subset.is_subset(&self.h) { panic!("K must be a subset of H!")};
        let k:u32=subset_as_u32(k_subset);
        let l  =subset_as_u32(l_subset);
        let onse_k=ones_positions(k, &self.h.len());
        let ones_l=ones_positions(l, &self.h.len());
        let mut indexes:Vec<(u32,u32)>=Vec::new();
        for a in onse_k{
            for b in &ones_l{
                    indexes.push((a,*b));
            }
        }
        println!("indexes are {:?}",indexes);
        indexes

    }
    
}
pub fn representation_random_hypercomposition_table(n:&u32)->Vec<Vec<u32>>{
    let vec: Vec<u32>=(0..*n).into_iter().map(|x|x).collect();
    let index_cartesian=cartesian_product(&vec);
    
    println!("indexes are {:?}",index_cartesian);
    let mut rng = rand::thread_rng();
let mut out=vec![vec![0u32;*n as usize];*n as usize];
 for item in index_cartesian {
        out[item.0 as usize][item.1 as usize]=rng.gen_range(0..2u32.pow(*n as u32))
} 
out


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
