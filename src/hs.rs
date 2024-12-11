use core::panic;
use std::{collections::HashSet, fmt::Display, vec};
extern crate nalgebra as na;
use itertools::Itertools;
use nalgebra::DMatrix;
use permutation::Permutation;
use rand:: Rng;
use crate::utilities::{cartesian_product, get_subset, ones_positions, representation_permutation_subset, representing_hypergroupoid, subset_as_u32, to_set};
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroupoid{
    pub h:HashSet<u32>,
    pub hyper_composition:Vec<Vec<u32>>,
    pub n:u32
}
#[derive(Debug, Clone,PartialEq)]
pub struct HyperGroupoidMat{
    pub h:HashSet<u32>,
    pub hyper_composition:DMatrix<u32>,
    pub n:u32
}
impl HyperGroupoidMat {
   pub fn new_random_from_cardinality(n:&u32)->Self{
    let h_vec=(0..*n as u32).into_iter().collect();
        let ht = get_random_hypercomposition_matrix(&n);
        HyperGroupoidMat { 
            h: h_vec, 
            hyper_composition: ht, 
            n: *n}
   }
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
pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
    let permutation_hypergroupoids = &self.hyper_composition;
    let alfa =DMatrix::from_iterator(self.n as usize, self.n as usize, 
        permutation_hypergroupoids.iter()
            .map(|x| representation_permutation_subset(x,&sigma)));
    
    HyperGroupoidMat { 
        h: self.h.clone(), 
        hyper_composition:alfa, 
        n: self.n}
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
pub fn get_subset_from_k(&self,k:&u32)->HashSet<u32>{
    /*
    k is a number in 0..2^n-1. We use its binary representation to build a set
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
    /*This function compute the value b\a={x in H s.t. a in bx} */
    
    let sub_a=to_set(&get_subset(&2u32.pow(*a), &self.n));
    let sub_b=2u32.pow(*b);
    self.get_singleton().iter()
    .filter(
        |x| sub_a.is_subset(
            &to_set(&get_subset(
                        &self.mul_by_representation(&sub_b, x), &self.n)
                    )
                )
            ).fold(0, |acc,t|acc|t)

   
}
pub fn right_division(&self,a:&u32,b:&u32)->u32{
        /*This function compute the value a/b={x in H s.t. a in xb} */
        let sub_a=to_set(&get_subset(&2u32.pow(*a), &self.n));
    let sub_b=2u32.pow(*b);
    self.get_singleton().iter()
    .filter(
        |x| sub_a.is_subset(
            &to_set(&get_subset(
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
        println!("xH = {:?}, Hx = {:?}",row_sum,col_sum);
        false
    }
   }
pub fn fix_reproductivity(&self)->Self{
    let h=2u32.pow(self.n)-1;
    if self.is_reproductive(){return self.clone();}
    let mut new_hypergroupoid_matrix = self.hyper_composition.clone();
    let mut missing_terms:u32;
    let row_sum:Vec<u32> = new_hypergroupoid_matrix.row_iter().map(|x|x.iter().fold(0u32, |acc,element|acc|element)).collect();
    
        for j in 0..self.n as usize /*here j is row index*/{
            if row_sum[j]<h {
                missing_terms=h-row_sum[j];
                let position_max = new_hypergroupoid_matrix.row(j).iter().position_max().unwrap();
                new_hypergroupoid_matrix[(j,position_max as usize)]|=missing_terms;
            }
            
        }            

    let new_hypergroupoid = HyperGroupoidMat{
        h:self.h.clone(),
        hyper_composition:new_hypergroupoid_matrix.clone(),
        n:self.n
    };
    if new_hypergroupoid.is_reproductive(){return new_hypergroupoid;}
    else {
        let col_sum:Vec<u32> = new_hypergroupoid_matrix.column_iter().map(|x|x.iter().fold(0u32, |acc,element|acc|element)).collect();
    
        for j in 0..self.n as usize{
            if col_sum[j]<h {
                missing_terms=h-col_sum[j];
                let position_max = new_hypergroupoid_matrix.column(j).iter().position_max().unwrap();
                new_hypergroupoid_matrix[(position_max as usize,j as usize)]|=missing_terms;
            }         
    }

    }
    HyperGroupoidMat{
        h:self.h.clone(),
        hyper_composition:new_hypergroupoid_matrix,
        n:self.n
    }
}
pub fn is_associative(&self)->bool{
    let mut associativity = true;
    for a in &self.get_singleton(){
        for b in &self.get_singleton(){
            for c in &self.get_singleton(){
                let ab_c=self.mul_by_representation(
                    &self.mul_by_representation(&a, &b),&c);
                let a_bc = self.mul_by_representation(&a, &self.mul_by_representation(&b, &c));
                if a_bc==ab_c{continue;}else {
                        println!("Not associative hypergroupoid:\n{:?}{:?}_{:?}= {:?}\n {:?}_{:?}{:?} = {:?}", to_set(&get_subset(a, &self.n)), to_set(&get_subset(b, &self.n)), to_set(&get_subset(c, &self.n)), to_set(&get_subset(&ab_c, &self.n)), to_set(&get_subset(a, &self.n)), to_set(&get_subset(b, &self.n)), to_set(&get_subset(c, &self.n)),  to_set(&get_subset(&a_bc, &self.n)));
                        associativity=false;
                    }
            }
        }
    }
associativity
}
pub fn get_singleton(&self)->DMatrix<u32>{
    DMatrix::from_row_iterator(1, self.n as usize, (0..self.n).into_iter().map(|i|2u32.pow(i)))
}
}
impl Display for HyperGroupoidMat{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table:DMatrix<String>=DMatrix::from_iterator(self.n as usize, self.n as usize, 
            self.hyper_composition.iter().map(|x|format!("{:?}",to_set(&get_subset(x, &self.n)))));
        
        write!(f, "\nH: {:?},\nHypercomposition table:\n{} It is represented by: {}Size:{}\n", self.h, table, self.hyper_composition, self.n )
    }
}




/*DELATE HERE! EVERYTHING IS ON HYPERGROUPOIDMAT

just give a look at fix_associativity(), sometime it works*/

impl Display for HyperGroupoid{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table=String::new();
        //let s:Vec<Vec<HashSet<u32>>> =self.hyper_composition.iter().map(|x|x.iter().map(|y|to_set(&get_subset(&y, &self.n))).collect()).collect();
        for item in &self.hyper_composition{
            let set_item:Vec<HashSet<u32>>=item.iter().map(|x|to_set(&get_subset(x, &self.n))).collect();
            table.push_str(&format!("{:?}\n",set_item));
        }
        write!(f, "H: {:?},\nHypercomposition table:\n{}", self.h, table) 
    }
}


impl HyperGroupoid {
    pub fn new_random_from_cardinality(n:&u32)->Self{
        let h_vec=(0..*n as u32).into_iter().collect();
        let ht = get_random_hypercomposition_table(&n);
        HyperGroupoid { 
            h: h_vec, 
            hyper_composition: ht, 
            n: *n}
    }
    pub fn get_subset_from_k(&self,k:&u32)->HashSet<u32>{
        /*
        k is a number in 0..2^n-1. We use its binary representation to build a set
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

    pub fn elements_multiplication(&self,k:&u32,l:&u32)->u32{
        if k>=&self.n||l>=&self.n {panic!("k and l must be elements of H. Size of H is {}, k is {} and l is {}",self.n,k,l)}
        let k_subset = to_set(&vec![*k]);
        let l_subset = to_set(&vec![*l]);
        self.subsets_multiplication(&k_subset, &l_subset)

    }
    pub fn subsets_multiplication(&self, k_subset:&HashSet<u32>,l_subset:&HashSet<u32>)->u32{
        if !k_subset.is_subset(&self.h)||!l_subset.is_subset(&self.h) { panic!("K must be a subset of H!")};
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
        indexes.iter().fold(0u32, |acc, x| acc|self.hyper_composition[x.0 as usize][x.1 as usize])

    }
    pub fn get_singleton(&self)->Vec<HashSet<u32>>{
        //self.h.iter().map(|x| vec![*x]).collect()
        let mut singletons = Vec::new();
        let sort_h:Vec<_>=self.h.clone().into_iter().sorted().collect();
        for item in sort_h {
            singletons.push(HashSet::from([item]));
        }
        singletons

        
    }
    pub fn is_associative(&self)->bool{
        let mut associativity:bool=false;
        for a in self.get_singleton() {
            for b in self.get_singleton(){
                for c in self.get_singleton(){
                    let ab_c=self.subsets_multiplication(&to_set(
                        &get_subset(
                            &self.subsets_multiplication(&a, &b),&self.n)),&c);
                    let a_bc = self.subsets_multiplication(&a, &to_set(&get_subset(&self.subsets_multiplication(&b, &c),&self.n)));
                if a_bc!=ab_c{
                    println!("a is {:?},b is {:?},c is {:?}, ab_c is {:?}, a_bc is {:?}",a,b,c,ab_c,a_bc);
                    return false}
                    else{associativity= true}
                }
            }
        }
        associativity
    }
    pub fn assert_associativity(&self) {
        for a in self.get_singleton() {
            for b in self.get_singleton(){
                for c in self.get_singleton(){
                    let ab_c=self.subsets_multiplication(&to_set(
                        &get_subset(
                            &self.subsets_multiplication(&a, &b),&self.n)),&c);
                    let a_bc = self.subsets_multiplication(&a, &to_set(&get_subset(&self.subsets_multiplication(&b, &c),&self.n)));
                assert_eq!(a_bc,ab_c, "Not an associative hyperstructure!\na(bc) = {:?}\n(ab)c= {:?})",to_set(&get_subset(&a_bc,&self.n)),to_set(&get_subset(&ab_c, &self.n)));
                }
            }
        }
    }
    pub fn check_associativity (&self){
        for a in 0..self.n{
            for b in 0..self.n{
                for c in 0..self.n{
                    let subset_a =to_set(&get_subset(&a, &self.n));
                    let subset_b =to_set(&get_subset(&b, &self.n));
                    let subset_c =to_set(&get_subset(&c, &self.n));
                    let ab= self.subsets_multiplication(&subset_a, &subset_b);
                    let ab_c=self.subsets_multiplication(&to_set(
                        &get_subset(
                            &ab,&self.n)),&subset_c);
                    let bc= self.subsets_multiplication(&subset_b, &subset_c);
                    let a_bc = self.subsets_multiplication(&subset_a, &to_set(&get_subset(&bc,&self.n)));
                if ab_c!=a_bc{
                    println!("entered in the if");
                    println!("a is {:?} b is {:?} c is {:?}, ab is {:?}, bc is {:?},ab_c is {:?},a_bc is {:?}",subset_a,subset_b,subset_c,to_set(&get_subset(&ab, &self.n)),to_set(&get_subset(&bc, &self.n)),to_set(&get_subset(&ab_c, &self.n)),to_set(&get_subset(&a_bc, &self.n)));
                }
                }
            }
        }
    }
    pub fn fix_associativity(&self)->Self{
        // TO BE FIXEDDDDDDD
        if self.is_associative() {return  self.clone();}
        let mut rng = rand::thread_rng();
        let mut new_hyperstructure = self.clone();
        let mut new_hypercomposition: Vec<Vec<u32>>=self.hyper_composition.clone();
        for a in self.get_singleton() {
            for b in self.get_singleton(){
                for c in self.get_singleton(){
                    let mut ab= new_hyperstructure.subsets_multiplication(&a, &b);
                    let ab_c=new_hyperstructure.subsets_multiplication(&to_set(
                        &get_subset(
                            &ab,&self.n)),&c);
                    let bc=new_hyperstructure.subsets_multiplication(&b, &c);
                    let a_bc = new_hyperstructure.subsets_multiplication(&a, &to_set(&get_subset(&bc,&self.n)));
                    if ab_c!=a_bc{
    

                    }
                /* if ab_c<a_bc {
                    println!("ab_c<a_bc");
                    let missing_elements=a_bc-ab_c;
                    let ab_ones=ones_positions(ab, &(self.n as usize));
                    if ab==0{
                        let index_a:u32=a.iter().sum();
                        let index_b:u32=b.iter().sum();
                        new_hyperstructure.hyper_composition[index_a as usize][index_b as usize]=missing_elements;
                    } else {
                        let c_ones=ones_positions(subset_as_u32(&c), &(self.n as usize));
                        let mut indexes:Vec<(u32,u32)>=Vec::new();
                            for a in ab_ones{
                                for b in &c_ones{
                                    indexes.push((a,*b));
                                }
                            }
                        let value_to_be_modified :Vec<((u32,u32),u32)>= indexes.iter().map(|x| (*x,new_hyperstructure.hyper_composition[x.0 as usize][x.1 as usize])).collect();
                        let min_value=value_to_be_modified.iter().min_by_key(|x|x.1);
                        new_hyperstructure.hyper_composition[min_value.unwrap().0.0 as usize][min_value.unwrap().0.1 as usize]|=missing_elements;  }
                        }
                if a_bc<ab_c {
                    println!("a_bc<ab_c");
                    let missing_elements=ab_c-a_bc;
                    let bc_ones=ones_positions(bc, &(self.n as usize));
                    if bc==0{
                        let index_b:u32=b.iter().sum();
                        let index_c:u32=c.iter().sum();
                        new_hyperstructure.hyper_composition[index_b as usize][index_c as usize]=missing_elements;
                    } else {
                        let a_ones=ones_positions(subset_as_u32(&a), &(self.n as usize));

                        let mut indexes:Vec<(u32,u32)>=Vec::new();
                            for a in a_ones{
                                for b in &bc_ones{
                                indexes.push((a,*b));
                            }
                        }
                    let value_to_be_modified :Vec<((u32,u32),u32)>= indexes.iter().map(|x| (*x,new_hyperstructure.hyper_composition[x.0 as usize][x.1 as usize])).collect();
                    let min_value=value_to_be_modified.iter().min_by_key(|x|x.1);
                    new_hyperstructure.hyper_composition[min_value.unwrap().0.0 as usize][min_value.unwrap().0.1 as usize]|=missing_elements;          }
                } */
                }
            }
        }
        new_hyperstructure

    }
    pub fn fix_reproductivity(&self)->Self{
        if self.is_reproductive() {return self.clone();}
        let h=subset_as_u32(&self.h);
        let mut new_hypercomposition: Vec<Vec<u32>>=self.hyper_composition.clone();
        for element in self.get_singleton() {
            let x_h=self.subsets_multiplication(&element, &self.h);
            if x_h!=h{
                
                let row_index: &u32=&element.into_iter().sum();//this is the element converted in a singleton.
                let missing_integer = h-x_h;
                let min_value_in_row=self.hyper_composition[*row_index as usize].iter().min().unwrap();
                let column_position=self.hyper_composition[*row_index as usize].iter().position(|x| x==min_value_in_row).unwrap();
                new_hypercomposition[*row_index as usize][column_position]=min_value_in_row|missing_integer;
            }
        }
        let modified_hyperstructure  =HyperGroupoid{
            h:self.h.clone(),
            hyper_composition:new_hypercomposition,
            n:self.n
        };
    let mut transpose: Vec<_>  =(0..modified_hyperstructure.hyper_composition[0].len())
                                .map(|i| modified_hyperstructure.hyper_composition.iter().map(|inner| inner[i].clone()).collect::<Vec<u32>>())
                                .collect();
    for element in self.get_singleton() {
        let h_x=modified_hyperstructure.subsets_multiplication(&self.h, &element);
        if h_x!=h{
            
            let row_index: &u32=&element.into_iter().sum();//this is the element converted in a singleton.
            let missing_integer = h-h_x;
            let min_value_in_row=transpose[*row_index as usize].iter().min().unwrap();
            let column_position=transpose[*row_index as usize].iter().position(|x| x==min_value_in_row).unwrap();

            transpose[*row_index as usize][column_position]=min_value_in_row|missing_integer;
        }
    }
    new_hypercomposition=(0..transpose[0].len())
                        .map(|i| transpose.iter().map(|inner| inner[i].clone()).collect::<Vec<u32>>())
                        .collect();  
    HyperGroupoid {h: self.h.clone(), 
                    hyper_composition: new_hypercomposition, 
                    n: self.n }     
       }
        
    pub fn is_reproductive(&self)->bool{
        let h=subset_as_u32(&self.h);
        let mut reproductive = false;
        for element in self.get_singleton() {
            let x_h=self.subsets_multiplication(&element, &self.h);
            let h_x=self.subsets_multiplication(&self.h, &element);
            if(x_h==h)&&(h_x==h)
                {reproductive=true}
                else {return false;}
        }
        reproductive
    }
    pub fn assert_reproductivity(&self){
        let h=subset_as_u32(&self.h);
        for element in self.get_singleton() {
            let x_h=self.subsets_multiplication(&element, &self.h);
            let h_x=self.subsets_multiplication(&self.h, &element);
            assert_eq!(x_h,h,"Not a reproductive hyperstructure:\nxH = {:?}, H = {:?}",to_set(&get_subset(&x_h, &self.n)),self.h);
            assert_eq!(h_x,h,"Not a reproductive hyperstructure:\nHx = {:?}, H = {:?}",to_set(&get_subset(&h_x, &self.n)),self.h);
    }
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
pub fn collect_hypergroupoid(cardinality:&u32)->Vec<u128>{
    //TO BE FIXED
    let size = cardinality.pow(3);
    let x = 2u128.pow(size-cardinality);
    let y = 2u128.pow(size);
    
        println!("size is {size}");
        println!("there are {} to be tested", y-x);

    (2u128.pow(size-cardinality)..2u128.pow(size)).into_iter().filter(|i|representing_hypergroupoid(&mut i.clone(),&cardinality)).collect()
        
    

}

