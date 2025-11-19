use core::panic;
use std::{collections::HashSet};
use nalgebra::{DMatrix, SimdRealField, iter};
use itertools::Itertools;
use permutation::Permutation;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator as _, ParallelIterator};
use crate::hs::hypergroupoids::HyperGroupoid;
use crate::hs::HyperStructureError;
use crate::utilities::{n_to_binary_vec, support, write};

use crate::{binary_relations::relation_matrix::RelationMatrix, utilities::{permutation_matrix_from_permutation, representation_permutation_subset}};

#[derive(Debug,Clone,PartialEq, Eq,)]
pub struct Relation {
    pub a: HashSet<u64>,
    pub b: HashSet<u64>,
    pub rel: Vec<(u64,u64)>,
}
impl PartialOrd for Relation {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
    
        self.rel.partial_cmp(&other.rel)
    }
}
impl Ord for Relation {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.rel.cmp(&other.rel)
    }
}
impl Relation {
    pub fn new_from_elements(cardinality:&u64,rel : &Vec<(u64,u64)>)->Self{
        let a: HashSet<u64>= (0..*cardinality).collect();
        let b: HashSet<u64> = (0..*cardinality).collect();
        Relation { a, b, rel: rel.to_vec() }
    }
/// Return an equivalence associated to a preorder.
/// If r is a preorder on X, then the relation k defined by 
/// `x k y if and only if x r y and y r x` is an equivalence on X.
/// 
/// #Example
///```
/// use hyperstruc::binary_relations::relations::Relation;
/// 
/// let cardinality = 3u64;
/// let r = [
///     (0,0),(0,1),
///           (1,1),
///                (2,2)
/// ].to_vec();
/// let preorder = Relation::new_from_elements(&cardinality,&r);
/// assert!(preorder.is_pre_order());
/// let e  = Relation::new_from_preorder(&preorder);
/// 
    pub fn new_from_preorder(preorder:&Self)->Result<Self,HyperStructureError>{
        match preorder.is_pre_order() {
            true => {
                let n =  preorder.a.len() as u64;
                let rel = (0..n).cartesian_product(0..n).into_iter().filter(|(a,b)|preorder.rel.contains(&(*a,*b))&&preorder.rel.contains(&(*b,*a))).collect_vec();
                Ok(Relation::new_from_elements(&n, &rel))
            },
            false => Err(HyperStructureError::NotPreOrder),
        }
    }
    pub fn is_reflexive(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        self.diagonal().rel.iter().all(|x|self.rel.contains(x))
    }
    pub fn is_symmetric(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        self.rel.iter().all(|(x,y)|self.rel.contains(&(*y,*x)))
    }
/// Checks whether the relation is transitive.
///
/// A relation `R` on a set `A` is transitive if for all `a, b, c ∈ A`,
/// whenever `(a, b) ∈ R` and `(b, c) ∈ R`, then `(a, c)` must also be in `R`.
///
/// # Examples
/// ```
/// use hyperstruc::binary_relations::relations::Relation;
/// use std::collections::HashSet;
/// 
/// let a:HashSet<u64>=(1..=3).into_iter().collect();
/// let rel = vec![(1, 1), (1, 2), (2, 2), (2, 3), (1, 3)];
/// let trans_rel = Relation{
/// a: a.clone(),
/// b: a,
/// rel: rel};
/// 
/// assert!(trans_rel.is_transitive());
///
/// ```
///
/// # Panics
/// Panics if the domain and codomain of the relation are not equal,
/// as transitivity is only defined for homogeneous relations.
///
    pub fn is_transitive(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        self.rel.iter()
            .all(|&(a, b)| {
                self.rel.iter().all(|&(c, d)| 
                    b != c || self.rel.contains(&(a, d))
                    )
                }
            )
    }
    pub fn diagonal(&self)->Self{
        assert_eq!(self.a,self.b);
        let id= DMatrix::identity(
            self.a.len(), 
            self.a.len());
        RelationMatrix(id).into_relation()
    }
/// Checks whether the relation is an **equivalence relation**.
///
/// A relation is an equivalence relation if it satisfies:
/// - **Reflexivity**: For every `x` in the set, `(x, x)` is in the relation.
/// - **Symmetry**: For every `(x, y)` in the relation, `(y, x)` is also in the relation.
/// - **Transitivity**: For all `(x, y)` and `(y, z)` in the relation, `(x, z)` is also in the relation.
///
/// # Returns
///
/// `true` if the relation is reflexive, symmetric, and transitive; otherwise `false`.
///
/// # Example
///
/// ```
/// use std::collections::HashSet;
/// use hyperstruc::binary_relations::relations::Relation;
///
/// let a: HashSet<u64> = [0,1,2].into_iter().collect();
/// let rel: Vec<(u64, u64)> = vec![
///     (0,0), (1,1), (2,2),
///     (0,1), (1,0),
///     (1,2), (2,1),
///     (0,2), (2,0)
///     ];
///
/// let r = Relation { a: a.clone(), b: a.clone(), rel };
/// assert!(r.is_equivalence());
/// 
/// ```
    pub fn is_equivalence(&self)->bool{
        self.is_reflexive()&&self.is_symmetric()&&self.is_transitive()
    }
/// Return true if the relation is regular on a hypercompositional structure
/// 
/// #Example 
///```
/// use hyperstruc::{binary_relations::relations::{self, Relation}, hg_3::tag_hypergroups_3::TAGS_HG_3, hs::hypergroups::HyperGroup};
/// use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
/// 
/// let cardinality =3u64;
/// TAGS_HG_3.par_iter().all(|tag| {
///     let hg = HyperGroup::new_from_tag_u128(tag,&cardinality);
///     hg.beta_relation().is_left_regular(&hg.0)
/// });
/// 

    pub fn is_left_regular(&self, h:&HyperGroupoid)->bool{
        match self.is_equivalence() {
            true => self.rel.iter().all(|(x,y)| 
            (0..h.n).into_iter().all(|a|
                {
                let ax = h.mul_by_representation(&(1<<a), &(1<<x));
                let ay = h.mul_by_representation(&(1<<a), &(1<<y));
                support(&ax, &h.n).iter().all(|u|
                    support(&ay, &h.n).iter().any(|v|
                        self.rel.contains(&(*u as u64,*v as u64))
                        )
                    )
                })
            ),
            false => panic!("The relation is not left regular. In particular, input relation is not an equivalence."),
        }
        
    }
    pub fn is_right_regular(&self, h:&HyperGroupoid)->bool{
        match self.is_equivalence() {
            true => self.rel.iter().all(|(x,y)| 
            (0..h.n).into_iter().all(|a|
                {
                let xa = h.mul_by_representation(&(1<<x), &(1<<a));
                let ya = h.mul_by_representation(&(1<<y), &(1<<a));
                support(&xa, &h.n).iter().all(|u|
                    support(&ya, &h.n).iter().any(|v|
                        self.rel.contains(&(*u as u64,*v as u64))
                        )
                    )
                })
            ),
            false => panic!("The relation is not right regular. In particular, input relation is not an equivalence."),
        }
            }
    pub fn is_regular(&self,h:HyperGroupoid)->bool{
        self.is_left_regular(&h)&&self.is_right_regular(&h)
    }
    pub fn get_class(&self, x:&u64)->(u64,Vec<u64>) {
        let class:Vec<u64> = self.a.iter()
            .filter(|y|self.are_in_relations(x, y))
            .map(|x|*x)
            .sorted()
            .collect();
        let representant = class.iter().min();
        (*representant.unwrap(),class)
    }
    pub fn are_in_relations(&self,x:&u64,y:&u64)->bool {
        self.rel.contains(&(*x,*y)) 
    }
    pub fn quotient_set(&self)->Vec<(u64,Vec<u64>)>{
    assert!(self.is_equivalence(), "Relation {:?} is not an equivalence!",self.rel);
    self.a.iter().map(|x|
        self.get_class(x)
        ).sorted().unique().collect()
    }
    pub fn zero_one_matrix(&self)-> RelationMatrix {
        let generating_function  =|a: usize,b: usize| {
            match self.rel.contains(&(a as u64,b as u64)) {
            true =>1u64,
            false => 0u64
        }
        };

        let rel_mat = DMatrix::from_fn(self.a.len(), self.b.len(), generating_function);
        RelationMatrix(rel_mat)
    }
    pub fn reflexive_closure(&self)->Self { 
        (self.zero_one_matrix()|self.diagonal().zero_one_matrix()).into_relation()
    }
    pub fn symmetric_closure(&self)->Self { 
        (self.zero_one_matrix()|self.zero_one_matrix().transpose_relation()).into_relation()
    }
    pub fn is_pre_order(&self)->bool{
        self.is_reflexive()&&self.is_transitive()
    }
    pub fn is_anti_symmetric(&self)->bool{
        todo!()

    }
///
/// Compute the transitive closure of a relation on a set H
/// 
/// #Example
/// ```
/// use hyperstruc::binary_relations::relations::Relation;
/// use std::collections::HashSet;
/// 
/// let a: HashSet<u64> = (0..3).into_iter().collect();
/// let rel: Vec<(u64, u64)> = vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)];
///  let relation  = Relation{
///     a: a.clone(),
///     b: a.clone(),
///     rel,
/// };
/// let cl_trans = relation.transitive_closure();
/// assert!(cl_trans.is_transitive());
/// 
/// 
    pub fn transitive_closure(&self)-> Self {
        if self.is_equivalence() {return self.clone();}
    //This is Algorithm 3.3
    let mut x  = self.clone();
    for _ in 2..=self.a.len() {
        for (a,b) in x.rel.clone(){
            for (c,d) in self.rel.clone(){
                if (b==c) && !x.rel.contains(&(a,d)) {
                    x.rel.push((a,d));
                    x.rel.sort();
                }
            }
        }
    }
    x
}
///
/// Compute the transitive closure of a relation on a set H
/// 
/// #Example
/// ```
/// use hyperstruc::binary_relations::relations::Relation;
/// use std::collections::HashSet;
/// 
/// let a: HashSet<u64> = (0..3).into_iter().collect();
/// let rel: Vec<(u64, u64)> = vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)];
///  let relation  = Relation{
///     a: a.clone(),
///     b: a.clone(),
///     rel,
/// };
/// let cl_trans = relation.transitive_closure();
/// assert!(cl_trans.is_transitive());
/// 
///
pub fn transitive_closure_warshall(&self)-> Self {
        if self.is_equivalence() {return self.clone();}

    let mut r = self.zero_one_matrix();
    for k in 0..r.0.ncols() {
        let col_k  =r.0.column(k);
        let row_k = r.0.row(k);
        let mul = col_k*row_k;
        let mul_rel = RelationMatrix(mul);
        r = r|mul_rel;
    }
    r.into_relation()
}
pub fn isomorphic_relation_from_permutation(&self,sigma: &Permutation)->Self{

    let perm_mat= permutation_matrix_from_permutation(&(self.a.len() as u64), &sigma.clone());
    let isomorphic_matrix=perm_mat.clone()*self.zero_one_matrix().0*perm_mat.transpose();
    RelationMatrix(isomorphic_matrix).into_relation()

    }
pub fn collect_isomorphism_class(&self)->(Relation,Vec<Relation>) {
    //We collect the isomorphism class of a relation on a set A.
    assert_eq!(self.a,self.b);

    let cardinality = self.a.len();
    let permutation_vec:Vec<Vec<usize>> = (0..cardinality).permutations(cardinality ).collect();
    let permutation:Vec<Permutation> = permutation_vec.par_iter()
        .map(|sigma| 
                Permutation::oneline(sigma.clone())
            ).collect();
    let mut isomorphism_classes: Vec<Relation>=permutation.par_iter()
        .map(|sigma|
                self.isomorphic_relation_from_permutation(&sigma).zero_one_matrix().into_relation()
            )
            .collect();
        isomorphism_classes.sort();
        isomorphism_classes.dedup();
   
    let representant_of_class= isomorphism_classes[0].clone();// As isomorphism_class is sorted, this is the minimum in the class
       (representant_of_class,isomorphism_classes)
    }
/* pub fn permutation_of_table(&self,sigma:&Permutation)->Self{
    let permutation_rel = &self.zero_one_matrix();
    let n = self.a.len();
    let alfa =DMatrix::from_iterator(
        n , 
        n, 
        permutation_rel.0
                .iter()
                .map(|x| 
                        sigma.apply_idx(*x as usize) as u64)
                    );
    RelationMatrix(alfa).into_relation()    
} */
pub fn is_isomorphic_to(&self,other:&Relation)->bool {
    self.collect_isomorphism_class().1.contains(&other)
}
}

pub fn enumerate_reflexive_relation(cardinality:&usize)->Vec<RelationMatrix>{

        // The number of reflexive relation from a set of `n` elements is 
        // 2^{n^2-n}.
        let free_positions:Vec<(usize,usize)> = (0..*cardinality)
            .cartesian_product(0..*cardinality)
            .filter(|(i,j)| i!=j).collect();
        let free = cardinality.pow(2u32)-cardinality;
        let total = 1usize<<free;

        (0..total)
        .into_par_iter()
        .map(|mask| {
            // Start with identity (reflexive diagonal = 1)
            let mut r = DMatrix::<u64>::identity(*cardinality, *cardinality);

            // Set off-diagonal entries according to mask bits
            for (bit_index, &(i, j)) in free_positions.iter().enumerate() {
                let bit = ((mask >> bit_index) & 1) as u64;
                r[(i, j)] = bit;
            }

            RelationMatrix(r)
        })
        .collect()
    

    }
    pub fn par_reflexive_relations(cardinality: &u64) -> impl ParallelIterator<Item = RelationMatrix> + use<'_> {
    let free_positions: Vec<(usize, usize)> = (0..*cardinality as usize)
        .cartesian_product(0..*cardinality as usize)
        .filter(|(i, j)| i != j)
        .collect();

    let free = cardinality * cardinality - cardinality;
    let total = 1usize << free;

    // Parallel iterator that generates each matrix independently
    (0..total).into_par_iter().map(move |mask| {
        let mut r = DMatrix::<u64>::identity(*cardinality as usize, *cardinality as usize);
        for (bit_index, &(i, j)) in free_positions.iter().enumerate() {
            let bit = ((mask >> bit_index) & 1) as u64;
            r[(i, j)] = bit;
        }
        RelationMatrix(r)
    })
}
pub fn pre_order_enumeration(cardinality:&u64)-> impl ParallelIterator<Item = RelationMatrix> + use<'_>{
    par_reflexive_relations(cardinality).into_par_iter().filter(|x|x.into_relation().is_transitive())
}
pub fn enumerate_pre_order_classes(cardinality: &u64)->(Vec<Vec<(Relation,Vec<Relation>)>>,Vec<usize>){
    let mut classes:Vec<(Relation,Vec<Relation>)> = pre_order_enumeration(cardinality)
        .map(|x|
            x.into_relation().collect_isomorphism_class()
        ).collect();
    classes.sort();
    classes.dedup();

    let permutation_len = (0..*cardinality).permutations(*cardinality as usize).count();

    let c_k =
        (1..=permutation_len).into_iter()
            .map(|k|
                classes.iter()
                .filter(
                    |x|
                    x.1.len()==k)
                    .sorted()
                    .map(|s|s.clone())
                    .collect()
        ).collect::<Vec<Vec<(Relation, Vec<Relation>)>>>();
    let c:Vec<usize>= c_k.iter().map(|x: &Vec<(Relation, Vec<Relation>)>|x.len()).collect();
    let s:String = c_k.iter().map(|x|format!("{:?}\n",x)).collect();
    let _ = write(s.clone(),&format!("enumeration_preorder_{}",cardinality));
    
    (c_k,c)
}