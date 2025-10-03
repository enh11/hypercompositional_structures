use std::{collections::HashSet};
use nalgebra::DMatrix;
use itertools::Itertools;

use crate::binary_relations::relation_matrix::RelationMatrix;

#[derive(Debug,Clone,PartialEq, Eq)]
pub struct Relation {
    pub a: HashSet<u64>,
    pub b: HashSet<u64>,
    pub rel: Vec<(u64,u64)>,
}
impl Relation {
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
            true =>1u8,
            false => 0u8
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
    //println!("r^{i} = {:?}",x.rel);
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
}