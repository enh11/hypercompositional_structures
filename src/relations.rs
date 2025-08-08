use std::collections::HashSet;

use itertools::Itertools;

#[derive(Debug,Clone,PartialEq, Eq)]
pub struct Relation {
    pub a: HashSet<u64>,
    pub b: HashSet<u64>,
    pub rel: Vec<(u64,u64)>
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
/// use hyperstruc::relations::Relation;
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
        self.rel.iter().all(|&(a, b)| {
            self.rel.iter().all(|&(c, d)| 
                b != c || self.rel.contains(&(a, d))
            )
})
    }
    pub fn diagonal(&self)->Self{
        assert_eq!(self.a,self.b);
        let diag:Vec<(u64,u64)> = 
                self.a.iter()
                    .zip(self.b.iter())
                    .map(|(x,y)| (*x,*y)).collect();
        Relation { a: self.a.clone(), b: self.b.clone(), rel: diag}

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
/// use hyperstruc::relations::Relation;
///
/// let a: HashSet<u64> = [1, 2, 3].into_iter().collect();
/// let rel: Vec<(u64, u64)> = vec![
///     (1, 1), (2, 2), (3, 3),
///     (1, 2), (2, 1),
///     (2, 3), (3, 2),
///     (1, 3), (3, 1)
/// ];
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
    assert!(self.is_equivalence(), "Relation is not an equivalence!");
    self.a.iter().map(|x|self.get_class(x)).sorted().unique().collect()
    }
}