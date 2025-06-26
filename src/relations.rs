use std::{collections::HashSet, ops::Mul};

use itertools::Itertools;

#[derive(Debug,Clone,PartialEq, Eq)]
pub struct Relation {
    pub a: HashSet<u64>,
    pub b: HashSet<u64>,
    pub rel: Vec<(u64,u64)>
}
impl <'a,'b>Mul<&'a Relation> for &'b Relation{
    type Output = Relation;
    fn mul(self, rhs: &'a Relation) -> Self::Output {
        assert_eq!(self.b,self.a,"product not define!");
        let pairs:Vec<(&u64, &u64)> =self.a.iter().cartesian_product(rhs.b.iter()).collect();
        let pairs=pairs.iter()
            .filter(|(x,z)|
                    self.b.iter()
                    .any(|y|self.rel.contains(&(**x,*y))&&rhs.rel.contains(&(*y,**z))))
                    .map(|p|(*p.0,*p.1)
                ).sorted().unique().collect();
    Relation{
        a:self.a.clone(),
        b:rhs.b.clone(),
        rel:pairs
    }
    }
}
impl Relation {
    pub fn is_reflexive(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        let pairs=self.a.iter().zip(self.b.clone())
        .map(|(x,y)|(*x,y))
            .collect_vec();
        pairs.iter().all(|x|self.rel.contains(x))

    }
    pub fn is_symmetric(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        self.rel.iter().all(|(x,y)|self.rel.contains(&(*y,*x)))
    }
    pub fn is_transitive(&self)->bool{
        assert_eq!(self.a,self.b,"Domain and codomain not coincede!");
        let rr=(self*self).rel;
        
            rr.iter().all(|x|self.rel.contains(x))
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
    pub fn get_class(&self, x:u64)->(u64,Vec<u64>) {
        let class:Vec<_> = self.a.iter().filter(|y|self.are_in_relations(x, **y)).map(|x|*x).sorted().collect();
        let representant = class.iter().min();        (*representant.unwrap(),class)
    }
    /// Checks if two elements `x` and `y` are directly related in the relation `R`,
    /// or if there exists a third element `z` such that `(x, z)` and `(z, y)` are in the relation.
    ///
    /// This is a **partial transitivity check**, useful when computing or testing transitive closure manually.
    ///
    /// # Arguments
    ///
    /// * `x` - The first element.
    /// * `y` - The second element.
    ///
    /// # Returns
    ///
    /// `true` if `(x, y)` ∈ R or if ∃z such that `(x, z)` ∈ R and `(z, y)` ∈ R.
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashSet;
    /// use hyperstruc::relations::Relation;
    ///
    /// let a: HashSet<u64> = [1, 2, 3].into_iter().collect();
    /// let rel: Vec<(u64,u64)> = vec![(1, 2), (2, 3)];
    /// let r = Relation { a: a.clone(), b: a.clone(), rel };
    ///
    /// assert_eq!(r.are_in_relations(1, 2), true); // directly in rel
    /// assert_eq!(r.are_in_relations(1, 3), true); // via 2 (1→2, 2→3)
    /// assert_eq!(r.are_in_relations(2, 1), false); // not related
    /// ```
    pub fn are_in_relations(&self,x:u64,y:u64)->bool {

        if self.rel.contains(&(x,y)) {return true}
        if self.a.iter().any(|z|self.rel.contains(&(x,*z))&&self.rel.contains(&(*z,y))) {return true;}
        false
    }
    pub fn quotient_set(&self)->Vec<(u64,Vec<u64>)>{
    assert!(self.is_equivalence(), "Relation is not an equivalence!");
    self.a.iter().map(|x|self.get_class(*x)).sorted().unique().collect()
    }
}