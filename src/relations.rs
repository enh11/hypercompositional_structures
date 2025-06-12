use std::{collections::HashSet, ops::Mul};

use itertools::{sorted, Itertools};

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
            .into_iter()
            .map(|(x,y)|(2u64.pow(*x as u32),2u64.pow(y as u32)))
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

    pub fn is_equivalence(&self)->bool{
        self.is_reflexive()&&self.is_symmetric()&&self.is_transitive()
    }
    pub fn get_class(&self, x:u64)->(u64,Vec<u64>) {
        let class:Vec<_> = self.a.iter().filter(|y|self.are_in_relations(x, **y)).map(|x|1<<*x).sorted().collect();
        let representant = class.iter().min().unwrap();
        (*representant,class)
    }
    pub fn are_in_relations(&self,x:u64,y:u64)->bool {
        let singletons = self.a.iter().map(|x|1<<x).sorted().collect_vec();
        let singleton_x = 1<<x;
        let singleton_y = 1<<y;
        if self.rel.contains(&(singleton_x,singleton_y)) {return true}
        if singletons.iter().any(|z|self.rel.contains(&(singleton_x,*z))&&self.rel.contains(&(*z,singleton_y))) {return true;}
        false
    }
    pub fn collect_classes(&self)->Vec<(u64,Vec<u64>)>{
    assert!(self.is_equivalence(), "Relation is not an equivalence!");
    self.a.iter().map(|x|self.get_class(*x)).sorted().unique().collect()
    }
}