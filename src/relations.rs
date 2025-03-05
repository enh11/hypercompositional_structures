use std::{collections::HashSet, ops::Mul, os::unix::process};

use itertools::Itertools;

#[derive(Debug)]
pub struct Relation {
    pub a: HashSet<u32>,
    pub b: HashSet<u32>,
    pub rel: Vec<(u32,u32)>
}
impl <'a,'b>Mul<&'a Relation> for &'b Relation{
    type Output = Relation;
    fn mul(self, rhs: &'a Relation) -> Self::Output {
        assert_eq!(self.b,self.a,"product not define!");
        let pairs:Vec<(&u32, &u32)> =self.a.iter().cartesian_product(rhs.b.iter()).collect();
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
            .map(|(x,y)|(2u32.pow(*x),2u32.pow(y)))
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
   pub fn collect_classes(&self)->Vec<(u32,Vec<u32>)>{
    assert!(self.is_equivalence());
    let a:Vec<u32>= self.a.iter().sorted().map(|x|2u32.pow(*x)).collect();
    let b:Vec<u32>= self.b.iter().sorted().map(|x|2u32.pow(*x)).collect();

        let mut processed:Vec<u32>=Vec::new();
        let mut classes:Vec<(u32,Vec<u32>)>= Vec::new();
        
        for representant in &a {
            if processed.contains(representant){continue;}
            let mut class:Vec<u32>=Vec::new();
                for element in &b {
                    if self.rel.iter().contains(&(*representant,*element)){
                        class.push(*element);
                        processed.push(*element);
                    }
                
                }
                classes.push((*representant,class.clone()));
            }

        
        classes
    }
}