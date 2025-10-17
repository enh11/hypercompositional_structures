//!
//! In this module we implement FuzzySubsets for HyperGroupoid.
//!  
//!  Examples and tests come from `Some Remarks on Hyperstructures their Connections with Fuzzy Sets and Extensions to Weak Structures` by `Piergiulio Corsini`.
//! 
use itertools::Itertools;
use num_rational::Rational64;
use num_traits::{One, Zero};

use crate::{hs::{HyperGroupoid}, hypergroups::{HyperGroup, HyperStructureError}};

#[derive(Debug, Copy, Clone)]
pub struct UnitInterval(Rational64);

impl UnitInterval {
    pub fn new(val: Rational64) -> Option<Self> {
        if (Rational64::zero()<=val) && (val<=Rational64::one()) {
            Some(UnitInterval(val))
        } else {
            None
        }
    }

    pub fn get(self) -> Rational64 {
        self.0
    }
}

type MembershipFunction = Box<dyn Fn(u64) -> Option<UnitInterval> + Send + Sync>;

pub struct FuzzySubset {
    pub h: HyperGroupoid,
    pub mu: MembershipFunction
}

impl FuzzySubset {
    pub fn new(h:HyperGroupoid,mu:MembershipFunction)->Self{
        FuzzySubset { h, mu }
    }
     pub fn get_corsini_join_space(&self)->Result<HyperGroup, HyperStructureError>{
        let function  = |a:usize,b:usize| 
            (0..self.h.n).into_iter().filter(|z|
                (self.mu)(*z).unwrap().0>=(self.mu)(a as u64).unwrap().0.min((self.mu)(b as u64).unwrap().0)
                &&
                (self.mu)(*z).unwrap().0<=(self.mu)(a as u64).unwrap().0.max((self.mu)(b as u64).unwrap().0)
).into_iter().fold(0,|acc,x| acc|1<<x);
HyperGroup::new_from_function(function, &self.h.n)
    }
}

impl HyperGroupoid {
    pub fn get_qq_u(&self,u:&u64)->Vec<(u64, u64)>{
        let u_singleton: u64 = 1u64<<u;
        (0..self.n).into_iter().cartesian_product(0..self.n)
                    .filter(|(x,y)|
                        u_singleton&self.mul_by_representation(&(1<<x), &(1<<y))==u_singleton).collect()

    }
/// 
/// Compute the value `q(u) = |Q(u)| = {(a, b)∈ HxH |u ∈ ab}`.
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// let cardinality =2u64;
/// let input_values  = vec![vec![0,1],vec![0,1],
///                            vec![1],vec![0]];
/// let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
/// 
/// let q_0=hs.get_q_u(&0);
/// assert_eq!(q_0,3);
/// 
/// let q_1=hs.get_q_u(&1);
/// assert_eq!(q_0,3);
/// 
/// 
    pub fn get_q_u(&self,u:&u64)->usize{
       self.get_qq_u(u).len()

    }
/// 
/// Compute the value `alpha(u) = sum_{(x,y)∈ Q(u)} 1/|xy|`, where `Q(u) = {(a, b)∈ HxH |u ∈ ab}`
/// 
/// # Example
/// ```
/// use hyperstruc::hs::HyperGroupoid;
/// use num_rational::Rational64;
/// 
/// let cardinality =2u64;
/// let input_values  = vec![vec![0,1],vec![0,1],
///                            vec![1],vec![0]];
/// let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);

/// 
/// let alpha_0=hs.get_alpha_u(&0);
/// assert_eq!(alpha_0,Rational64::from(2));
/// 
/// let alpha_1=hs.get_alpha_u(&1);
/// assert_eq!(alpha_1,Rational64::from(2));
/// 
/// 
    pub fn get_alpha_u(&self,u:&u64)->Rational64{
        self.get_qq_u(u).into_iter().map(|(x,y)|
            Rational64::new(1,self.mul_by_representation(&(1u64<<x), &(1u64<<y)).count_ones() as i64)).sum()
    }
    pub fn get_mu_u(&self,u:&u64)->Rational64{
         self.get_alpha_u(u)/Rational64::from(self.get_q_u(u) as i64)
    } 
    pub fn get_corsini_membership_function(&self) -> MembershipFunction {
        let h = self.clone();  // or use Arc<Self> to avoid cloning big data
        Box::new(move |u: u64| UnitInterval::new(h.get_mu_u(&u)))
    }
    pub fn get_fuzzy_grade(&self)->usize{
        let mut h_0:HyperGroup;
        let mut h_1:HyperGroup;
        let mut grade;
        if self.is_hypergroup(){
            h_0=HyperGroup::new_from_tag_u1024(&self.get_integer_tag_u1024(),&self.n);
            h_1=h_0.get_next_corsini_joinspace();
            grade=0;
        }
        else{
            h_0=self.get_corsini_fuzzysubset().get_corsini_join_space().unwrap();
            h_1 =h_0.get_next_corsini_joinspace();

            grade=1;
        }
        while !h_0.is_isomorphic_to(&h_1) {
            h_0=h_1;
            h_1=h_0.get_next_corsini_joinspace();
            grade+=1;            
        }
        grade
       
    }
    pub fn get_strong_fuzzy_grade(&self)->usize{
        let mut h_0:HyperGroup;
        let mut h_1:HyperGroup;
        let mut degree;
        if self.is_hypergroup(){
            h_0=HyperGroup::new_from_tag_u1024(&self.get_integer_tag_u1024(),&self.n);
            h_1=h_0.get_next_corsini_joinspace();
            degree=0;
        }
        else{
            h_0=self.get_corsini_fuzzysubset().get_corsini_join_space().unwrap();
            h_1 =h_0.get_next_corsini_joinspace();

            degree=1;
        }
        while h_0!=h_1 {
            h_0=h_1;
            h_1=h_0.get_next_corsini_joinspace();
            degree+=1;            
        }
        degree
       
    }

}
impl HyperGroup {
    pub fn get_next_corsini_joinspace(&self)->Self{
        self.get_corsini_fuzzysubset().get_corsini_join_space().unwrap()
}
}