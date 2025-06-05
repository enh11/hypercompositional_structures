use itertools::Itertools;

use crate::{hs::{circumference_radius_d_filtered, hg_in_circumference_radius_one, HyperGroupoid}, hypergroups::{HyperGroup, HyperStructureError}, utilities::{get_complement_subset, ones_positions, U1024}};

#[derive(Debug, Copy, Clone)]
pub struct UnitInterval(f64);

impl UnitInterval {
    pub fn new(val: f64) -> Option<Self> {
        if (0f64<=val) && (val<=1f64) {
            Some(UnitInterval(val))
        } else {
            None
        }
    }

    pub fn get(self) -> f64 {
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
        let function  = |a:u64,b:u64| 
            (0..self.h.n).into_iter().filter(|z|
                (self.mu)(*z).unwrap().0>=(self.mu)(a).unwrap().0.min((self.mu)(b).unwrap().0)
                &&
                (self.mu)(*z).unwrap().0<=(self.mu)(a).unwrap().0.max((self.mu)(b).unwrap().0)
).into_iter().fold(0,|acc,x| acc|1<<x);
HyperGroup::new_from_function(function, &self.h.n)
    }
}

impl HyperGroupoid {
    pub fn get_Q_u(&self,u:&u64)->Vec<(u64, u64)>{
        let u_singleton: u64 = 1u64<<u;
        (0..self.n).into_iter().cartesian_product(0..self.n)
                    .filter(|(x,y)|
                        u_singleton&self.mul_by_representation(&(1<<x), &(1<<y))==u_singleton).collect()

    }
    pub fn get_q_u(&self,u:&u64)->usize{
       self.get_Q_u(u).len()

    }
    pub fn get_alpha_u(&self,u:&u64)->f64{
        self.get_Q_u(u).into_iter().map(|(x,y)|1f64/(ones_positions(&self.mul_by_representation(&(1u64<<x), &(1u64<<y)), &self.n).len() as f64)).sum()
    }
    pub fn get_mu_u(&self,u:&u64)->f64{
        self.get_alpha_u(u)/(self.get_q_u(u) as f64)
    } 
    pub fn get_corsini_membership_function(&self) -> MembershipFunction {
        let h = self.clone();  // or use Arc<Self> to avoid cloning big data
        Box::new(move |u: u64| UnitInterval::new(h.get_mu_u(&u)))
    }
}