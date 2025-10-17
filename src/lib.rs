//! This library contains utilities for Hyperstructures.
//!
#[macro_use]
pub mod hs;
pub mod hypergroups;
pub mod utilities;
pub mod hg_2;
pub mod hg_3;
pub mod binary_relations;
pub mod enumeration;
pub mod group_cayley_tables;
pub mod fuzzy;
pub mod generating_functions;
pub mod quotient_hg;
#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use nalgebra::DMatrix;
    use crate::{hs::HyperGroupoid, hypergroups::HyperGroup};
    #[test]
    fn aboutotorab_4_3() {
        let cardinality=4u64;
        let hs = HyperGroupoid::new_from_matrix(
                &DMatrix::from_row_slice(
                    cardinality as usize,
                    cardinality as usize,
                    &[6,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]));
        assert!(hs.is_associative());

        let ph = hs.collect_all_finite_hyperproducts().0;
        let expected_ph  = vec![1,2,4,6,8,10];
        assert_eq!(ph,expected_ph);

        let beta = hs.beta_relation();
        let expected_beta: Vec<(u64, u64)>  = vec![(0, 0), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (3, 1), (3, 3)];
        assert_eq!(beta.rel,expected_beta);
    }
    #[test]
    fn aboutotorab_4_4() {
        let cardinality=  3u64;
        let hs = HyperGroupoid::new_from_matrix(
            &DMatrix::from_row_slice(
                cardinality as usize, 
                cardinality as usize, 
                &[1,3,5,1,2,5,1,3,4]));
        assert!(hs.is_associative());

        let ph = hs.collect_all_finite_hyperproducts().0;
        let expected_ph  = vec![1,2,3,4,5];
        assert_eq!(ph,expected_ph);

        let beta = hs.beta_relation().rel;
        let expected_beta: Vec<(u64, u64)> = vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)];
        assert_eq!(beta,expected_beta);
    }

    #[test]
    fn reza_ameri_2_1(){
        let cardinality=  3u64;
        let hs = HyperGroup::new_from_matrix(
            &DMatrix::from_row_slice(
                cardinality as usize, 
                cardinality as usize, 
                &[1,6,6,6,1,1,6,1,1]));
        let beta = hs.beta_relation();
        assert!(beta.is_equivalence());
        let eq_classes=hs.collect_beta_classes();
        let expected_beta_classes = vec![(0u64, vec![0]), (1u64, vec![1, 2])];
        assert_eq!(eq_classes,expected_beta_classes)
    }
    #[test]
    fn reza_ameri_4(){
        let cardinality  =5u64;
        let input_array  = vec![
            vec![0],vec![1,2],vec![1,2],vec![3],vec![4],
            vec![1,2],vec![3],vec![3],vec![4],vec![0],
            vec![1,2],vec![3],vec![3],vec![4],vec![0],
            vec![3],vec![4],vec![4],vec![0],vec![1,2],
            vec![4],vec![0],vec![0],vec![1,2],vec![3]];

        let hg = HyperGroup::new_from_elements(&input_array, &cardinality);
        let beta = hg.beta_relation().quotient_set();
        let expected_beta = vec![(0_u64,vec![0]),(1_u64,vec![1,2]),(3_u64,vec![3]),(4_u64,vec![4])];
        assert_eq!(beta,expected_beta);
        let cardinality = 4u64;
        let fundamental_group = hg.get_isomorphic_fundamental_group();
        let input_array = vec![
            vec![0],vec![1],vec![2],vec![3],
            vec![1],vec![2],vec![3],vec![0],
            vec![2],vec![3],vec![0],vec![1],
            vec![3],vec![0],vec![1],vec![2]];
        let expected_fundamental_group = HyperGroup::new_from_elements(&input_array, &cardinality);
    assert_eq!(expected_fundamental_group,fundamental_group);
        let heart = hg.heart();
        let expected_heart = Some(HashSet::from([0]));
    assert_eq!(heart,expected_heart);
    }
    #[test]
    fn reza_armeri_6(){
        let cardinality  =7u64;
        let input_array  = vec![
            vec![0],vec![1],vec![2],vec![3],vec![4],vec![5,6], vec![5,6],
            vec![1],vec![0],vec![4],vec![5,6],vec![2],vec![3],vec![3],
            vec![2],vec![5,6],vec![0],vec![4],vec![3],vec![1],vec![1],
            vec![3],vec![4],vec![5,6],vec![0],vec![1],vec![2],vec![2],
            vec![4],vec![3],vec![1],vec![2],vec![5,6],vec![0],vec![0],
            vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
            vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
            ];

        let hg = HyperGroup::new_from_elements(&input_array, &cardinality);
        let beta = hg.beta_relation().quotient_set();
        let expected_beta = vec![(0_u64,vec![0]),(1_u64,vec![1]),(2_u64,vec![2]),(3_u64,vec![3]),(4_u64,vec![4]),(5_u64,vec![5,6])];
        assert_eq!(beta,expected_beta);
        let cardinality = 6u64;
        let fundamental_group = hg.get_isomorphic_fundamental_group();
        let input_array = vec![
            vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],
            vec![1],vec![0],vec![4],vec![5],vec![2],vec![3],
            vec![2],vec![5],vec![0],vec![4],vec![3],vec![1],
            vec![3],vec![4],vec![5],vec![0],vec![1],vec![2],
            vec![4],vec![3],vec![1],vec![2],vec![5],vec![0],
            vec![5],vec![2],vec![3],vec![1],vec![0],vec![4],
            ];
        let expected_fundamental_group = HyperGroup::new_from_elements(&input_array, &cardinality);
    assert_eq!(expected_fundamental_group,fundamental_group);
        let heart = hg.heart();
        let expected_heart = Some(HashSet::from([0]));
    assert_eq!(heart,expected_heart);
    }
    #[test]
    fn reza_armeri_7(){
let cardinality  =8u64;
        let input_array  = vec![
            vec![0],vec![1],vec![2],vec![3],vec![4],vec![5], vec![6],vec![7],
            vec![1,5],vec![0,4],vec![3,7],vec![2,6],vec![1,5],vec![0,4],vec![3,7],vec![2,6],
            vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![0],vec![1],
            vec![3,7],vec![2,6],vec![1,5],vec![0,4],vec![3,7],vec![2,6],vec![1,5],vec![0,4],
            vec![4],vec![5],vec![6],vec![7],vec![0],vec![1],vec![2],vec![3],
            vec![1,5],vec![0,4],vec![3,7],vec![2,6],vec![1,5],vec![0,4],vec![3,7],vec![2,6],
            vec![6],vec![7],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],
            vec![3,7],vec![2,6],vec![1,5],vec![0,4],vec![3,7],vec![2,6],vec![1,5],vec![0,4],
            ];

        let hg = HyperGroup::new_from_elements(&input_array, &cardinality);
        let beta = hg.beta_relation().quotient_set();
        let expected_beta = vec![(0_u64,vec![0,4]),(1_u64,vec![1,5]),(2_u64,vec![2,6]),(3_u64,vec![3,7])];
        assert_eq!(beta,expected_beta);
        let cardinality = 4u64;
        let fundamental_group = hg.get_isomorphic_fundamental_group();
        let input_array = vec![
            vec![0],vec![1],vec![2],vec![3],
            vec![1],vec![0],vec![3],vec![2],
            vec![2],vec![3],vec![0],vec![1],
            vec![3],vec![2],vec![1],vec![0]
            ];
        let expected_fundamental_group = HyperGroup::new_from_elements(&input_array, &cardinality);
    assert_eq!(expected_fundamental_group,fundamental_group);
        let heart = hg.heart();
        let expected_heart = Some(HashSet::from([0,4]));
    assert_eq!(heart,expected_heart);

    }
    #[test]
    fn reza_ameri_3_1(){
        let cardinality=  4u64;
        let hg = HyperGroup::new_from_matrix(
                &DMatrix::from_row_slice(
                    cardinality as usize, 
                    cardinality as usize, 
                    &[3,3,12,12,3,3,12,12,12,12,1,2,12,12,2,1]));
        let beta  =hg.beta_relation();
        assert!(beta.is_equivalence());
        let eq_classes=hg.collect_beta_classes();        
        let expected_beta_classes = vec![(0, vec![0, 1]), (2, vec![2, 3])];
        assert_eq!(eq_classes,expected_beta_classes)
    }
    #[test]
    fn reza_ameri_3_2(){
        let cardinality=  4u64;
        let hs = HyperGroup::new_from_matrix(
            &DMatrix::from_row_slice(
                cardinality as usize, 
                cardinality as usize, 
                &[1,6,6,8,6,8,8,1,6,8,8,1,8,1,1,6]));
        let beta = hs.beta_relation();
        assert!(beta.is_equivalence());
        let eq_classes=hs.collect_beta_classes();
        let expected_beta_classes = vec![(0, vec![0]), (1, vec![1, 2]), (3, vec![3])];
        assert_eq!(eq_classes,expected_beta_classes)
    }
    #[test]
    fn example_corsini_134() {
        let cardinality  =3u64;
        let hg=HyperGroup::new_from_matrix(
                &DMatrix::from_row_slice(
                    cardinality as usize, 
                    cardinality as usize, 
                    &[1,2,4,1,2,4,7,7,7]));
        let id = hg.collect_identities();
        
        assert!(id.is_none())
    }

    #[test]
    fn example_corsini_112_1() {
        let cardinality = 4u64;
        let matrix = 
            DMatrix::from_row_slice(cardinality as usize, cardinality as usize, 
                &[1,1,7,11,1,1,7,11,7,7,7,12,11,11,12,11]);
        let hg = HyperGroup::new_from_matrix(&matrix);
        let subset_a = 1u64;
        let subset_b = 7u64;
        let subset_c = 11u64;

        assert!(!hg.subhypergroup_is_closed(&subset_a));
        assert!(!hg.subhypergroup_is_closed(&subset_b));
        assert!(!hg.subhypergroup_is_closed(&subset_c));
    }
    #[test]
    fn example_corsini_112_3() {
        let cardinality = 4u64;
        let matrix = 
        DMatrix::from_row_slice(cardinality as usize, cardinality as usize, 
            &[1,4,4,14,4,2,4,13,4,4,4,15,14,13,15,15]);
        let hg = HyperGroup::new_from_matrix(&matrix);
        let subset_a = 1u64;
        let subset_b = 2u64;

        assert!(hg.subhypergroup_is_closed(&subset_a));
        assert!(hg.subhypergroup_is_closed(&subset_b));

        assert!(!hg.subhypergroup_is_invertible(&subset_a));
        assert!(!hg.subhypergroup_is_invertible(&subset_b)); 
    }
        #[test]
    fn example_corsini_112_4() {
        let cardinality = 5u64;
        let matrix = DMatrix::from_row_slice(cardinality as usize, cardinality as usize, &[1,2,4,8,24,2,3,24,28,28,4,24,5,26,26,8,28,26,31,31,24,28,26,31,31]);
        let hg = HyperGroup::new_from_matrix(&matrix);
        let subset_ab = 3u64;
        let subset_ac = 5u64;

        assert!(hg.subhypergroup_is_invertible(&subset_ab));
        assert!(hg.subhypergroup_is_invertible(&subset_ac));

        let subset_a=1u64;
        assert!(hg.subhypergroup_is_closed(&subset_a));
        assert!(!hg.subhypergroup_is_invertible(&subset_a));

    }
    #[test]
    ///
    /// Generate the transposition hypergroup with hyperoperation defined as 
    /// `x+y = {max(x,y)} if x != y`
    /// `x+y = {z ∈ H | z≤x} if x = y`.
    /// 
    fn example_transposition_hg() {
        let cardinality = 5;
        let function = {
            |a:usize,b:usize| 
                if a!=b {return 1<<a.max(b)}
                else {return (0..=a).into_iter().fold(0, |acc,x|acc|1<<x)}
            };
        let hg = HyperGroup::new_from_function(function, &cardinality).unwrap();
        assert!(hg.is_transposition());

    }
    #[test]
    ///
    /// Tests here come from `Some Remarks on Hyperstructures their Connections with Fuzzy Sets and Extensions to Weak Structures` by `Piergiulio Corsini`.
    /// These Hv-hypergroupoids of order 2, which are not associative.
    /// 
    fn example_fuzzy_grade() {
    //H12
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![0,1],
                                            vec![1],vec![0]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,1usize);
    //H13
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![1],
                                            vec![0,1],vec![0]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,1usize);
    //H15
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![0],
                                            vec![1],vec![0,1]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,1usize);
    //H9
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![1],
                                            vec![1],vec![0]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);
            
    //H10
    let cardinality =2u64;
    let input_values  = vec![vec![0],   vec![0,1],
                                            vec![1],vec![0]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);

    //H11
    let cardinality =2u64;
    let input_values  = vec![vec![1],   vec![0,1],
                                            vec![0],vec![1]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);
    //H14
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![0],
                                            vec![0],vec![0,1]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);
    //H16
    let cardinality =2u64;
    let input_values  = vec![vec![0],   vec![0,1],
                                            vec![0,1],vec![0]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);
   //H17
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![0,1],
                                            vec![0],vec![0,1]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);
   //H18
    let cardinality =2u64;
    let input_values  = vec![vec![0,1],   vec![0,1],
                                            vec![1],vec![0,1]];
    let hs = HyperGroupoid::new_from_elements(&input_values, &cardinality);
    let degree = hs.get_fuzzy_grade();
    assert_eq!(degree,2usize);

    }
}



