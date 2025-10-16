use hyperstruc::{hs::HyperGroupoid, hypergroups::HyperGroup};

fn main() {
    let cardinality = 7u64;
    let input_array = vec![
        vec![0],vec![1],vec![2],vec![3],vec![4],vec![5,6],vec![5,6],
        vec![1],vec![0],vec![4],vec![5,6],vec![2],vec![3],vec![3],
        vec![2],vec![5,6],vec![0],vec![4],vec![3],vec![1],vec![1],
        vec![3],vec![4],vec![5,6],vec![0],vec![1],vec![2],vec![2],
        vec![4],vec![3],vec![1],vec![2],vec![5,6],vec![0],vec![0],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
    ];
    let hypergroupoid = HyperGroupoid::new_from_elements(&input_array, &cardinality);
    let hg = HyperGroup::new_from_hypergroupoid(hypergroupoid);
    let beta = hg.collect_beta_classes();
    println!("beta {:?}",beta);
    let fundamental = hg.get_fundamental_group();
    println!("f {}",fundamental);
    let heart = hg.heart();
    println!("heart is {:?}",heart);
    let isomorphic_fg = hg.get_isomorphic_fundamental_group();
    println!("{}",isomorphic_fg);
}