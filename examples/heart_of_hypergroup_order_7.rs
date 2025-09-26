use hyperstruc::hs::HyperGroupoid;

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
    let hg = match hypergroupoid.is_hypergroup() {
        true => hyperstruc::hypergroups::HyperGroup::new_from_hypergroupiod(&hypergroupoid),
        false => panic!()
    };
    let beta = hg.collect_beta_classes();
    println!("beta {:?}",beta);
    let fundamental = hg.get_isomorphic_fundamental_group();
    println!("f {}",fundamental);
    let heart = hg.heart();
    println!("heart is {:?}",heart);
}