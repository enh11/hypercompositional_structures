use std::time::Instant;
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
    let now = Instant::now();
    let heart1 = hg.heart_fast();
    let end1 = now.elapsed();
        println!("time {:?}",now.elapsed());
    let now = Instant::now();
    let heart2 = hg.heart();
    let end2= now.elapsed();
        println!("time {:?}",now.elapsed());
    assert_eq!(heart1,heart2.unwrap());
    println!("ratio {}",end2.as_secs_f64()/end1.as_secs_f64());
}