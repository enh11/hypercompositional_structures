use hyperstruc::hs::HyperGroupoid;

fn main() {
let cardinality = 10u64;
let input_array1 = vec![
        vec![0],vec![1],vec![2],vec![3],vec![4],vec![0,5],vec![6],vec![7],vec![8],vec![9],
        vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],
        vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],
        vec![3],vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],
        vec![4],vec![5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],
        vec![0,5],vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],
        vec![6],vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],
        vec![7],vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],
        vec![8],vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],
        vec![9],vec![0],vec![1],vec![2],vec![3],vec![4],vec![5],vec![6],vec![7],vec![8]
    ];
let hs1 = HyperGroupoid::new_from_elements(&input_array1, &cardinality);
    println!("{}",hs1);
    
for x in 0..cardinality{
    println!("inverses of {:?} are {:?}",x,hs1.show_inverses_of_x(&x));
    }
}