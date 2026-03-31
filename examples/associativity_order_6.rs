use hyperstruc::hs::hypergroupoids::HyperGroupoid;

fn main() {
let cardinality = 6u64;
let input_array = vec![
vec![0,1],vec![0,1],vec![2,3],vec![2,3],vec![4,5],vec![4,5],
vec![0,1],vec![0,1],vec![2,3],vec![2,3],vec![4,5],vec![4,5],
vec![2,3],vec![2,3],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],
vec![3],vec![3],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],
vec![4,5],vec![4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],
vec![5],vec![5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],vec![0,1,2,3,4,5],
];
let h = HyperGroupoid::new_from_elements(&input_array, &cardinality);
h.show();
println!("Associative: {}", h.is_associative());
}