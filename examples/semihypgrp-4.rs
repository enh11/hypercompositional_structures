use hyperstruc::hs::hypergroupoids::HyperGroupoid;
use num_traits::PrimInt;


fn main() {

let cardinality = 5u64;
let input_array1 = vec![
        vec![0],vec![0,1,3],vec![0],vec![0,1,3],vec![0,1,3],
        vec![0],vec![1],vec![0],vec![0,1,3], vec![0,1,3],
        vec![0], vec![0,1,3],vec![0,2],vec![0,1,3],vec![0,1,2,3,4],
        vec![0],vec![0,1,3],vec![0],vec![0,1,3],vec![0,1,3],
        vec![0],vec![0,1,3],vec![0,2],vec![0,1,3],vec![0,1,2,3,4],
];
let hs = HyperGroupoid::new_from_elements(&input_array1, &cardinality);
    println!("{}",hs);
    hs.show_left_identities();
    hs.show_right_identities();
    hs.show_identities();

    hs.show_right_scalars();
    hs.show_left_scalars();
    hs.show_scalars();

    hs.show_partial_left_identities();

    let rs = hs.collect_right_scalars();
    println!("{:?}",rs);
for x in 0..cardinality{
    hs.show_left_inverses_of_x(&x);
    hs.show_right_inverses_of_x(&x);
    hs.show_inverses_of_x(&x);
    }
}