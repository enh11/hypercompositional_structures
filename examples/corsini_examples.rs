use hyperstruc::{hs::HyperGroupoid, hypergroups::HyperGroup};

fn main() {
// We define a hypergroupoid of order 3 and we look for its identities.
// First, we set the cardinality.
let cardinality = 7u64;
// Then, we define its cayley table.
let cayley_table_array = vec![
        vec![0],vec![1],vec![2],vec![3],vec![4],vec![5,6],vec![5,6],
        vec![1],vec![0],vec![4],vec![5,6],vec![2],vec![3],vec![3],
        vec![2],vec![5,6],vec![0],vec![4],vec![3],vec![1],vec![1],
        vec![3],vec![4],vec![5,6],vec![0],vec![1],vec![2],vec![2],
        vec![4],vec![3],vec![1],vec![2],vec![5,6],vec![0],vec![0],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
        vec![5,6],vec![2],vec![3],vec![1],vec![0],vec![4],vec![4],
    ]; 
// Now we can initialize the Hypergroupoid.
let hg = HyperGroupoid::new_from_elements(&cayley_table_array, &cardinality);
//Show the hypergruopoid
hg.show(); 
//Show the left_identities of hg 
hg.show_left_identities();
//Show the right identities of hg
hg.show_right_identities();
//Show the identities of hg
hg.show_identities();
//Show left scalar elements
hg.show_left_scalars();
//Show right scalar elements
hg.show_right_scalars();
// Show scalar elements
hg.show_scalars();
// Find left invertible elements in H
for x in 0..cardinality {
    hg.show_left_inverses_of_x(&x);
}
// Find right invertible elements in H
for x in 0..cardinality {
    hg.show_right_inverses_of_x(&x);
}
// We can check the hypergroupoid is a hypergroup
println!("Is hypergroup: {}", hg.is_hypergroup());
// Now we build the hypergroup from the hypergroupoid, in order to analyze subshypergroups
let hg = HyperGroup::new_from_hypergroupiod(&hg);
// We collect all non trivial subhypergroup
let sub_hg = hg.collect_proper_subhypergroups();
print!("{:?}",sub_hg);





}