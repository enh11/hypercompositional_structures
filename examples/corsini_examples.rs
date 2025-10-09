use hyperstruc::{hs::HyperGroupoid, hypergroups::HyperGroup};

fn main() {
// We define a hypergroupoid of order 3 and we look for its identities.
// First, we set the cardinality.
 let cardinality = 3u64; 
 // Then, we define its cayley table.
 let cayley_table_array = vec![
    vec![0],     vec![1],     vec![2],
    vec![0],     vec![1],     vec![2],
    vec![0,1,2], vec![0,1,2], vec![0,1,2]
    ];
// Now we can initialize the Hypergroupoid.
let hg = HyperGroupoid::new_from_elements(&cayley_table_array, &cardinality);
//Show the left_identities of hg 
hg.show_left_identities();
//Show the right identities of hg
hg.show_right_identities();
//Show the identities of hg
hg.show_identities();
// We can check if the hypergroupoid is a hypergroup
println!("Is hypergroup: {}", hg.is_hypergroup());





}