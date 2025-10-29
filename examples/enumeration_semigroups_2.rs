use hyperstruc::enumeration::enumeration_hyperstructure;


fn main() {

let cardinality =2u64;
// This produce a .txt file containing alla semigroups of order 2 up to isomorphism.
let t = enumeration_hyperstructure("semigroups", &cardinality);
println!("t = {:?}",t);

}