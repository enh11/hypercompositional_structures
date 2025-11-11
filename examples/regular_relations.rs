use hyperstruc::{binary_relations::relations::{self, Relation}, hg_3::tag_hypergroups_3::TAGS_HG_3, hs::hypergroups::HyperGroup};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};


fn main() {
let cardinality =3u64;

TAGS_HG_3.par_iter().all(|tag| {
    let hg = HyperGroup::new_from_tag_u128(tag,&cardinality);
    hg.beta_relation().is_left_regular(&hg.0)

}
);
println!("true");
TAGS_HG_3.par_iter().all(|tag| {
    let hg = HyperGroup::new_from_tag_u128(tag,&cardinality);
    hg.beta_relation().is_right_regular(&hg.0)

}
);
println!("true");
}