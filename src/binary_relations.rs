pub mod relations;
pub mod relation_matrix;

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use crate::binary_relations::relations::Relation;
    #[test]
    fn exercize_9_4_1() {
        let set_a:HashSet<u64>  = (0..=3u64).into_iter().collect();
        let set_b:HashSet<u64>  = (0..=3u64).into_iter().collect();
        let relation_vec = vec![
            (0,1),(1,1),(1,2),(2,0),(2,2), (3,0)
        ];
        let relation = Relation{
            a: set_a,
            b: set_b,
            rel: relation_vec,
        };
        let ref_cl = relation.reflexive_closure();
        let sym_cl = relation.symmetric_closure();

        let expected_reflexive = vec![ (0,0),(0,1),(1,1),(1,2),(2,0),(2,2),(3,0),(3,3)];
        let expected_sym = vec![ (0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0)];
        assert_eq!(expected_reflexive,ref_cl.rel);
        assert_eq!(expected_sym,sym_cl.rel);
    }
}