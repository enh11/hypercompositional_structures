use std::fmt;

pub mod hypergroupoids;
pub mod hypergroups;
pub mod quotient_hg;
pub mod fuzzy;
pub mod ordered_semigroup;

#[derive(Debug, Clone)]
pub enum HyperStructureError {
    NotHypergroup,
    NotPreOrderedSemigroup,
    NotPreOrder,
    NotAssociative,
    NotSemigroup,
    NotRegularRelation,
    ListOfSemigroupsNotAvailable
}
impl fmt::Display for HyperStructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HyperStructureError::NotHypergroup => write!(f, "The structure is not a valid hypergroup."),
            HyperStructureError::NotPreOrderedSemigroup => write!(f, "The preorder is not compatible with the semigroup operation."),
            HyperStructureError::NotPreOrder => write!(f,"The relation is not a preorder."),
            HyperStructureError::NotAssociative => write!(f,"The operation is not associative"),
            HyperStructureError::NotSemigroup => write!(f,"The structure is not a semigroup."),
            HyperStructureError::NotRegularRelation => write!(f,"The relation is not regular."),
            HyperStructureError::ListOfSemigroupsNotAvailable => write!(f, "The lists of semigroups are available up to order 5.")
                    }
    }
}