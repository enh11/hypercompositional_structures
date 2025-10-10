//! Example: Analyzing a Hypergroupoid and Its Substructures
//!
//! This example demonstrates how to use the `hyperstruc` crate to:
//! 1. Define a **hypergroupoid** from a Cayley table.
//! 2. Explore its **algebraic properties** (identities, scalars, inverses).
//! 3. Check whether it forms a **hypergroup**.
//! 4. Analyze its **subhypergroups** if applicable.
//!
//! Run this example with:
//! ```bash
//! cargo run --example hypergroupoid_analysis
//! ```

use hyperstruc::{hs::HyperGroupoid, hypergroups::HyperGroup};

fn main() {
    // === 1. Define the hypergroupoid ===
    // We'll create a hypergroupoid with 7 elements (0..6).
    let cardinality = 7u64;

    // Cayley table for the hyperoperation (*)
    // Each entry represents the possible results of a * b as a set.
    // Example: cayley_table_array[a][b] = {result(s) of a * b}
    let cayley_table_array = vec![
        vec![0], vec![1], vec![2], vec![3], vec![4], vec![5, 6], vec![5, 6],
        vec![1], vec![0], vec![4], vec![5, 6], vec![2], vec![3], vec![3],
        vec![2], vec![5, 6], vec![0], vec![4], vec![3], vec![1], vec![1],
        vec![3], vec![4], vec![5, 6], vec![0], vec![1], vec![2], vec![2],
        vec![4], vec![3], vec![1], vec![2], vec![5, 6], vec![0], vec![0],
        vec![5, 6], vec![2], vec![3], vec![1], vec![0], vec![4], vec![4],
        vec![5, 6], vec![2], vec![3], vec![1], vec![0], vec![4], vec![4],
    ];

    // Initialize the hypergroupoid from the Cayley table
    let hg = HyperGroupoid::new_from_elements(&cayley_table_array, &cardinality);

    println!("=== Hypergroupoid Analysis ===");
    hg.show();

    // === 2. Inspect algebraic properties ===
    println!("\n=== Identities ===");
    hg.show_left_identities();
    hg.show_right_identities();
    hg.show_identities();

    println!("\n=== Scalars ===");
    hg.show_left_scalars();
    hg.show_right_scalars();
    hg.show_scalars();

    println!("\n=== Inverses ===");
    for x in 0..cardinality {
        println!("Element {x}:");
        hg.show_left_inverses_of_x(&x);
        hg.show_right_inverses_of_x(&x);
    }

    // === 3. Check if it forms a hypergroup ===
    let is_hypergroup = hg.is_hypergroup();
    println!("\nIs hypergroup: {is_hypergroup}");

    // === 4. Analyze subhypergroups (if hypergroup) ===
    if is_hypergroup {
        let hypergroup = HyperGroup::new_from_hypergroupiod(&hg);

        println!("\n=== Subhypergroup Analysis ===");
        hypergroup.show_proper_subhypergroups();
        hypergroup.show_proper_invertible_subhypergroups();
        hypergroup.show_proper_closed_subhypergroups();
        hypergroup.show_proper_reflexive_subhypergroups();
        hypergroup.show_proper_normal_subhypergroups();
    } else {
        println!("Not a hypergroup, skipping subhypergroup analysis.");
    }
}
