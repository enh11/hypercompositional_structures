use hyperstruc::hs::hypergroupoids::HyperGroupoid;
use rayon::prelude::*;

fn main() {
    for card in 3..=12 {
        let ms: Vec<u64> = (0..500).into_par_iter()
            .map(|_| {
                let hg = HyperGroupoid::new_random_from_cardinality(&card);
                hg.collect_all_finite_hyperproducts().1
            })
            .collect();

        let count = ms.len() as f64;
        let sum: u64 = ms.iter().sum();
        let mean = sum as f64 / count;

        let variance = ms.iter()
            .map(|&x| {
                let dx = x as f64 - mean;
                dx * dx
            })
            .sum::<f64>() / count;
        let std_dev = variance.sqrt();

        let min = ms.iter().min().unwrap();
        let max = ms.iter().max().unwrap();
        let ratio = mean / card as f64;

        println!(
            "card {:2} mean {:6.2} std {:6.2} min {:2} max {:2} ratio {:.3}",
            card, mean, std_dev, min, max, ratio
        );
    }
}
