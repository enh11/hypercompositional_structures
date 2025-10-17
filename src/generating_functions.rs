pub fn b_hypercomposition() -> impl Fn(usize, usize) -> u64 {
    Box::new(move |a:usize,b:usize| {1<<a|1<<b})
}
pub fn tropical_hypergroup() -> impl Fn(usize, usize) -> u64 {
    Box::new(move |a: usize,b: usize| {
                    if a!=b {return 1<<a.max(b)}
                    else {return (0..=a).into_iter().fold(0, |acc,x|acc|1<<x)}
                
    })
}
pub fn genetics_hypergroup(cardinality: &u64) -> impl Fn(usize, usize) -> u64 {
    let h = (1u64 << cardinality) - 1;
    move |a: usize, b: usize| {
        let lower_elements = (0..a.min(b)).fold(0u64, |acc, x| acc | (1 << x));
        h - lower_elements
    }
}

