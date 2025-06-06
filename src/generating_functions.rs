pub fn b_hypercomposition() -> impl Fn(u64, u64) -> u64 {
    Box::new(move |a:u64,b:u64| {1<<a|1<<b})
}
pub fn tropical_hypergroup() -> impl Fn(u64, u64) -> u64 {
    Box::new(move |a: u64,b: u64| {
                    if a!=b {return 1<<a.max(b)}
                    else {return (0..=a).into_iter().fold(0, |acc,x|acc|1<<x)}
                
    })
}
pub fn genetics_hypergroup(cardinality: &u64) -> impl Fn(u64, u64) -> u64 {
    let h = (1u64 << cardinality) - 1;
    move |a: u64, b: u64| {
        let lower_elements = (0..a.min(b)).fold(0u64, |acc, x| acc | (1 << x));
        h - lower_elements
    }
}

