use nalgebra::base::{MatrixMN, DVector};
use nalgebra::base::dimension::{U6, Dynamic};

// column major
type Matrix6N = MatrixMN<f32, U6, Dynamic>;
type MatrixN6 = MatrixMN<f32, Dynamic, U6>;

fn main() {
    let b = DVector::new_random(50_000);
    let m = MatrixN6::new_random_generic(Dynamic::new(50_000), U6);
    let t = std::time::Instant::now();
    let x = m.svd(true, true).solve(&b, 0.0).unwrap();
    let dt = t.elapsed();
    println!("{:?} {}", dt, x);
}