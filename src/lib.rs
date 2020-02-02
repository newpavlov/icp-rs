//use num_traits::identities::Zero;
use nalgebra::base::{Matrix3, MatrixMN, DVector};
use nalgebra::base::dimension::{U6, Dynamic};

mod voxel_bucket;

type Point = nalgebra::base::Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = nalgebra::base::Vector3<f32>;
type Normal = nalgebra::base::Vector3<f32>;

type Matrix6N = MatrixMN<f32, Dynamic, U6>;

const NORMALS_THRESH: usize = 10;

pub struct Icp {
    vbt: voxel_bucket::VoxelBucket,
    max_iter: u32,
    dist_delta: f32,
    normals: Vec<Point>,
}

fn calc_normal(points: &[Point]) -> Normal {
    let n = points.len() as f32;
    let cm: Point = points.iter()
        .fold(Point::zeros(), |a, p| a + p)/n;

    let cov = points.iter()
        .map(|p| p - cm)
        .map(|p| p*p.transpose())
        .fold(Matrix3::zeros(), |a, v| a + v);

    let svd_res = cov.svd(true, false);

    let u = svd_res.u.unwrap();
    let sing = svd_res.singular_values;
    let i = if sing[0] < sing[1] && sing[0] < sing[2] {
        0
    } else if sing[1] < sing[0] && sing[1] < sing[2] {
        1
    } else {
        2
    };

    Normal::new(u[(0, i)], u[(1, i)], u[(2, i)])
}

impl Icp {
    pub fn new(
        ref_scan: &[Point], search_radius: f32,
        max_iter: u32, dist_delta: f32,
    ) -> Result<Self, voxel_bucket::ConversionError> {
        let vbt = voxel_bucket::VoxelBucket::new(ref_scan, search_radius)?;
        let normals = vec![Normal::zeros(); ref_scan.len()];
        Ok(Self { vbt, max_iter, dist_delta, normals })
    }

    pub fn new_with_normals(
        ref_scan: &[Point], normals: &[Normal],
        search_radius: f32, max_iter: u32, dist_delta: f32,
    ) -> Result<Self, voxel_bucket::ConversionError> {
        assert_eq!(ref_scan.len(), normals.len());
        let vbt = voxel_bucket::VoxelBucket::new(ref_scan, search_radius)?;
        Ok(Self { vbt, max_iter, dist_delta, normals: normals.to_vec() })
    }

    pub fn get_points(&self) -> Vec<Point> {
        self.vbt.iter_points().map(|&(p, _)| p).collect()
    }

    pub fn get_points_normals(&self) -> Vec<(Point, Normal)> {
        self.vbt.iter_points()
            .map(|&(p, idx)| (p, self.normals[idx as usize]))
            .collect()
    }

    pub fn calc_normals(&mut self) {
        let mut remove_idxs = Vec::new();
        let mut points_buf = Vec::new();
        for &(p, idx) in self.vbt.iter_points() {
            self.vbt
                .inside_radius(p, |p, _, _| points_buf.push(p))
                .unwrap();
            if points_buf.len() >= NORMALS_THRESH {
                self.normals[idx as usize] = calc_normal(&points_buf);
            } else {
                remove_idxs.push(idx);
            }

            points_buf.clear();
        }
        // remove points for which normal vector was not calculated
        //println!("removing: {:?}", remove_idxs.len());
        self.vbt.remove_points(&remove_idxs);
    }

    pub fn register(
        &self, scan: &[Point], mut r: Rot, mut t: Trans,
    ) -> (Rot, Trans, u32, f32) {
        let mut corresp = 0u32;
        let mut sum_dist = 0.0;

        let mut a = Vec::new();
        let mut b = Vec::new();
        for _ in 0..self.max_iter {
            /*
            use rayon::prelude::*;

            #[derive(Default)]
            struct IterData {
                a: Vec<f32>, b: Vec<f32>,
                sum_dist: f32, corresp: u32,
            }

            let res = scan.par_iter()
                .fold(|| IterData::default(), |mut accum, s| {
                    let s = r*s + t;
                    let (d, idx, dist) = match self.vbt.search_closest(s) {
                        Some(v) => v,
                        None => return accum,
                    };

                    let n = self.normals[idx as usize];
                    accum.sum_dist += dist;
                    accum.corresp += 1;

                    accum.a.extend_from_slice(&[
                        n.z*s.y - n.y*s.z,
                        n.x*s.z - n.z*s.x,
                        n.y*s.x - n.x*s.y,
                        n.x,
                        n.y,
                        n.z,
                    ]);
                    let p1 = n.x*d.x + n.y*d.y + n.z*d.z;
                    let p2 = n.x*s.x + n.y*s.y + n.z*s.z;
                    accum.b.push(p1 - p2);
                    accum
                })
                .reduce(|| IterData::default(), |mut accum, val| {
                    accum.sum_dist += val.sum_dist;
                    accum.corresp += val.corresp;
                    accum.a.extend_from_slice(&val.a);
                    accum.b.extend_from_slice(&val.b);
                    accum
                });

            sum_dist = res.sum_dist;
            corresp = res.corresp;
            */

            sum_dist = 0.0;
            corresp = 0;
            a.clear();
            b.clear();
            for &s in scan.iter() {
                let s = r*s + t;
                let (d, idx, dist) = match self.vbt.search_closest(s) {
                    Some(v) => v,
                    None => continue,
                };

                let n = self.normals[idx as usize];
                sum_dist += dist;
                corresp += 1;

                a.extend_from_slice(&[
                    n.z*s.y - n.y*s.z,
                    n.x*s.z - n.z*s.x,
                    n.y*s.x - n.x*s.y,
                    n.x,
                    n.y,
                    n.z,
                ]);
                let p1 = n.x*d.x + n.y*d.y + n.z*d.z;
                let p2 = n.x*s.x + n.y*s.y + n.z*s.z;
                b.push(p1 - p2);
            }

            let a = Matrix6N::from_row_slice(&a);
            let b = DVector::from_row_slice(&b);

            let x = nalgebra::linalg::SVD::new(a, true, true)
                .solve(&b, 1e-10)
                .unwrap();

            let new_t = Trans::new(x[3], x[4], x[5]);
            let (sin_a, cos_a) = x[0].sin_cos();
            let (sin_b, cos_b) = x[1].sin_cos();
            let (sin_y, cos_y) = x[2].sin_cos();
            let new_r = Rot::new(
                cos_y*cos_b,
                -sin_y*cos_a + cos_y*sin_b*sin_a,
                sin_y*sin_a + cos_y*sin_b*cos_a,
                sin_y*cos_b,
                cos_y*cos_a + sin_y*sin_b*sin_a,
                -cos_y*sin_a + sin_y*sin_b*cos_a,
                -sin_b,
                cos_b*sin_a,
                cos_b*cos_a,
            );

            let prev_t = t;
            t += r*new_t;
            r = r*new_r;
            let delta = (t - prev_t).norm();

            //println!("{:?}", sum_dist/(corresp as f32));
            /*println!(
                "start error: {}\ncorresp: {}\nt_delta: {}\nr: {}t: {}",
                sum_dist/(corresp as f32), corresp, delta, r, t
            );*/
            if delta < self.dist_delta {
                break;
            }
        }
        (r, t, corresp, sum_dist)
    }
}
