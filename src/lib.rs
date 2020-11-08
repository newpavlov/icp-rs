//use num_traits::identities::Zero;
use nalgebra::base::{Vector3, Unit, Matrix3, MatrixMN, DVector};
use nalgebra::geometry::{Rotation3, Translation3, Point3};
use nalgebra::base::dimension::{U6, Dynamic};
use std::iter::ExactSizeIterator;

mod voxel_bucket;

type Matrix6N = MatrixMN<f32, Dynamic, U6>;
pub type Normal = Unit<Vector3<f32>>;
pub type Rotation = Rotation3<f32>;
pub type Translation = Translation3<f32>;
pub type Point = Point3<f32>;

const NORMALS_THRESH: usize = 5;

pub struct Icp {
    vbt: voxel_bucket::VoxelBucket,
    max_iter: u32,
    dist_delta: f32,
    normals: Vec<Normal>,
}

pub fn calc_normal(points: impl Iterator<Item=Point> + ExactSizeIterator + Clone) -> Normal {
    let n = points.len() as f32;
    let cm: Vector3<f32> = points
        .clone()
        .fold(Vector3::zeros(), |a, p| a + p.coords)/n;

    let cov = points
        .map(|p| p.coords - cm)
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

    Normal::new_normalize(Vector3::new(u[(0, i)], u[(1, i)], u[(2, i)]))
}

impl Icp {
    pub fn new(
        ref_scan: &[Point], search_radius: f32,
        max_iter: u32, dist_delta: f32,
    ) -> Result<Self, voxel_bucket::ConversionError> {
        let vbt = voxel_bucket::VoxelBucket::new(ref_scan, search_radius)?;
        let dn = Normal::new_normalize(Vector3::new(1.0, 0.0, 0.0));
        let normals = vec![dn; ref_scan.len()];
        let mut s = Self { vbt, max_iter, dist_delta, normals };
        s.calc_normals();
        return Ok(s)
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
                self.normals[idx as usize] = calc_normal(points_buf.iter().cloned());
            } else {
                remove_idxs.push(idx);
            }

            points_buf.clear();
        }
        // remove points for which normal vector was not calculated
        self.vbt.remove_points(&remove_idxs);
    }

    pub fn register(
        &self, scan: &[Point], mut r: Rotation, mut t: Translation,
    ) -> (Rotation, Translation, u32, f32) {
        let mut corresp = 0u32;
        let mut sum_dist = 0.0;

        let mut a = Vec::new();
        let mut b = Vec::new();
        for _ in 0..self.max_iter {
            sum_dist = 0.0;
            corresp = 0;
            a.clear();
            b.clear();
            for &s in scan.iter() {
                let s = t*r*s;
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
            sum_dist /= a.len() as f32;


            // TODO: add error
            if corresp < 100 {
                return (r, t, 0, 0.0);
            }

            let a = Matrix6N::from_row_slice(&a);
            let b = DVector::from_row_slice(&b);

            let res = nalgebra::linalg::SVD::new(a, true, true)
                .solve(&b, 1e-10)
                .unwrap();

            let new_t = Translation::new(res[3], res[4], res[5]);
            let new_r = Rotation::from_euler_angles(res[0], res[1], res[2]);

            t = Translation::from(new_r*t.vector + new_t.vector);
            r = new_r*r;
            let delta = new_t.vector.norm();

            if delta < self.dist_delta {
                break;
            }
        }
        (r, t, corresp, sum_dist)
    }
}
