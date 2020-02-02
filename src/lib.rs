//use num_traits::identities::Zero;
use nalgebra::base::Matrix3;

mod voxel_bucket;

type Point = nalgebra::base::Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = nalgebra::base::Vector3<f32>;
type Normal = nalgebra::base::Vector3<f32>;

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
        println!("removing: {:?}", remove_idxs.len());
        self.vbt.remove_points(&remove_idxs);
    }

    pub fn register(
        &self, scan: &[Point], mut r: Rot, mut t: Trans,
    ) -> (Rot, Trans, u32, f32) {
        let mut corresp = 0u32;
        let mut sum_dist = 0.0;
        let mut prev_t = t;
        let mut buf = vec![None; scan.len()].into_boxed_slice();

        use rayon::prelude::*;

        for _ in 0..self.max_iter {
            //println!("\n\n===========================\nICP iteration: {}", n);

            /*let mut ref_accum = Point::zeros();
            let mut accum = Point::zeros();
            sum_dist = 0.0;
            corresp = 0;
            for (&p, buf_ref) in scan.iter().zip(buf.iter_mut()) {
                let p2 = r*p + t;
                match self.vbt.search_closest(&p2) {
                    Some((p_ref, dist)) => {
                        *buf_ref = Some(p_ref);
                        sum_dist += dist;
                        corresp += 1;
                        ref_accum += *p_ref;
                        accum += p;
                    }
                    None => *buf_ref = None,
                }
            }
            */

            struct IterData {
                ref_accum: Point, accum: Point, sum_dist: f32, corresp: u32,
            }

            impl Default for IterData {
                fn default() -> Self {
                    let z = Point::zeros();
                    Self { ref_accum: z, accum: z, sum_dist: 0.0, corresp: 0 }
                }
            }

            let itd: IterData = scan.par_iter().zip(buf.par_iter_mut())
                .map(|(&p, buf_ref)| {
                    let p2 = r*p + t;
                    match self.vbt.search_closest(p2) {
                        Some((p_ref, idx, dist)) => {
                            *buf_ref = Some(p_ref);
                            IterData {
                                sum_dist: dist,
                                corresp: 1,
                                ref_accum: p_ref,
                                accum: p,
                            }
                        }
                        None => {
                            *buf_ref = None;
                            IterData::default()
                        }
                    }
                })
                .reduce(|| IterData::default(),  |mut a, v| {
                    a.sum_dist += v.sum_dist;
                    a.corresp += v.corresp;
                    a.ref_accum += v.ref_accum;
                    a.accum += v.accum;
                    a
                });
            sum_dist = itd.sum_dist;
            corresp = itd.corresp;

            let corresp_f32 = corresp as f32;
            let ref_cm = itd.ref_accum/corresp_f32;
            let cm = itd.accum/corresp_f32;

            // calculate W
            let mut w = Matrix3::zeros();
            for (&ref_p, p) in buf.iter().zip(scan.iter()) {
                let ref_p = match ref_p {
                    Some(v) => v,
                    None => continue,
                };
                let ref_p = ref_p - ref_cm;
                let p = p - cm;
                //ref_p[2] = 0.0;
                //p[2] = 0.0;
                w += ref_p*p.transpose();
            }

            let svd_res = w.svd(true, true);
            let u = svd_res.u.unwrap();
            let v_t = svd_res.v_t.unwrap();

            r = u*v_t;
            //println!("new_r: {}", new_r);
            let new_r_det = r.determinant();
            //println!("new_r_det: {}", new_r_det);
            if new_r_det < 0.0 {
                let mut m = Rot::identity();
                m.m33 = -1.;
                r = u*m*v_t;
            }
            t = ref_cm - r*cm;

            let delta = (t - prev_t).norm();
            prev_t = t;
            /*println!(
                "start error: {}\ncorresp: {}\nt_delta: {}\nr: {}t: {}",
                sum_dist/corresp_f32, corresp, delta, r, t
            );
            */
            if delta < self.dist_delta {
                break;
            }
        }
        (r, t, corresp, sum_dist)
    }
}
