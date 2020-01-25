//use num_traits::identities::Zero;
use rayon::prelude::*;
use log::debug;
use nalgebra::base::Matrix3;

mod voxel_bucket;

type Point = nalgebra::base::Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = nalgebra::base::Vector3<f32>;

type KdTree = kdtree::KdTree<f32, usize, [f32; 3]>;


pub struct Icp<'a> {
    kdt: KdTree,
    ref_scan: &'a [Point],
    dist_thresh: f32,
    max_iter: u32,
    dist_delta: f32,
}

fn mass_centeres(
    ref_scan: &[Point], scan: &[Point], corresp: &[[u32; 2]],
) -> [Point; 2] {
    let mut n = 0u32;
    let mut ref_accum = Point::zeros();
    let mut accum = Point::zeros();
    for &[i, j] in corresp {
        ref_accum += ref_scan[i as usize];
        accum += scan[j as usize];
        n += 1;
    }
    let n = n as f32;
    [ref_accum/n, accum/n]
}

fn calc_w_mc(
    ref_scan: &[Point], scan: &[Point], corresp: &[[u32; 2]],
) -> (Matrix3<f32>, Point, Point) {
    let mut w = Matrix3::zeros();
    let [ref_mc, mc] = mass_centeres(ref_scan, scan, corresp);
    for &[i, j] in corresp {
        let ref_p = ref_scan[i as usize] - ref_mc;
        let p = scan[j as usize] - mc;
        w += ref_p*p.transpose();
    }
    (w, ref_mc, mc)
}

impl<'a> Icp<'a> {
    pub fn new(
        ref_scan: &'a [Point], dist_thresh: f32, max_iter: u32, dist_delta: f32,
    ) -> Self {
        let n = ref_scan.len();
        let mut kdt = KdTree::with_capacity(3, n);
        debug!("Constructing KD-tree.");
        for (i, p) in ref_scan.iter().enumerate() {
            kdt.add([p[0], p[1], p[2]], i)
                .expect("failed to add a point to KD-tree");            
        }
        debug!("KD-tree completed.");
        Self { kdt, ref_scan, dist_thresh: dist_thresh*dist_thresh, max_iter, dist_delta }
    }

    pub fn register(
        &self, scan: &[Point], mut r: Rot, mut t: Trans,
    ) -> (Rot, Trans, u32) {
        let mut corresp_u32 = 0;
        let mut prev_t = t;
        for n in 0..self.max_iter {
            println!("\n\n===========================\nICP iteration: {}", n);

            let (corresp, sum_dist) = scan.par_iter()
                .enumerate()
                .fold(|| (Vec::new(), 0.0f32), |mut a, (j, p)| {
                    use kdtree::distance::squared_euclidean;
                    let p = r*p + t;
                    let v = self.kdt
                        .nearest(&[p[0], p[1], p[2]], 1, &squared_euclidean)
                        .expect("kdtree::nearest failure");
                    for (d, &i) in v {
                        if d > self.dist_thresh {
                            continue;
                        }

                        a.0.push([i as u32, j as u32]);
                        a.1 += d;
                    }
                    a
                })
                .reduce(|| (Vec::new(), 0.0f32), |mut a, v| {
                    a.0.extend_from_slice(&v.0);
                    a.1 += v.1;
                    a
                });

            let (w, ref_mc, mc) = calc_w_mc(self.ref_scan, scan, &corresp);

            let svd_res = w.svd(true, true);
            let u = svd_res.u.unwrap();
            let v_t = svd_res.v_t.unwrap();

            r = u*v_t;
            //println!("new_r: {}", new_r);
            let new_r_det = r.determinant();
            println!("new_r_det: {}", new_r_det);
            if new_r_det < 0.0 {
                let mut m = Rot::identity();
                m.m33 = -1.;
                r = u*m*v_t;
            }
            t = ref_mc - r*mc;
            
            corresp_u32 = corresp.len() as u32;
            let corresp_f32 = corresp.len() as f32;
            let delta = (t - prev_t).norm();
            prev_t = t;
            println!(
                "start error: {}\ncorresp: {}\nt_delta: {}\nr: {}t: {}",
                sum_dist/corresp_f32, corresp.len(), delta, r, t
            );
            if delta < self.dist_delta {
                break;
            }
        }
        (r, t, corresp_u32)
    }
}
