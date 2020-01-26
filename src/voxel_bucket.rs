use std::fmt;
use fxhash::FxHashMap;
use super::Point;

pub struct VoxelBucket {
    map: FxHashMap<[i16; 3], Vec<Point>>,
    radius: f32,
}

#[derive(Debug)]
pub struct ConversionError;

impl fmt::Display for ConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("failed to convert f32 into the desired range")
    }
}

impl std::error::Error for ConversionError {}

impl VoxelBucket {
    pub fn new(points: &[Point], radius: f32) -> Result<Self, ConversionError> {
        let mut s = Self { map: Default::default(), radius };
        for point in points {
            let key = convert2key(point, radius)?;
            s.map.entry(key).or_default().push(*point);
        }
        Ok(s)
    }

    /// Returns closest point inside the radius if any
    /// and square distance to it
    pub fn search_closest(&self, p: &Point) -> Option<(&Point, f32)> {
        let r = self.radius;
        let key = match convert2key(p, r) {
            Ok(k) => k,
            Err(_) => return None,
        };
        let mut min_dist2 = r*r;
        let mut closest = None;

        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let new_key = [key[0] + i, key[1] + j, key[2] + k];
                    let points = match self.map.get(&new_key) {
                        Some(v) => v,
                        None => continue,
                    };
                    for p2 in points.iter() {
                        let dx = p2.x - p.x;
                        let dy = p2.y - p.y;
                        let dz = p2.z - p.z;
                        let d2 = dx*dx + dy*dy + dz*dz;
                        if d2 < min_dist2 {
                            min_dist2 = d2;
                            closest = Some(p2);
                        }
                    }
                }
            }
        }
        closest.map(|p| (p, min_dist2))
    }
}

fn convert(v: f32, r: f32) -> Result<i16, ConversionError> {
    let v = (v/r).round();
    if std::i16::MIN as f32 <= v && v <= std::i16::MAX as f32 {
        Ok(v as i16)
    } else {
        Err(ConversionError)
    }
}

fn convert2key(p: &Point, r: f32) -> Result<[i16; 3], ConversionError> {
    Ok([convert(p.x, r)?, convert(p.y, r)?, convert(p.z, r)?])
}

#[cfg(test)]
mod tests {
    use super::{Point, VoxelBucket};

    #[test]
    fn random_points() {
        use rand::{SeedableRng, Rng};
        const RADIUS: f32 = 1.0;

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let mut points = Vec::new();
        // get points inside two boxes
        for _ in 0..500 {
            points.push(Point::new(
                rng.gen_range(-10.0, 0.0),
                rng.gen_range(-3.0, 5.0),
                rng.gen_range(1.0, 3.0),
            ));
            points.push(Point::new(
                rng.gen_range(2.0, 10.0),
                rng.gen_range(-10.0, -5.0),
                rng.gen_range(3.0, 8.0),
            ));
        }

        let vbt = VoxelBucket::new(&points, RADIUS).unwrap();
        // test 5k points
        for _ in 0..5000 {
            let p = Point::new(
                rng.gen_range(-10.0, 10.0),
                rng.gen_range(-10.0, 10.0),
                rng.gen_range(-10.0, 10.0),
            );
            let mut closest = None;
            let mut min_dist2 = RADIUS*RADIUS;
            for p2 in points.iter() {
                let dx = p2.x - p.x;
                let dy = p2.y - p.y;
                let dz = p2.z - p.z;
                let d2 = dx*dx + dy*dy + dz*dz;
                if d2 < min_dist2 {
                    min_dist2 = d2;
                    closest = Some((p2, min_dist2));
                }
            }

            let res = vbt.search_closest(&p);
            match (res, closest) {
                (Some((p1, d1)), Some((p2, d2))) => {
                    assert_eq!(d1, d2);
                    for i in 0..3 {
                        assert_eq!(p1[i], p2[i]);
                    }
                }
                (None, None) => {}
                v => panic!("result mismatch for point {}: {:?}", p, v),
            }
        }
    }
}