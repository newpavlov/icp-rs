use std::fmt;
use fxhash::FxHashMap;
use super::Point;

type Map = FxHashMap<[i16; 3], Vec<(Point, u32)>>;

pub struct VoxelBucket {
    map: Map,
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
        let mut map = Map::default();
        for (idx, &point) in points.iter().enumerate() {
            let key = convert2key(point, radius)?;
            let val = (point, idx as u32);
            map.entry(key).or_default().push(val);
        }
        Ok(Self { map, radius })
    }

    pub fn inside_radius(
        &self, p: Point, mut f: impl FnMut(Point, u32, f32),
    ) -> Result<(), ConversionError> {
        let r = self.radius;
        let key = match convert2key(p, r) {
            Ok(k) => k,
            Err(_) => return Err(ConversionError),
        };
        let max_dist2 = r*r;

        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let new_key = [key[0] + i, key[1] + j, key[2] + k];
                    let points = match self.map.get(&new_key) {
                        Some(v) => v,
                        None => continue,
                    };
                    for &(p2, idx) in points.iter() {
                        let d = p2 - p;
                        let d2 = d.x*d.x + d.y*d.y + d.z*d.z;
                        if d2 < max_dist2 {
                            f(p2, idx, d2)
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Returns closest point inside the radius if any
    /// and square distance to it
    pub fn search_closest(&self, p: Point) -> Option<(Point, u32, f32)> {
        let r = self.radius;
        let key = match convert2key(p, r) {
            Ok(k) => k,
            Err(_) => return None,
        };
        let mut min_dist2 = r*r;
        let mut res = None;

        self.inside_radius(p, |p2, idx, d2| {
            if d2 < min_dist2 {
                min_dist2 = d2;
                res = Some((p2, idx, d2));
            }
        });

        res
    }

    pub fn iter_points(&self) -> impl Iterator<Item=&(Point, u32)> {
        self.map.values().flat_map(|v| v.iter())
    }

    /// Collect stored points into a vector
    pub fn get_points(&self) -> Vec<Point> {
        let mut buf = Vec::new();
        for points in self.map.values() {
            buf.extend(points.iter().map(|p| p.0));
        }
        buf
    }

    pub fn remove_points(&mut self, idxs: &[u32]) {
        // TODO optimize
        for val in self.map.values_mut() {
            val.retain(|(_, idx)| idxs.iter().find(|&i| i == idx).is_none())
        }
        self.map.retain(|_, v| !v.is_empty());
        self.map.shrink_to_fit();
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

fn convert2key(p: Point, r: f32) -> Result<[i16; 3], ConversionError> {
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
            for (i, p2) in points.iter().enumerate() {
                let dx = p2.x - p.x;
                let dy = p2.y - p.y;
                let dz = p2.z - p.z;
                let d2 = dx*dx + dy*dy + dz*dz;
                if d2 < min_dist2 {
                    min_dist2 = d2;
                    closest = Some((p2, i as u32, min_dist2));
                }
            }

            let res = vbt.search_closest(p);
            match (res, closest) {
                (Some((p1, i1, d1)), Some((p2, i2, d2))) => {
                    assert_eq!(d1, d2);
                    assert_eq!(i1, i2);
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
