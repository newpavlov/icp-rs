use fxhash::FxHashMap;
use super::Point;

pub struct VoxelBucket {
    map: FxHashMap<[i16; 3], Vec<Point>>,
    radius: f32,
}

#[derive(Debug)]
pub struct ConversionError;

impl VoxelBucket {
    pub fn new(points: &[Point], radius: f32) -> Result<Self, ConversionError> {
        let mut s = Self { map: Default::default(), radius };
        for &point in points {
            let key = convert2key(point, radius)?;
            s.map.entry(key).or_default().push(point);
        }
        Ok(s)
    }

    /// Returns closest point inside the radius if any
    /// and square distance to it
    pub fn search_closest(&self, point: Point) -> Option<(&Point, f32)> {
        let key = match convert2key(point, self.radius) {
            Ok(k) => k,
            Err(_) => return None,
        };
        let mut min_dist2 = self.radius*self.radius;
        let mut closest = None;
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let new_key = [key[0] + i, key[1] + j, key[2] + k];
                    let points = match self.map.get(&key) {
                        Some(v) => v,
                        None => continue,
                    };
                    for p2 in points.iter() {
                        let dx = p2.x - point.x;
                        let dy = p2.y - point.y;
                        let dz = p2.z - point.z;
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
    let v = v/r;
    if std::i16::MIN as f32 <= v && v <= std::i16::MAX as f32 {
        Ok( v as i16)
    } else {
        Err(ConversionError)
    }
}

fn convert2key(p: Point, r: f32) -> Result<[i16; 3], ConversionError> {
    Ok([convert(p.x, r)?, convert(p.x, r)?, convert(p.x, r)?])
}
