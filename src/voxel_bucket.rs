use std::convert::TryFrom;
use fxhash::FxHashMap;
use super::Point;

const VOXEL_SIZE: f32 = 0.2;
const SEARCH_RADIUS: f32 = VOXEL_SIZE;

pub struct VoxelBucket {
    map: FxHashMap<[i16; 3], Vec<Point>>,
}

fn convert(v: f32) -> Result<i16, ConversionError> {
    let v = v/VOXEL_SIZE;
    if std::i16::MIN as f32 <= v && v <= std::i16::MAX as f32 {
        Ok( v as i16)
    } else {
        Err(ConversionError)
    }
}

fn convert2key(p: Point) -> Result<[i16; 3], ConversionError> {
    Ok([convert(p.x)?, convert(p.x)?, convert(p.x)?])
}

#[derive(Debug)]
pub struct ConversionError;

impl VoxelBucket {
    pub fn new(points: &[Point]) -> Result<Self, ConversionError> {
        let mut s = Self { map: Default::default() };
        for &point in points {
            s.map.entry(convert2key(point)?).or_default().push(point);
        }
        Ok(s)
    }

    /// Returns closest point inside the radius if any
    /// and square distance to it
    pub fn search_closest(&self, point: Point) -> Option<(&Point, f32)> {
        let key = match convert2key(point) {
            Ok(k) => k,
            Err(_) => return None,
        };
        let mut min_dist2 = SEARCH_RADIUS*SEARCH_RADIUS;
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
