use std::{io, fs};
use std::io::Read;
use fxhash::FxHashMap;

type Point = nalgebra::base::Vector3<f32>;

const VOXEL_SIZE: f32 = 0.05;
const PATTERN: &[u8] = b"end_header\n";
const FILTER_PERCENTILE: usize = 20;

fn process_point(map: &mut FxHashMap<[i32; 3], usize>, p: [f32; 3]) {
    let p = [
        (p[0]/VOXEL_SIZE) as i32,
        (p[1]/VOXEL_SIZE) as i32,
        (p[2]/VOXEL_SIZE) as i32,
    ];
    *(map.entry(p).or_default()) += 1;
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path_in = "/media/newpavlov/DATA/ouster_ply/accum_all.ply";
    let path_out = "/media/newpavlov/DATA/ouster_ply/accum_filter.ply";

    let mut data = Vec::new();
    fs::File::open(path_in)?.read_to_end(&mut data)?;
    let n = PATTERN.len();

    let mut map = FxHashMap::default();
    for i in n..data.len() {
        if &data[i - n..i] == PATTERN {
            const POINT_LEN: usize = std::mem::size_of::<[f32; 3]>();
            let payload = &data[i..];
            if payload.len() % POINT_LEN != 0 {
                Err(io::Error::new(io::ErrorKind::InvalidData,
                    "unexpected file length".to_string()))?;
            }
            for chunk in payload.chunks_exact(POINT_LEN) {
                let p: [f32; 3] = unsafe {
                    std::ptr::read_unaligned(chunk.as_ptr() as *const [f32; 3])
                };
                
                process_point(&mut map, p);
            }
        }
    }
    drop(data);

    println!("collected hash map: {}", map.len());

    /*let mut counts: Vec<usize> = map.values().cloned().collect();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let counts_n = counts.len();
    let count_thresh = counts[(FILTER_PERCENTILE*counts_n)/100];
    println!("count_thresh:  {}", count_thresh);
    */
    let count_thresh = 5;


    let points: Vec<Point> = map.into_iter()
        .filter(|(_, c)| *c > count_thresh)
        .map(|(p, _)| {
            Point::new(
                (p[0] as f32)*VOXEL_SIZE,
                (p[1] as f32)*VOXEL_SIZE,
                (p[2] as f32)*VOXEL_SIZE,
            )
        })
        .collect();

    save_scan(path_out, &points)?;

    Ok(())
}

fn save_scan(
    out_path: &str, scan: &[Point],
) -> std::io::Result<()> {
    use std::io::Write;
    let header = format!("ply\n\
        format binary_little_endian 1.0\n\
        element vertex {}\n\
        property float x\n\
        property float y\n\
        property float z\n\
        end_header\n\
    ", scan.len());
    println!("Saving scan with {} points.", scan.len());
    let mut writer = io::BufWriter::new(fs::File::create(out_path)?);
    writer.write_all(header.as_bytes())?;
    for p in scan {
        writer.write_all(&p[0].to_le_bytes())?;
        writer.write_all(&p[1].to_le_bytes())?;
        writer.write_all(&p[2].to_le_bytes())?;
    }
    Ok(())
}
