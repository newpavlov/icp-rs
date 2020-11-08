use std::{io, fs};
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use indicatif::{ProgressBar, ProgressStyle, ProgressIterator};
use icp::{Point, Translation, Rotation};

use structopt::StructOpt;

pub const PBAR_TEMPLATE: &str = "\
    {wide_bar} {percent:>3}% {pos:>7}/{len} \
    Elapsed: {elapsed_precise} ETA: {eta_precise}\
";
pub type DynError = Box<dyn std::error::Error>;

type PointsMap = fxhash::FxHashMap<[i32; 3], smallvec::SmallVec<[u32; 4]>>;

const VOXEL_SIZE: f32 = 0.10;
const COUNT_THRESH: usize = 15;

fn voxel_filter(points: &[Point]) -> Vec<Point> {
    let mut map = PointsMap::default();
    map.reserve(1 << 26);
    let bar = ProgressBar::new(points.len() as u64);
    bar.set_draw_delta((points.len() as u64)/100);
    bar.set_style(ProgressStyle::default_bar().template(PBAR_TEMPLATE));

    assert!(points.len() <= std::u32::MAX as usize);
    for (i, &p) in points.iter().enumerate().progress_with(bar) {
        let key = [
            (p[0]/VOXEL_SIZE).round() as i32,
            (p[1]/VOXEL_SIZE).round() as i32,
            (p[2]/VOXEL_SIZE).round() as i32,
        ];
        map.entry(key).or_default().push(i as u32);
    }

    println!("collected hash map: {}", map.len());

    let bar = ProgressBar::new(map.len() as u64);
    bar.set_draw_delta((map.len() as u64)/100);
    bar.set_style(ProgressStyle::default_bar().template(PBAR_TEMPLATE));

    let res: Vec<_> = map.into_iter()
        .progress_with(bar)
        .filter(|(_, c)| c.len() >= COUNT_THRESH)
        .map(|(p, _)| {
            Point::new(
                (p[0] as f32)*VOXEL_SIZE,
                (p[1] as f32)*VOXEL_SIZE,
                (p[2] as f32)*VOXEL_SIZE,
            )
        })
        .collect();
    println!("filtered points: {}", res.len());
    res
}

fn main() -> Result<(), DynError> {
    let cli = Cli::from_args();

    let rt = load_rt(&cli.traj_path)?;
    println!("loaded poses: {}", rt.len());
    let points = load_scans(&cli.scans_dir, &rt)?;
    println!("loaded points: {}", points.len());

    let fitlered = voxel_filter(&points);
    save_scan(&cli.out_path, &fitlered)?;

    Ok(())
}

#[derive(StructOpt)]
#[structopt(
    name = "assemble",
    about = "Assemble scans into map using provided trajectory")
]
pub struct Cli {
    #[structopt(default_value = "15")]
    /// Voxel filter threshold
    pub voxel_thresh: usize,
    #[structopt(default_value = "0.1")]
    /// Map voxel size
    pub voxel_size: f32,
    #[structopt(short = "o", default_value = "map.ply")]
    /// Output path
    pub out_path: PathBuf,
    /// Path to a folder with scans
    pub scans_dir: PathBuf,
    /// Path to a trajectory file in CSV format
    pub traj_path: PathBuf,
}

fn save_scan(
    out_path: &Path, scan: &[Point],
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
    let mut writer = io::BufWriter::new(fs::File::create(out_path)?);
    writer.write_all(header.as_bytes())?;
    for p in scan {
        writer.write_all(&p[0].to_le_bytes())?;
        writer.write_all(&p[1].to_le_bytes())?;
        writer.write_all(&p[2].to_le_bytes())?;
    }
    Ok(())
}

fn load_rt(path: &Path) -> Result<Vec<(Rotation, Translation)>, DynError> {
    BufReader::new(fs::File::open(path)?);
    panic!();
}

fn load_scans(dir: &Path, rt: &[(Rotation, Translation)]) -> Result<Vec<Point>, DynError> {
    let mut buf = Vec::new();
    for (i, (r, t)) in rt.iter().enumerate() {
        let mut scan = load_scan(dir, i as u32)?;
        for p in scan.iter_mut() {
            *p = t*r*(*p);
        }
        buf.extend_from_slice(&scan);
    }
    Ok(buf)
}

fn load_scan(dir: &Path, n: u32) -> io::Result<Vec<Point>> {
    const PATTERN: &[u8] = b"end_header\n";

    let mut path = dir.to_path_buf();
    path.push(&format!("{}.ply", n));

    let mut data = Vec::new();
    fs::File::open(path)?.read_to_end(&mut data)?;
    let n = PATTERN.len();
    for i in n..data.len() {
        if &data[i - n..i] == PATTERN {
            const POINT_LEN: usize = std::mem::size_of::<[f32; 3]>();
            let payload = &data[i..];
            if payload.len() % POINT_LEN != 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    "unexpected file length".to_string()));
            }
            let mut scan = Vec::with_capacity(payload.len() / POINT_LEN);
            for chunk in payload.chunks_exact(POINT_LEN) {
                let [x, y, z]: [f32; 3] = unsafe {
                    std::ptr::read_unaligned(chunk.as_ptr() as *const [f32; 3])
                };
                scan.push(Point::new(x, y, z));
            }
            return Ok(scan);
        }
    }
    Err(io::Error::new(io::ErrorKind::InvalidData,
        "failed to read PLY header".to_string()))
}