use std::{fs, io, fmt, f32};
use icp::{Point, Translation, Rotation};
use std::path::{Path, PathBuf};
use std::convert::TryInto;
use structopt::StructOpt;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::from_args();

    let s1 = load_scan(&cli.scans_dir, cli.scan1)?;
    let mut s2 = load_scan(&cli.scans_dir, cli.scan2)?;

    let mut icp = icp::Icp::new(&s1, cli.corresp_dist, cli.iters_thresh, cli.dist_thresh)?;
    icp.calc_normals();
    let p = cli.initial_pose.0;
    let t = Translation::new(p[0], p[1], p[2]);
    let r = Rotation::from_euler_angles(p[3], p[4], p[5]);

    let (r, t, _, _) = icp.register(&s2, r, t);
    let (roll, pitch, yaw) = r.euler_angles();
    println!(
        "{},{},{},{},{},{}",
        t.vector[0], t.vector[1], t.vector[2],
        roll, pitch, yaw
    );

    if cli.save_result {
        save_scan("scan1.ply", &s1)?;
        save_scan("scan2_before.ply", &s2)?;
        for p in s2.iter_mut() {
            *p = t*r*(*p);
        }
        save_scan("scan2_after.ply", &s2)?;
    }

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
    let mut writer = io::BufWriter::new(fs::File::create(out_path)?);
    writer.write_all(header.as_bytes())?;
    for p in scan {
        writer.write_all(&p[0].to_le_bytes())?;
        writer.write_all(&p[1].to_le_bytes())?;
        writer.write_all(&p[2].to_le_bytes())?;
    }
    Ok(())
}

fn load_scan(dir: &Path, n: u32) -> io::Result<Vec<Point>> {
    use std::io::Read;
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

#[derive(StructOpt)]
#[structopt(
    name = "icp",
    about = "Iterative Closest Point example application")
]
pub struct Cli {
    #[structopt(short = "p", long = "pose", default_value = "0.0,0.0,0.0,0.0,0.0,0.0")]
    /// Initial position estimate (x, y, z, roll, pitch, yaw), distance in meters,
    /// angles in radians
    pub initial_pose: Position,
    #[structopt(short = "d", default_value = "0.25")]
    /// Correspondencies search distance (meters)
    pub corresp_dist: f32,
    #[structopt(short = "i", default_value = "50")]
    /// Maximum number of iterations
    pub iters_thresh: u32,
    #[structopt(short = "e", default_value = "0.0001")]
    /// Distance change threshold
    pub dist_thresh: f32,
    #[structopt(short = "s")]
    /// Save resulting scans?
    pub save_result: bool,
    /// Path to a folder with scans
    pub scans_dir: PathBuf,
    /// First (reference) scan ID
    pub scan1: u32,
    /// Second (matched) scan ID
    pub scan2: u32,
}

#[derive(Copy, Clone, Debug)]
pub struct Position(pub [f32; 6]);

#[derive(Copy, Clone, Debug)]
pub struct PoseParseError;

impl fmt::Display for PoseParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::str::FromStr for Position {
    type Err = PoseParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.split(',')
            .map(|v| f32::from_str(v))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| PoseParseError)?[..]
            .try_into()
            .map(|v| Self(v))
            .map_err(|_| PoseParseError)
    }
}
