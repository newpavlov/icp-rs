use nalgebra::base::{Vector3, Matrix3};
use std::{fs, io};

type Point = Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = Vector3<f32>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s1 = load_scan(800)?;
    let mut s2 = load_scan(810)?;

    save_scan("scan1.ply", &s1)?;
    save_scan("scan2_before.ply", &s2)?;

    let icp = icp::Icp::new(&s1, 0.25, 50, 0.0005)?;
    let r = Rot::identity();
    let mut t = Trans::zeros();
    t[1] = 3.0;

    let (r, t, _, _) = icp.register(&s2, r, t);

    for p in s2.iter_mut() {
        *p = r*(*p) + t;
    }

    save_scan("scan2_after.ply", &s2)?;

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

fn load_scan(n: u32) -> io::Result<Vec<Point>> {
    use std::io::Read;
    const PATTERN: &[u8] = b"end_header\n";

    let path = format!("clouds/{}.ply", n);
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
                scan.push(Vector3::new(x, y, z));
            }
            return Ok(scan);
        }
    }
    Err(io::Error::new(io::ErrorKind::InvalidData,
        "failed to read PLY header".to_string()))
}
