use nalgebra::base::Vector3;
use std::{fs, io};

type Point = Vector3<f32>;

static SCANS_PATH: &str = "/media/newpavlov/DATA/ouster_ply/moscow3/";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = load_scan(1)?;
    let normals = vec![Point::new(0.0, 0.0, 1.0); s.len()];
    let mut icp_state = icp::Icp::new(&s, 0.5, 50, 0.001)?;
    icp_state.calc_normals();

    let r = icp_state.get_points_normals();

    let points: Vec<Point> = r.iter().map(|&(p, _)| p).collect();
    let norm_points: Vec<Point> = r.iter().map(|&(p, n)| p + 0.1*n).collect();
    save_scan("test.ply", &points)?;
    save_scan("test_normals.ply", &norm_points)?;
    Ok(())
}

fn load_scan(n: u32) -> io::Result<Vec<Point>> {
    use std::io::Read;
    const PATTERN: &[u8] = b"end_header\n";

    let path = format!("{}{}.ply", SCANS_PATH, n);
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

/*fn save_scan_normals(
    out_path: &str, scan: &[(Point, Point)]
) -> std::io::Result<()> {
    use std::io::Write;
    let header = format!("ply\n\
        format binary_little_endian 1.0\n\
        element vertex {}\n\
        property float x\n\
        property float y\n\
        property float z\n\
        property float nx\n\
        property float ny\n\
        property float nz\n\
        end_header\n\
    ", scan.len());
    println!("Saving scan with normals with {} points.", scan.len());
    let mut writer = io::BufWriter::new(fs::File::create(out_path)?);
    writer.write_all(header.as_bytes())?;
    for (p, n) in scan.iter() {
        writer.write_all(&p[0].to_le_bytes())?;
        writer.write_all(&p[1].to_le_bytes())?;
        writer.write_all(&p[2].to_le_bytes())?;

        writer.write_all(&n[0].to_le_bytes())?;
        writer.write_all(&n[1].to_le_bytes())?;
        writer.write_all(&n[2].to_le_bytes())?;
    }
    Ok(())
}
*/
