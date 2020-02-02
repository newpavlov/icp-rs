use std::{fs, io};
use nalgebra::base::{Vector3, Matrix3};
use indicatif::{ProgressBar, ProgressStyle};
//use rayon::prelude::*;

use icp;

type Point = Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = Vector3<f32>;

/*fn convert(p: &[ouster::FullPoint]) -> Vec<[f32; 3]> {
    p.iter().map(|p| p.xyz).collect()
}
*/

// 2x u32 pair indices, 9x f32 rot, 3x f32 trans,
// 1x u32 corresp, 1x f32 sum distance
//const RES_LOG_LEN: usize = 2*4 + 9*4 + 3*4 + 4 + 4;

fn save_result_entry(
    mut w: impl io::Write, entry: IcpResultEntry,
) -> io::Result<()> {
    w.write_all(&entry.ref_idx.to_le_bytes())?;
    w.write_all(&entry.scan_idx.to_le_bytes())?;
    for i in 0..3 {
        for j in 0..3 {
            w.write_all(&entry.r[(i, j)].to_le_bytes())?;
        }
    }
    for i in 0..3 {
        w.write_all(&entry.t[i].to_le_bytes())?;
    }
    w.write_all(&entry.corresp.to_le_bytes())?;
    w.write_all(&entry.sum_dist.to_le_bytes())?;
    Ok(())
}

const PBAR_TEMPLATE: &str = "\
    {wide_bar} {percent:>3}% {pos:>7}/{len} \
    Elapsed: {elapsed_precise} ETA: {eta_precise}\
";

#[derive(Copy, Clone, Debug)]
struct IcpResultEntry {
    ref_idx: u32,
    scan_idx: u32,
    r: Rot,
    t: Trans,
    corresp: u32,
    sum_dist: f32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //let (mut f, start_i, start_j) = get_result_file()?;

    const SCANS: u32 = 4222;
    const MIN_POINTS: usize = 45_000;

    let scan_idxs: Vec<u32> = (1..=SCANS).collect();
    let bar = ProgressBar::new(scan_idxs.len() as u64);
    bar.set_style(ProgressStyle::default_bar().template(PBAR_TEMPLATE));
    bar.tick();

    let path = "/media/newpavlov/DATA/ouster_ply/moscow1_window.bin";
    let mut f = io::BufWriter::new(fs::File::create(path)?);
    for i in 1..=SCANS {
        let ref_scan = load_scan(i)?;
        if ref_scan.len() < MIN_POINTS {
            continue
        }
        let icp = icp::Icp::new(&ref_scan, 0.25, 50, 0.0005)?;
        let mut r = Matrix3::<f32>::identity();
        let mut t = Vector3::<f32>::zeros();

        for j in (i + 1)..=std::cmp::min(i + 50, SCANS) {
            //println!("Processing pair: {} {}", i, j);
            let scan = load_scan(j)?;
            if scan.len() < MIN_POINTS {
                continue;
            }
            let (new_r, new_t, corresp, sum_dist) = icp.register(&scan, r, t);
            r = new_r;
            t = new_t;
            save_result_entry(&mut f, IcpResultEntry {
                ref_idx: i, scan_idx: j,
                r, t, corresp, sum_dist,
            })?;
            if corresp < MIN_POINTS as u32 {
                break;
            }
        }
        bar.inc(1)
    }

    /*let results: Vec<Vec<IcpResultEntry>> = scan_idxs
        .par_iter()
        .progress_with(bar)
        .map(|&i| {
            let ref_scan = load_scan(i).unwrap();
            let mut buf = Vec::new();
            if ref_scan.len() < MIN_POINTS {
                return buf;
            }
            let icp = icp::Icp::new(&ref_scan, 0.20, 50, 0.0005).unwrap();
            let mut r = Matrix3::<f32>::identity();
            let mut t = Vector3::<f32>::zeros();

            for j in (i + 1)..=std::cmp::min(i + 50, SCANS) {
                //println!("Processing pair: {} {}", i, j);
                let scan = load_scan(j).unwrap();
                if scan.len() < MIN_POINTS {
                    continue;
                }
                let (new_r, new_t, corresp, sum_dist) = icp.register(&scan, r, t);
                r = new_r;
                t = new_t;
                buf.push(IcpResultEntry {
                    ref_idx: i, scan_idx: j,
                    r, t, corresp, sum_dist,
                });
                if corresp < MIN_POINTS as u32 {
                    break;
                }
            }
            buf
        })
        .collect();


    let path = "/media/newpavlov/DATA/ouster_ply/moscow1_window.bin";
    let mut f = io::BufWriter::new(fs::File::create(path)?);
    for result in results {
        for entry in result {
            save_result_entry(&mut f, entry)?;
        }
    }
    */

    Ok(())
}

#[allow(dead_code)]
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

    let path = format!("/media/newpavlov/DATA/ouster_ply/moscow1/{}.ply", n);
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
