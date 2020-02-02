use std::{fs, io};
use nalgebra::base::{Vector3, Matrix3};
use indicatif::{ProgressBar, ProgressStyle, ProgressIterator};

use icp;

type Point = Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = Vector3<f32>;

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
    const START_SCAN: u32 = 2;
    const END_SCAN: u32 = 3392;
    const SCANS: u32 = END_SCAN - START_SCAN + 1;
    const KF_DIST_THRESH: f32 = 2.0;

    let s0 = load_scan(START_SCAN)?;

    let mut accum_cloud = Vec::with_capacity(30_000*(SCANS as usize));
    accum_cloud.extend_from_slice(&s0);

    let mut kf_idx = START_SCAN;
    let mut kf_r = Rot::identity();
    let mut kf_t = Trans::zeros();
    let mut kf_icp = icp::Icp::new(&s0, 0.25, 50, 0.0005)?;
    kf_icp.calc_normals();
    let out_path = "/media/newpavlov/DATA/ouster_ply/skolkovo_kf_results.bin";
    let mut f = io::BufWriter::new(fs::File::create(out_path)?);
    let mut r = kf_r;
    let mut t = kf_t;

    drop(s0);

    let bar = ProgressBar::new((SCANS - 1) as u64);
    bar.set_style(ProgressStyle::default_bar().template(PBAR_TEMPLATE));

    for i in ((START_SCAN + 1)..=END_SCAN).progress_with(bar) {
        let mut scan = load_scan(i)?;
        let (new_r, new_t, corresp, sum_dist) = kf_icp.register(&scan, r, t);
        r = new_r;
        t = new_t;
        let global_r = kf_r*r;
        let global_t = kf_t + kf_r*t;

        save_result_entry(&mut f, IcpResultEntry {
            ref_idx: kf_idx, scan_idx: i,
            r, t, corresp, sum_dist,
        })?;

        //println!("================\npair: {}-{}, kf_t: {}, t: {}, global_t: {}",
        //    kf_idx, i, kf_t, t, global_t);

        if t.norm() > KF_DIST_THRESH {
            //println!("\n\nkf update\n\n");
            kf_idx = i;
            kf_r = global_r;
            kf_t = global_t;
            r = Rot::identity();
            t = Trans::zeros();
            kf_icp = icp::Icp::new(&scan, 0.25, 50, 0.0005)?;
            kf_icp.calc_normals();
        }

        for p in scan.iter_mut() {
            *p = global_r*(*p) + global_t;
        }
        accum_cloud.extend_from_slice(&scan);
    }

    save_scan("/media/newpavlov/DATA/ouster_ply/skolkovo_accum_all.ply", &accum_cloud)?;

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

    let path = format!("/media/newpavlov/DATA/ouster_ply/skolkovo/{}.ply", n);
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
