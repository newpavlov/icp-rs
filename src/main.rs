use std::{fs, io};
use nalgebra::base::{Vector3, Matrix3};

use icp;

type Point = Vector3<f32>;
type Rot = Matrix3<f32>;
type Trans = Vector3<f32>;

/*fn convert(p: &[ouster::FullPoint]) -> Vec<[f32; 3]> {
    p.iter().map(|p| p.xyz).collect()
}
*/

// 2x u32 pair indices, 9x f32 rot, 3x f32 trans,
// 1x u32 corresp, 4x u8 sync bytes
const RES_LOG_LEN: usize = 2*4 + 9*4 + 3*4 + 4 + 4;
const SYNC_BYTES: [u8; 4] = [0xFF; 4];

fn save_result(
    mut w: impl io::Write,
    i: u32, j: u32,
    r: Rot, t: Trans,
    corresp: u32,
) -> io::Result<()> {
    w.write_all(&i.to_le_bytes())?;
    w.write_all(&j.to_le_bytes())?;
    for i in 0..3 {
        for j in 0..3 {
            w.write_all(&r[(i, j)].to_le_bytes())?;
        }
    }
    for i in 0..3 {
        w.write_all(&t[i].to_le_bytes())?;
    }
    w.write_all(&corresp.to_le_bytes())?;
    // sync bytes
    w.write_all(&SYNC_BYTES)?;
    Ok(())
}

fn get_result_file() -> io::Result<impl io::Write> {
    use std::io::{Read, Write};

    let mut old_data = Vec::new();
    fs::File::open("icp_results.bin")?.read_to_end(&mut old_data)?;
    let n = old_data.len()/RES_LOG_LEN;
    for chunk in old_data.chunks_exact(RES_LOG_LEN) {
        if &chunk[RES_LOG_LEN - 4..] != &SYNC_BYTES {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "didn't get sync bytes".to_string()));
        }
    }
    let mut f = io::BufWriter::new(fs::File::create("icp_results.bin")?);
    f.write_all(&old_data[..n*RES_LOG_LEN])?;
    Ok(f)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut f = get_result_file()?;
    for i in 1..=850 {
        let ref_scan = load_scan(i)?;
        let icp = icp::Icp::new(&ref_scan, 0.25, 50, 0.001);
        let mut r = Matrix3::<f32>::identity();
        let mut t = Vector3::<f32>::zeros();
        for j in (i + 1)..=850 {
            println!("-------------------------\n Processing pair: {} {}", i, j);
            let scan = load_scan(j)?;
            let (new_r, new_t, corresp) = icp.register(&scan, r, t);
            r = new_r;
            t = new_t;
            if corresp < 40_000 {
                break;
            }
            save_result(&mut f, i, j, r, t, corresp)?;
        }
    }

    /*
    let s800 = load_scan(800)?;
    //let mut s2 = load_scan(825)?;

    //let (s1, mut s2) = get_ellipses();
    //println!("got scans: {} and {}", s1.len(), s2.len());

    let mut r = Matrix3::<f32>::identity();
    let mut t = Vector3::<f32>::zeros();
    //t[1] = 7.0;

    let icp = icp::Icp::new(&s800, 0.5, 50, 0.001);
    save_scan(&format!("test800.ply"), &s800)?;

    /*let res = icp.register(&s2, r, t);
    r = res.0;
    t = res.1;

    for p in s2.iter_mut() {
        *p = r*(*p) + t;
    }

    save_scan(&format!("test1.ply"), &s1)?;
    save_scan(&format!("test2.ply"), &s2)?;
    */



    for i in 801..=850 {
        println!("-------------------------\n Scan: {}", i);
        let mut s = load_scan(i)?;

        let res = icp.register(&s, r, t);
        r = res.0;
        t = res.1;

        for p in s.iter_mut() {
            *p = r*(*p) + t;
        }

        save_scan(&format!("test{}.ply", i), &s)?;
    }


    /*
    let s2 = &mut ss[n-1];
    let mut r = Matrix3::<f32>::identity();
    let mut t = Vector3::<f32>::zeros();
    t[1] = 15.65;

    println!("{} {}", r, t);

    for p in s2.iter_mut() {
        let pt = nalgebra::base::Vector3::new(p[0], p[1], p[2]);
        let pt = r*pt + t;
        *p = [pt[0], pt[1], pt[2]];
    }

    let path = format!("test2.ply");
    save_scan(&path, &s2)?;
    */
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

/*
fn get_ellipses() -> (Vec<Point>, Vec<Point>) {
    let s1: Vec<Point> = [
       [  5.00000000e-01,   0.00000000e+00],
       [  4.99747763e-01,   3.17599050e-02],
       [  4.98985182e-01,   6.36800766e-02],
       [  4.97694392e-01,   9.59227316e-02],
       [  4.95844775e-01,   1.28653943e-01],
       [  4.93391648e-01,   1.62045444e-01],
       [  4.90274322e-01,   1.96276230e-01],
       [  4.86413429e-01,   2.31533805e-01],
       [  4.81707398e-01,   2.68014798e-01],
       [  4.76027881e-01,   3.05924543e-01],
       [  4.69213988e-01,   3.45474940e-01],
       [  4.61065124e-01,   3.86879575e-01],
       [  4.51332356e-01,   4.30344532e-01],
       [  4.39708411e-01,   4.76052575e-01],
       [  4.25816853e-01,   5.24137415e-01],
       [  4.09201876e-01,   5.74643629e-01],
       [  3.89321641e-01,   6.27466842e-01],
       [  3.65550666e-01,   6.82268894e-01],
       [  3.37200296e-01,   7.38365655e-01],
       [  3.03570472e-01,   7.94594157e-01],
       [  2.64048155e-01,   8.49184483e-01],
       [  2.18262863e-01,   8.99691775e-01],
       [  1.66289623e-01,   9.43075313e-01],
       [  1.08850489e-01,   9.76015514e-01],
       [  4.74211434e-02,   9.95492311e-01],
       [ -1.58599763e-02,   9.99496796e-01],
       [ -7.85137436e-02,   9.87594233e-01],
       [ -1.38178828e-01,   9.61054861e-01],
       [ -1.93021761e-01,   9.22480569e-01],
       [ -2.41946689e-01,   8.75126962e-01],
       [ -2.84576930e-01,   8.22231040e-01],
       [ -3.21087206e-01,   7.66558560e-01],
       [ -3.51992280e-01,   7.10215276e-01],
       [ -3.77964473e-01,   6.54653671e-01],
       [ -3.99707172e-01,   6.00780083e-01],
       [ -4.17881606e-01,   5.49090022e-01],
       [ -4.33072606e-01,   4.99792430e-01],
       [ -4.45778559e-01,   4.52908274e-01],
       [ -4.56414440e-01,   4.08342301e-01],
       [ -4.65320741e-01,   3.65932276e-01],
       [ -4.72774244e-01,   3.25481268e-01],
       [ -4.78998541e-01,   2.86777946e-01],
       [ -4.84173389e-01,   2.49608730e-01],
       [ -4.88442601e-01,   2.13764595e-01],
       [ -4.91920496e-01,   1.79044418e-01],
       [ -4.94697041e-01,   1.45256157e-01],
       [ -4.96841880e-01,   1.12216687e-01],
       [ -4.98407415e-01,   7.97508264e-02],
       [ -4.99431095e-01,   4.76898873e-02],
       [ -4.99937032e-01,   1.58699588e-02],
       [ -4.99937032e-01,  -1.58699588e-02],
       [ -4.99431095e-01,  -4.76898873e-02],
       [ -4.98407415e-01,  -7.97508264e-02],
       [ -4.96841880e-01,  -1.12216687e-01],
       [ -4.94697041e-01,  -1.45256157e-01],
       [ -4.91920496e-01,  -1.79044418e-01],
       [ -4.88442601e-01,  -2.13764595e-01],
       [ -4.84173389e-01,  -2.49608730e-01],
       [ -4.78998541e-01,  -2.86777946e-01],
       [ -4.72774244e-01,  -3.25481268e-01],
       [ -4.65320741e-01,  -3.65932276e-01],
       [ -4.56414440e-01,  -4.08342301e-01],
       [ -4.45778559e-01,  -4.52908274e-01],
       [ -4.33072606e-01,  -4.99792430e-01],
       [ -4.17881606e-01,  -5.49090022e-01],
       [ -3.99707172e-01,  -6.00780083e-01],
       [ -3.77964473e-01,  -6.54653671e-01],
       [ -3.51992280e-01,  -7.10215276e-01],
       [ -3.21087206e-01,  -7.66558560e-01],
       [ -2.84576930e-01,  -8.22231040e-01],
       [ -2.41946689e-01,  -8.75126962e-01],
       [ -1.93021761e-01,  -9.22480569e-01],
       [ -1.38178828e-01,  -9.61054861e-01],
       [ -7.85137436e-02,  -9.87594233e-01],
       [ -1.58599763e-02,  -9.99496796e-01],
       [  4.74211434e-02,  -9.95492311e-01],
       [  1.08850489e-01,  -9.76015514e-01],
       [  1.66289623e-01,  -9.43075313e-01],
       [  2.18262863e-01,  -8.99691775e-01],
       [  2.64048155e-01,  -8.49184483e-01],
       [  3.03570472e-01,  -7.94594157e-01],
       [  3.37200296e-01,  -7.38365655e-01],
       [  3.65550666e-01,  -6.82268894e-01],
       [  3.89321641e-01,  -6.27466842e-01],
       [  4.09201876e-01,  -5.74643629e-01],
       [  4.25816853e-01,  -5.24137415e-01],
       [  4.39708411e-01,  -4.76052575e-01],
       [  4.51332356e-01,  -4.30344532e-01],
       [  4.61065124e-01,  -3.86879575e-01],
       [  4.69213988e-01,  -3.45474940e-01],
       [  4.76027881e-01,  -3.05924543e-01],
       [  4.81707398e-01,  -2.68014798e-01],
       [  4.86413429e-01,  -2.31533805e-01],
       [  4.90274322e-01,  -1.96276230e-01],
       [  4.93391648e-01,  -1.62045444e-01],
       [  4.95844775e-01,  -1.28653943e-01],
       [  4.97694392e-01,  -9.59227316e-02],
       [  4.98985182e-01,  -6.36800766e-02],
       [  4.99747763e-01,  -3.17599050e-02],
       [  5.00000000e-01,  -1.22464680e-16],
    ].iter().map(|p| Point::new(p[0], p[1], 0.0)).collect();
    let s2: Vec<Point> = [
       [ 0.23879128,  0.06028723],
       [ 0.25379643,  0.0882801 ],
       [ 0.26843055,  0.11665829],
       [ 0.28275573,  0.14557271],
       [ 0.29682471,  0.17518381],
       [ 0.31068063,  0.2056637 ],
       [ 0.32435603,  0.23719857],
       [ 0.33787116,  0.26999101],
       [ 0.35123115,  0.30426228],
       [ 0.36442181,  0.34025412],
       [ 0.37740352,  0.37822961],
       [ 0.39010266,  0.41847237],
       [ 0.40239956,  0.4612826 ],
       [ 0.4141122 ,  0.506968  ],
       [ 0.42497431,  0.55582638],
       [ 0.43460726,  0.6081154 ],
       [ 0.44248551,  0.66400322],
       [ 0.44789802,  0.72349296],
       [ 0.44991245,  0.78631439],
       [ 0.44735688,  0.85178254],
       [ 0.43884478,  0.91863807],
       [ 0.4228789 ,  0.98491302],
       [ 0.39806726,  1.04790296],
       [ 0.36345205,  1.10434849],
       [ 0.31888041,  1.15089179],
       [ 0.26526585,  1.18474464],
       [ 0.2045756 ,  1.20433697],
       [ 0.13949091,  1.20965145],
       [ 0.07286821,  1.20209242],
       [ 0.00723002,  1.18399158],
       [-0.05554119,  1.15800907],
       [-0.11427278,  1.12665583],
       [-0.16840695,  1.09202663],
       [-0.21783734,  1.05571847],
       [-0.26274673,  1.01886395],
       [-0.30347783,  0.98221494],
       [-0.34044371,  0.94623519],
       [-0.3740717 ,  0.91118203],
       [-0.40477163,  0.87717082],
       [-0.43292009,  0.84422243],
       [-0.4588544 ,  0.81229673],
       [-0.4828721 ,  0.78131546],
       [-0.50523332,  0.75117736],
       [-0.5261645 ,  0.72176794],
       [-0.54586238,  0.69296551],
       [-0.56449799,  0.66464467],
       [-0.58222022,  0.63667809],
       [-0.59915907,  0.60893718],
       [-0.61542827,  0.58129184],
       [-0.63112756,  0.55360978],
       [-0.64634448,  0.52575538],
       [-0.66115577,  0.49758821],
       [-0.67562824,  0.46896131],
       [-0.68981932,  0.43971928],
       [-0.70377701,  0.40969612],
       [-0.71753932,  0.37871299],
       [-0.73113292,  0.34657578],
       [-0.74457092,  0.31307282],
       [-0.75784944,  0.27797281],
       [-0.77094246,  0.24102336],
       [-0.78379465,  0.20195086],
       [-0.79631108,  0.16046266],
       [-0.80834328,  0.11625322],
       [-0.81967022,  0.06901695],
       [-0.82997339,  0.01847129],
       [-0.83880536, -0.0356043 ],
       [-0.84555272, -0.09330682],
       [-0.84939763, -0.15451845],
       [-0.84928828, -0.21878102],
       [-0.84393831, -0.28514218],
       [-0.83188641, -0.35200074],
       [-0.81165327, -0.4170133 ],
       [-0.78201757, -0.47715853],
       [-0.74238019, -0.52905398],
       [-0.69310273, -0.56953728],
       [-0.63564847, -0.5963616 ],
       [-0.57240147, -0.6087199 ],
       [-0.50620152, -0.60734994],
       [-0.43979153, -0.5941946 ],
       [-0.37539667, -0.57182092],
       [-0.31454058, -0.54286141],
       [-0.25807025, -0.50963926],
       [-0.20629624, -0.47400161],
       [-0.15916175, -0.4373047 ],
       [-0.1163904 , -0.40047906],
       [-0.07759542, -0.36412133],
       [-0.04235133, -0.32858288],
       [-0.01023675, -0.29404312],
       [ 0.01914276, -0.26056516],
       [ 0.0461445 , -0.22813595],
       [ 0.07108573, -0.19669397],
       [ 0.09424487, -0.16614794],
       [ 0.11586472, -0.13638905],
       [ 0.13615636, -0.10729863],
       [ 0.15530318, -0.07875281],
       [ 0.17346474, -0.05062511],
       [ 0.19078011, -0.02278752],
       [ 0.20737084,  0.00488924],
       [ 0.22334341,  0.03253622],
       [ 0.23879128,  0.06028723],
    ].iter().map(|p| Point::new(p[0], p[1] - 0.2, 0.1)).collect();
    (s1, s2)
}
*/
