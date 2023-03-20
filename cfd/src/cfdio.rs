use ndarray::{Array2, Array3};
use std::fs::File;
use std::io::Write;

pub fn writedatafiles(psi: &Array2<f64>, m:&i64, n:&i64, scalefactor:&i64) {
    let mut vel: Array3<f64> = Array3::zeros((*m as usize, *n as usize,2));
    let mut rgb: Array3<i64> = Array3::zeros((*m as usize, *n as usize,3));

    for i in 0..*m {
        let i_u: usize = i as usize;
        for j in 0..*n {
            let j_u: usize = j as usize;
            vel[[i_u, j_u, 0]] = (psi[[i_u+1, j_u+2]]-psi[[i_u+1, j_u]]) / 2.0;
            vel[[i_u, j_u, 1]] = - (psi[[i_u+2, j_u+1]]-psi[[i_u, j_u+1]]) / 2.0;

            let modvsq = vel[[i_u,j_u, 0]].powi(2) + vel[[i_u, j_u, 1]].powi(2);
            let hue = modvsq.powf(0.4);
            //println!("hue: {}", hue);

            let colours = hue2rgb(hue);
            rgb[[i_u,j_u, 0]] = colours.0;
            rgb[[i_u,j_u, 1]] = colours.1;
            rgb[[i_u,j_u, 2]] = colours.2;
            //println!("rgb: {}, {}, {}", colours.0, colours.1, colours.2);
        }
    }

    let mut colour_file = File::create("colourmap.dat").unwrap();
    let mut velocity_file = File::create("velocity.dat").unwrap();
    for i in 0..*m {
        let i_u: usize = i as usize;
        for j in 0..*n {
            let j_u: usize = j as usize;
            writeln!(&mut colour_file, "{} {} {} {} {}", i+1, j+1, rgb[[i_u,j_u,0]], rgb[[i_u,j_u,1]], rgb[[i_u,j_u,2]]).expect("Unable to write file");
            //fs::write("colourmap.dat", "{},{},{},{},{}", i, j, rgb[[0,i_u,j_u]], rgb[[1,i_u,j_u]], rgb[[2,i_u,j_u]]).expect("Unable to write file");
            if (i%scalefactor == (scalefactor-1)/2) && (j%scalefactor == (scalefactor-1)/2) {
                writeln!(&mut velocity_file, "{} {} {} {}", i+1, j+1, vel[[i_u,j_u,0]], vel[[i_u,j_u, 1]]).expect("Unable to write file");
            }
        }
    }
}

fn hue2rgb(hue:f64) ->  (i64,i64,i64) {
    let rgbmax: f64 = 255.0;

    let r: i64 = (rgbmax*(colfunc(hue-1.0))) as i64;
    let g: i64 = (rgbmax*(colfunc(hue-0.5))) as i64;
    let b: i64 = (rgbmax*(colfunc(hue))) as i64;

    (r,g,b)
}

fn colfunc(x:f64) -> f64 {
    let x1: f64 = 0.2;
    let x2: f64 = 0.5;
    let mut val: f64 = 0.0;

    if x.abs() > x2 {
        val = 0.0;
    } else if x.abs() < x1 {
        val = 1.0;
    } else {
        let frac = (x.abs()-x1)/(x2-x1);
        val = 1.0 - frac.powi(2);
    }
    //println!("Val: {}", val );
    val
}
