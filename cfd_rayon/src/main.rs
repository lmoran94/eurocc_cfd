extern crate ndarray;
use std::env;
use ndarray::Array2;
use std::time::{Instant, Duration};
mod cfdio;
use crate::cfdio::writedatafiles;
use rayon::prelude::*;
use itertools::{iproduct, Itertools};

fn main() {
    //simulation sizes
    let bbase: f64 = 10.0;
    let hbase: f64 = 15.0;
    let wbase: f64 = 5.0;
    let mbase: i64 = 32;
    let nbase: i64 = 32;
    let mut re: f64 = 0.0;

    //command line parameters parsing
    let args: Vec<String> = env::args().collect();
    if args.len() > 3 {
        re = args[3].parse().unwrap();
    }
    let scalefactor: i64 = args[1].parse().unwrap();
    let niter: i64 = args[2].parse().unwrap();

    println!("Scale Factor = {}, iterations = {}", scalefactor, niter);
    println!("Reynolds number = {}", re);

    //calculate b, h, w, m & n
    let b: f64 = bbase*(scalefactor as f64);
    let h: f64 = hbase*(scalefactor as f64);
    let w: f64 = wbase*(scalefactor as f64);
    let m: i64 = mbase*scalefactor;
    let n: i64 = nbase*scalefactor;
    let re = re / (scalefactor as f64);

    println!("Running CFD on {} x {} grid in parallel", m, n);

    let mdim: usize = (m + 2) as usize;
    let ndim: usize = (n + 2) as usize;
    let m_u: usize = m as usize;
    let n_u: usize = n as usize;

    //allocate arrays of zeros
    let mut psi: Array2<f64> = Array2::zeros((mdim, ndim));
    let mut psi_temp: Array2<f64> = Array2::zeros((mdim, ndim));

    //set boundary conditions
    for i in (b as i64)+1..(b+w) as i64 {
        let i_u: usize = i as usize;
        psi[[i_u, 0]] = (i as f64) - b;
    }
    for i in (b+w) as i64..m+1 {
        let i_u: usize = i as usize;
        psi[[i_u, 0]] = w;
    }
    for j in 1..(h as i64)+1 {
        let j_u: usize = j as usize;
        let m_u: usize = m as usize;
        psi[[m_u+1, j_u]] = w;
    }
    for j in (h as i64)+1..(h+w) as i64 {
        let j_u: usize = j as usize;
        let m_u: usize = m as usize;
        psi[[m_u+1, j_u]] = w-(j as f64)+h;
    }

    //begin iterative Jacobi loop
    println!("Starting main loop...");

    let start = Instant::now();

    for _k in 1..niter+1 {

        //println!("psi before: {}", psi);
        jacobistep(&mut psi_temp, &psi, mdim, ndim);
        //println!("psi_temp: {}", psi_temp);


//        for i in 1..m_u+1 {
//            for j in 1..n_u+1 {
//                psi[[i,j]] = psi_temp[[i,j]];
//            }
//        }
        psi.outer_iter_mut().into_par_iter().enumerate().for_each(|(i, mut view)| {
            if (1..m_u+1).contains(&i) {
                for j in 1..n_u+1 {
                    view[j] = psi_temp[[i,j]]
                }
            }
        });
        //println!("psi: {}", psi);
        if _k%1000==0{
            println!("Completed iteration {}", _k);
        }
    }

    let time_taken = start.elapsed();
    let titer = time_taken.as_secs_f64() / niter as f64;

    //print timings stats
    println!("... finished!");
    println!("Time for {} iterations was {:?}", niter, time_taken);
    println!("Each iteration took {:?}", Duration::from_secs_f64(titer));

    writedatafiles(&psi, &m, &n, &scalefactor);



}

fn jacobistep<'a>(psi_temp: &'a mut Array2<f64>, psi: &Array2<f64>, m:usize, n:usize) -> &'a Array2<f64> {
    psi_temp.outer_iter_mut().into_par_iter().enumerate().for_each(|(i, mut view)| {
        if (1..m-1).contains(&i) {
            for j in 1..n-1 {
                view[j] = 0.25*(psi[[i-1,j]] + psi[[i+1,j]] + psi[[i,j-1]] + psi[[i,j+1]]);
            }
        }
    });
    psi_temp
}

//fn jacobistep<'a>(psi_temp: &'a mut Array2<f64>, psi: &Array2<f64>, m:usize, n:usize) -> &'a Array2<f64> {
//    iproduct!(1..=m, 1..=n).par_bridge().map(|(i,j)| {
//        psi_temp[[i,j]] = 0.25*(psi[[i-1,j]] + psi[[i+1,j]] + psi[[i,j-1]] + psi[[i,j+1]])
//    })
//}

//fn jacobistep<'a>(psi_temp: &'a mut Array2<f64>, psi: & Array2<f64>, m:usize, n:usize) -> &'a Array2<f64>{
//    for i in 1..m-1 {
//        //let j_u: usize = j as usize;
//        for j in 1..n-1 {
//            //let i_u: usize = i as usize;
//            psi_temp[[i,j]] = 0.25*(psi[[i-1,j]]+psi[[i+1,j]]+psi[[i,j-1]]+psi[[i,j+1]]);
//        }
//    }
//    psi_temp
//}

//    let row_major = psi_temp.iter().enumerate();
//    for elem in row_major {
//        println!("elem! {:?}", elem);
        //psi_temp[elem] = 0.25*psi[elem-1];
//    }
//    psi_temp
    //println!("row_major: {:?}", row_major)
