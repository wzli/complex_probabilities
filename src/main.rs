use clap::Parser;
use nalgebra as na;
use plotters::prelude::*;

#[derive(Parser)]
struct Args {
    #[clap(short, long, default_value_t = 0.5)]
    a_self_transition: f64,
    #[clap(short, long, default_value_t = 0.1)]
    b_self_transition: f64,
    #[clap(short, long, default_value_t = 0.1)]
    initial_state: f64,
    #[clap(short, long, default_value_t = 32)]
    nth_root: i32,
    #[clap(short, long, default_value_t = 5)]
    cycles: i32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let m = na::Matrix2::new(
        args.a_self_transition,
        1. - args.b_self_transition,
        1. - args.a_self_transition,
        args.b_self_transition,
    );
    println!("Stochastic Matrix {m}");

    let d = na::Vector2::new(1., args.a_self_transition + args.b_self_transition - 1.);
    println!("Eigenvalues {d}");

    let e01 = (d[0] - m[(0, 0)]) / m[(0, 1)];
    let e0 = na::Vector2::new(1., e01).unscale(1. + e01);
    // second eigen vector is always [1, -1]
    let e1 = na::Vector2::new(1., -1.);
    println!("Eigenvectors {e0} {e1}");
    assert!(((m * e0) - e0 * d[0]).norm() < 0.00000000000001);
    assert!(((m * e1) - e1 * d[1]).norm() < 0.00000000000001);

    let e = na::Matrix2::from_columns(&[e0, e1]);
    let e_inv = e.try_inverse().unwrap();

    // initial state in eigen basis
    let p_init = e_inv * na::Vector2::new(args.initial_state, 1. - args.initial_state);
    let p_init: na::Vector2<na::Complex<f64>> = na::convert(p_init);
    // compute probability state trajectory
    let mut traj = Vec::new();
    for k in 0..(args.nth_root * args.cycles + 1) {
        // compute the diagonal matrix to the power of k
        let theta = if d[1] < 0. {
            (k as f64) * std::f64::consts::PI / (args.nth_root as f64)
        } else {
            0.
        };
        let dk = na::Complex::new(theta.cos(), theta.sin())
            .scale(d[1].abs().powf((k as f64) / (args.nth_root as f64)));
        // compute k'th state in eigen basis
        let p = na::Matrix2::from_diagonal(&[1.0.into(), dk].into()) * p_init;
        // store imaginary part
        let im = p[1].im;
        // transform real part back into original basis
        let p = e * na::Vector2::new(p[0].re, p[1].re);
        let p = (p[0], p[1], im);
        println!("{k} {p:?}");
        traj.push(p);
    }

    // create chart
    let root = SVGBackend::new("plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Complex Probabilities", ("sans-serif", 20).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_3d(0.0..1.0, 0.0..1.0, -0.2..0.2)?;

    chart.configure_axes().draw()?;

    // x axis
    chart
        .draw_series(LineSeries::new([(0., 0., 0.), (0., 1., 0.)], &BLACK))
        .unwrap();

    // y axis
    chart
        .draw_series(LineSeries::new([(0., 0., 0.), (1., 0., 0.)], &BLACK))
        .unwrap();

    // eigen vector 0
    chart
        .draw_series(LineSeries::new([(1., 0., 0.), (0., 1., 0.)], &RED))
        .unwrap();

    // eigen vector 1 / physical probablity space
    chart
        .draw_series(LineSeries::new([(0., 0., 0.), (e0[0], e0[1], 0.)], &GREEN))
        .unwrap();

    // complex probability trajectory line
    chart
        .draw_series(LineSeries::new(traj.clone(), &BLUE))
        .unwrap();

    // complex probability trajectory points
    let i: std::cell::Cell<i32> = std::cell::Cell::new(-1);
    chart.draw_series(PointSeries::of_element(traj, 2, &BLUE, &|c, s, st| {
        i.set(i.get() + 1);
        return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                + Text::new(format!("{:?}", i.get()), (10, 0), ("sans-serif", 8).into_font());
    }))?;

    Ok(())
}
