use plotters::prelude::*;
use ndarray::{s};

extern crate reactdiff; // import src/lib.rs

const OUT_FILE_NAME: &'static str = "output/matshow.png";
const DT: f32 = 0.8;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // Take input
    let mut buffer = String::new();
    let mut stdin = std::io::stdin();
    println!("Enter number of time steps");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let time_steps: usize = buffer.trim().parse().unwrap();
        buffer.clear();
    println!("Enter diffusion rate (0,1)");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let relaxation: f32 = buffer.trim().parse().unwrap();
        buffer.clear();
    println!("Enter k1");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let k1: f32 = buffer.trim().parse().unwrap();
        buffer.clear();
    println!("Enter k2");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let k2: f32 = buffer.trim().parse().unwrap();
        buffer.clear();
    println!("Enter k3");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let k3: f32 = buffer.trim().parse().unwrap();
        buffer.clear();
    println!("Diffusion only? (true, false)");
        stdin.read_line(&mut buffer).expect("failed to read input");
        let diff_only: bool = buffer.trim().parse().unwrap();
        buffer.clear();

    // Simulation (this is the key part)
    let mut grid = reactdiff::ReacDiffGrid::random(300, 3);
    for t in 0..time_steps {
        grid.diffuse(relaxation);
        if !diff_only {grid.ball_model([k1, k2, k3, DT]);}
        println!("Time: {}", t);
    }

    // Plotting
    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 1024)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Species A, t={}", time_steps), ("sans-serif", 40))
        .margin(2)
        .build_cartesian_2d(0i32..grid.width as i32, grid.width as i32..0i32)?;

    chart
        .configure_mesh()
        .draw()?;

    let concentrations = grid.a.slice_mut(s![.., .., 0]);

    chart.draw_series(
        concentrations
            .rows()
            .into_iter()
            .zip(0..) // (value: &f32, index: usize)
            .map(|(row, y)| row.into_iter().zip(0..).map(move |(val, x)| (x as i32, y as i32, val)))
            .flatten()
            .map(|(x, y, v)| {
                Rectangle::new(
                    [(x, y), (x + 1, y + 1)],
                    HSLColor(
                        240.0/360.0 - 240.0 / 360.0 * (*v as f64 / 2.0),
                        0.7,
                        0.1 + 0.4 * *v as f64,
                    )
                    .filled(),
                )
            }),
    )?;

    root.present().expect("Unable to write result to file, please make sure 'output' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);

    Ok(())
}
