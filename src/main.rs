use plotters::prelude::*;
use ndarray::{s};

mod automata;

const OUT_FILE_NAME: &'static str = "output/matshow.png";
const TIME_STEPS: usize = 10;
const RELAXATION: f32 = 0.5;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Simulation
    let mut grid = automata::ReacDiffGrid::new(50, 1);
    grid.write(10..20, 10..20, 0, 1.0);
    for _t in 0..TIME_STEPS {
        grid.diffuse(RELAXATION);
    }

    // Plotting
    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 1024)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Species A", ("sans-serif", 40))
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
                        240.0 / 360.0 - 240.0 / 360.0 * (*v as f64 / 20.0),
                        0.7,
                        0.1 + 0.4 * *v as f64,
                    )
                    .filled(),
                )
            }),
    )?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'output' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);

    Ok(())
}
