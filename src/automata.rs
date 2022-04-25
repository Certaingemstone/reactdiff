use rand::Rng;
use ndarray::{array, s, Array2, Array3};

// Store the concentrations of A,B,C on a 2D arrays
pub struct ReacDiffGrid {
    pub width: usize,
    pub species: usize,
    pub a: Array3<f32>, // first two dimensions are position, third is species
}

// Initialize
impl ReacDiffGrid {
    // create with zeroes
    pub fn new(width: usize, species: usize) -> Self {
        Self{a: Array3::<f32>::zeros([width, width, species]), width: width, species: species}
    }

    // create with random values between 0 and 1
    pub fn random(width: usize, species: usize) -> Result<Self, ()> {
        let mut rng = rand::thread_rng();
        let mut grid = ReacDiffGrid::new(width, species);
        for elem in grid.a.iter_mut() {
            *elem += rng.gen::<f32>();
        }
        Ok(grid)
    }

    // fill a block
    pub fn write(&mut self, xrange: std::ops::Range<usize>, yrange: std::ops::Range<usize>, spec_idx: usize, value: f32) {
        for (x, y) in xrange.zip(yrange) {
            self.a[(x,y, spec_idx)] = value
        }
    }
}

// Evolution
impl ReacDiffGrid {

    // Apply ball reaction model to each cell
    // params: k1, k2, k3, dt
    pub fn ball_model(&mut self, params: [f32; 4]) {
        // du/dt = u(k1*v - k3*w)
        // dv/dt = v(k2*w - k1*u)
        // dw/dt = w(k3*u - k2*v)
        assert_eq!(self.a.shape()[2], 3);
        // for each cell, calculate evolution and add it
        for mut x in self.a.rows_mut() {
            // TODO: vectorize
            let du = x[0] * (params[0]*x[1] - params[2]*x[2]) * params[3];
            let dv = x[1] * (params[1]*x[2] - params[2]*x[0]) * params[3];
            let dw = x[2] * (params[2]*x[0] - params[2]*x[1]) * params[3];
            x += &ndarray::ArrayView::from(&array![du, dv, dw])
        }
    }

    // Apply diffusion, with relaxation towards local average defined by relax
    // e.g. relax 0.5 -> if currently 2 and average is 1, next value is 1.5
    pub fn diffuse(&mut self, relax: f32) {
        assert!(relax > 0. && relax < 1.);
        // calculate averages with hard boundary condition
        let mut d: Array2<f32> = Array2::<f32>::zeros((self.width, self.width));
        for i in 0..self.species {
            // TODO: Get a convolution to work with periodic BC and get rid of the naive implementation
            let avg = self.a.sum() / (self.width * self.width) as f32;
            let mut species_slice = self.a.slice_mut(s![.., .., i]);
            // Calculate averages over 3x3 windows except on edges 
            for j in 0..self.width-2 {
                for k in 0..self.width-2 {
                    let window = species_slice.slice(s![j..j+3, k..k+3]);
                    d[(j+1, k+1)] = window.sum() / 9.
                }
            }
            // Pad averages on the edges with the overall average
            d.slice_mut(s![0, ..]).fill(avg);
            d.slice_mut(s![-1, ..]).fill(avg);
            d.slice_mut(s![.., 0]).fill(avg);
            d.slice_mut(s![.., -1]).fill(avg);
            // Do the step
            species_slice += &ndarray::ArrayView::from(&d);
            d = Array2::<f32>::zeros((self.width, self.width)); // reset
        }    
    }
}

// Misc
/*
impl ReacDiffGrid {
    // Compute the average of each species
    pub fn means(&self) -> Vec<f32> {

    }

    // Compute the sum of each species
    pub fn sums(&self) -> Vec<f32> {
        for elem in self.a {

        }
    }
}
*/

fn wrap(idx: usize, width: usize) -> usize {
    ((idx % width) + width) % width
}