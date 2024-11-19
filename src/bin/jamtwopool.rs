/* JAM first attempt at a toy two pool model with HMM kinetics in
 Rust.  This generalised model should be expanadable to any number of
 pools and interactions.  I'll use rk4 integration algorithm only
 because it is what I have used historically and it worked!  These
 biological systems, comprised of Henri-Michaelis-Menten (HMM) kinetic
 equations are usually not stiff.  This model follows the structure of
the accompanying diagram called "Two Pool Model.pdf */

/* Aknowlegements to Dr. Pavel Sountsov for help with the graphing,
 * Dr. Nikhil Suresh for help with compartmental model structures in
 * Rust, and Dr. Sylvain Renevey for the ode_solver crate */

// First produced by JAM in Norwich UK on 14/11/2024
// Last updated on 19/11/2024

// declare the external Rust crates required
use ode_solvers::rk4::*;
use ode_solvers::*;
use gnuplot::*;
use std::thread::sleep;
use std::time::Duration;

//Declaration of initial pool sizes and flux rates

//SA=20.0, SB=25.0,QA0= 6.0,QB0=9.0,
//QT0=15.0, FOA=7.0
static I:[f64;6] = [20.0, 25.0,6.0,9.0,15.0,7.0];

/* Declaration of kinetic constants for the HMM equations.  A V12
  variable refers to a VMax for the equation describing a flux from
  pool "1" to pool "2".  A K12 variable refers to the "affinity
  constant" for the equation describing a flux from pool "1" to pool
  "2". */

// VAB = 18.0, VBA = 13.0, VBO = 8.0, KAB = 0.32,
// KBA = 0.36, KBO = 0.31
static C: [f64; 6] = [18.0,13.0,8.0,0.32,0.36,0.31];

/* declare the vector for the State (dependent) variables, and the
 independent variable, usually time (t).*/
type State = Vector3<f64>;
type Time = f64;

fn main() {
    // Initial state. State values of  QA0, QB0, QT0
        let y0 = State::new(I[2], I[3], I[4]);

    /* declare the system structure (function) name containing the
    system of ODEs */
    let system = TwoPool;

    /* the five arguments required for 4th-order Runge-Kutta (RK4) algorithm
    f - Structure implementing the System trait
    x - Initial value of the independent variable (usually time)
    y - Initial value of the dependent variable(s)
    x_end - Final value of the independent variable
    step_size - step size used in method */

    /* declare the stepper function and input the necessary RK4Cr
     method arguments .*/
    let mut stepper = Rk4::new(system, 0.0,  y0, 1.0e1, 1.0e-1);

    // Execute the stepper function, do integrations and produce results
    let  results = stepper.integrate();

    // declare vector to collect output data, primarily for graphing 
    let (t, y_out) = stepper.results().get();
    // create new figure (graph)
    let mut fg = Figure::new();

    // collect initial values to build graph
    let y_min = y_out.iter().fold(y_out[0].min(), |x, y| x.min(y.min()));
    let y_max = y_out.iter().fold(y_out[0].max(), |x, y| x.max(y.max()));

    // start iteratively producing figure using gnuplot code
    for i in 0..t.len()
    {
	fg.clear_axes();
	fg.axes2d()
            .set_x_range(Fix(0.), Fix(t[t.len() - 1]))
            .set_y_range(Fix(y_min), Fix(y_max))
            .lines_points(&t[..i], y_out[..i].iter().map(|y| y[0]),&[Caption("A")])
            .lines_points(&t[..i], y_out[..i].iter().map(|y| y[1]), &[Caption("B")])
            .lines_points(&t[..i], y_out[..i].iter().map(|y| y[2]), &[Caption("Total")]);	
	fg.show_and_keep_running().unwrap();
		sleep(Duration::from_millis(10));
    }

    // Handle result.
    match results {
        Ok(stats) => println!("{}", stats),
        Err(e) => println!("An error occured: {}", e),
    }
}

//declare system function, used by the stepper above
struct TwoPool;

impl ode_solvers::System<f64, State> for TwoPool {
    fn system(&self, time: Time, y: &State, dy: &mut State) {

	/*calculate the concentration of metabolite in each pool at
	 each iteration*/
	let  con_a = y[0] / I[0];
	let  con_b = y[1] / I[1];	

	/*calculate fluxes corresponding to arrows on diagram, using a
	 HMM equation*/
	let fab  = C[0] /  (1.0 + (C[3] / con_a));
	let fba =  C[1] /  (1.0 + (C[4] / con_b));
	let fbo = C[2] /  (1.0 + (C[5] / con_b));

	/* specify the differential equations for each of the state
	variables*/
	dy[0] = I[5] + fba - fab;
	dy[1] = fab - fba - fbo;

	//This is a very crude way of of getting a value for the size
	// of the total system, don't like recalculating values that
	// are already calculated!
	dy[2] = I[5] - fbo;

	// print the outputs of the model at each integration iteration 
	println!("PoolSizes t={:.3}, A={:.3}, B={:.3}, Tot={:.3}", time, y[0], y[1], y[2]);	
    }
}
