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
// Last updated on 14/11/2024

//Declaration of initial pool sizes and flux rates

const SA:f64=20.0;
const SB:f64= 25.0;
const QA0:f64= 6.0;
const QB0:f64=9.0;
const QT0:f64=15.0;
 const FOA:f64=7.0;

/* Declaration of kinetic constants for the HMM equations.  A V12
  variable refers to a VMax for the equation describing a flux from
  pool "1" to pool "2".  A K12 variable refers to the "affinity
  constant" for the equation describing a flux from pool "1" to pool
  "2".
 */

const VAB: f64 = 18.0;
const VBA: f64 = 13.0;
const VBO: f64 = 8.0;
const KAB: f64 = 0.32;
const KBA: f64 = 0.36;
const KBO: f64 = 0.31;

// declare the external Rust crates required
use ode_solvers::rk4::*;
use ode_solvers::*;
use gnuplot::*;
use std::thread::sleep;
use std::time::Duration;

/* declare the vector for the State (dependent) variables, and the
 independent variable, usually time (t).*/
type State = Vector3<f64>;
type Time = f64;

fn main() {
    // Initial state. State values of  X, Y, and Z and t0
    let y0 = State::new(QA0, QB0, QT0);

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
            .lines_points(&t[..i], y_out[..i].iter().map(|y| y[1]), &[Caption("B")]);
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
    fn system(&self, _: Time, y: &State, dy: &mut State) {

	/*calculate the concentration of metabolite in each pool at
	 each iteration*/
	let  con_a = y[0] / SA;
	let  con_b = y[1] / SB;	

	/*calculate fluxes corresponding to arrows on diagram, using a
	 HMM equation*/
	let fab  = VAB /  (1.0 + (KAB / con_a));
	let fba = VBA /  (1.0 + (KBA / con_b));
	let fbo = VBO /  (1.0 + (KBO / con_b));

	/* specify the differential equations for each of the state
	variables*/
	dy[0] = FOA + fba - fab;
	dy[1] = fab - fba - fbo;
	/*create a third state varible, only to monitor the size of
	 the total system */
	let total = y[0] + y[1];

	// print the outputs of the model at each integration iteration 
	println!("PoolSizes A={:.3}, B={:.3}, Tot={:.3}", y[0], y[1], total);	
    }
}
