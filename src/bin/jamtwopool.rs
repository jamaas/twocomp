// JAM first attempt at a two pool model with HMM kinetics in Rust
// I'll continue to use  dop853 algorithm although it is likely overkill, it works .
// This is an explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince
//(with stepsize control and dense output). like a heavy duty rk4!

/* This model follows the structure shown in the diagram called
 "Two Pool Model.pdf
Size of Pool A = 20
Size of Pool B  = 25
Qantity of solute in pool A at t0 = 6
Qantity of solute in pool B at t0 = 9
Flux into pool A from outside FOA = 7
For HMM equations:
Flux from Pool A to B
FAB  = VAB /  (1 + (KAB / (y[0] / 20)))
FBA = VBA /  (1 + (KBA / (y[1] / 25)))
FBO = VBO /  (1 + (KBO / (y[1] / 25)))

VAB = 18
VBA = 13
VBO = 8

KAB = 0.32
KBA =  0.36
KB0  =  0.31
 */

/* arguments for dop853 algorithm
 f - Structure implementing the System trait
x - Initial value of the independent variable (usually time)
x_end - Final value of the independent variable
dx - Increment in the dense output. This argument has no effect if the output type is Sparse
y - Initial value of the dependent variable(s)
rtol - Relative tolerance used in the computation of the adaptive step size
atol - Absolute tolerance used in the computation of the adaptive step size
 */

use ode_solvers::dop853::*;
use ode_solvers::*;

type State = Vector2<f64>;
type Time = f64;

fn main() {
    // Initial state. State values of  X, Y, and Z and t0
    let y0 = State::new(6.0, 9.0);

    // Create the structure containing the ODEs.
    let system = TwoPool;

    // Create a stepper and run the integration.
    let mut stepper = Dop853::new(system, 0., 10.0, 0.01, y0, 1.0e-2, 1.0e-6);
    let results = stepper.integrate();

    // Handle result.
    match results {
        Ok(stats) => println!("{}", stats),
        Err(e) => println!("An error occured: {}", e),
    }
}

struct TwoPool;

impl ode_solvers::System<f64, State> for TwoPool {
    fn system(&self, _: Time, y: &State, dy: &mut State) {
	// dy[0] = FOA - (FAB  = VAB /  (1.0 + (KAB / (y[0] / 20)))) + (FBA = VBA /  (1.0 + (KBA / (y[1] / 25))))
	dy[0] = 7.0 - ( 18.0 /  (1.0 + (0.32 / (y[0] / 20.0)))) + (13.0 / ( 1.0 + ( 0.36 /(y[1]/25.0)))) ;
	// dy[1] = (FAB  = VAB /  (1.0 + (KAB / (y[0] / 20)))) - (FBA = VBA /  (1.0 + (KBA / (y[1] / 25)))) -
	//     (FBO  = VBO /  (1 + (KBO / (y[1] / 25.0)))) 
	dy[1] =  (18.0 /  (1.0 + (0.32 / (y[0] / 25.0)))) -   (13.0 / ( 1.0 + ( 0.36 /(y[1]/25.0)))) - (8.0 / ( 1.0 + ( 0.31 /(y[1]/25.0)))) ;  		let total = y[0]+y[1];
	println!("PoolSizes A={:.3}, B={:.3}, Tot={:.3}", y[0], y[1], total);
    }
}
