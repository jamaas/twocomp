// JAM first attempt at a two pool model with HMM kinetics in Rust
// I'll continue to use  dop853 algorithm although it is likely overkill, it works .
// This is an explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince
//(with stepsize control and dense output). like a heavy duty rk4!

// This model follows the structure shown in the diagram called
// "Two Pool Model.pdf
const SA:f64=20.0;
const SB:f64= 25.0;
const QA0:f64= 6.0;
const QB0:f64=9.0;
const FOA:f64=7.0;

/*
For HMM equations:
Flux from Pool A to B
FAB  = VAB /  (1 + (KAB / (y[0] / 20)))


*/

const VAB: f64 = 18.0;
const VBA: f64 = 13.0;
const VBO: f64 = 8.0;
const KAB: f64 = 0.32;
const KBA: f64 = 0.36;
const KBO: f64 = 0.31;

/* arguments for Rk4 algorithm
 f - Structure implementing the System trait
x - Initial value of the independent variable (usually time)
y - Initial value of the dependent variable(s)
x_end - Final value of the independent variable
step_size - step size used in method
 */

//use ode_solvers::dop853::*;
use ode_solvers::rk4::*;
use ode_solvers::*;

type State = Vector2<f64>;
type Time = f64;

fn main() {
    // Initial state. State values of  X, Y, and Z and t0
    let y0 = State::new(QA0, QB0);

    // Create the structure containing the ODEs.
    let system = TwoPool;

    // Create a stepper and run the integration.
//    let mut stepper = Dop853::new(system, 0., 10.0, 0.01, y0, 1.0e-2, 1.0e-6);
    let mut stepper = Rk4::new(system, 0.0,  y0, 1.0e1, 1.0e-2);    
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
	//FAB  = VAB /  (1 + (KAB / (y[0] / SA)))
	//FBA = VBA /  (1 + (KBA / (y[1] / SB)))
	//FBO = VBO /  (1 + (KBO / (y[1] / SB)))

	let  con_a = y[0] / SA;
	let  con_b = y[1] / SB;	

	// dA/dt = FOA + FBA  - FAB
	//dy[0] = FOA + hmm(  VBA, KBA, y[1],SB) - hmm(VAB,KAB,y[0],SA);
	dy[0] = FOA + hmm(  VBA, KBA, con_b) - hmm(VAB,KAB,con_a);

	//dB/dt = FAB - FBA - FBO
	dy[1] = hmm(VAB,KAB,con_a) - hmm(VBA,KBA,con_b) - hmm(VBO,KBO,con_b);

	let total = y[0]+y[1];
	println!("PoolSizes A={:.3}, B={:.3}, Tot={:.3}", y[0], y[1], total);
    }
}

// HMM equation function for fluxes
fn hmm(vm:f64,km:f64,con:f64) -> f64{
    let flux:f64 = vm / (1.0 + (km /con ));
    return flux;
}
