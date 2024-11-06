/* JAM first attempt at a two pool model with HMM kinetics in Rust
I'll  use  rk4 integration algorithm only because it is what I used to
use and it worked!  These systems are usually not stiff, I suspect
 of the forgiving nature of the HMM equation*/

// This model follows the structure shown in the diagram called
// "Two Pool Model.pdf
const SA:f64=20.0;
const SB:f64= 25.0;
const QA0:f64= 6.0;
const QB0:f64=9.0;
const FOA:f64=7.0;

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
step_size - step size used in method */

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
    let mut stepper = Rk4::new(system, 0.0,  y0, 1.0e1, 1.0e-1);

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
	let  con_a = y[0] / SA;
	let  con_b = y[1] / SB;	

	let FAB  = VAB /  (1.0 + (KAB / con_a));
	let FBA = VBA /  (1.0 + (KBA / con_b));
	let FBO = VBO /  (1.0 + (KBO / con_b));

	// dA/dt = FOA + FBA  - FAB
	//dy[0] = FOA + hmm(  VBA, KBA, y[1],SB) - hmm(VAB,KAB,y[0],SA);
//	dy[0] = FOA + hmm(  VBA, KBA, con_b) - hmm(VAB,KAB,con_a);
	dy[0] = FOA + FBA - FAB;	

	//dB/dt = FAB - FBA - FBO
//	dy[1] = hmm(VAB,KAB,con_a) - hmm(VBA,KBA,con_b) - hmm(VBO,KBO,con_b);
	dy[1] = FAB - FBA - FBO;	

	let total = y[0]+y[1];
	println!("PoolSizes A={:.3}, B={:.3}, Tot={:.3}", y[0], y[1], total);
    }
}

/* HMM equation function for fluxes
fn hmm(vm:f64,km:f64,con:f64) -> f64{
    let flux:f64 = vm / (1.0 + (km /con ));
    return flux;
}
 */

// Initial graph drawing , again from SIR model
/*fn init_graph() -> Figure {
    let mut  fg = Figure::new();
    fg.axes2d()
        .set_title("Pool A Conc", &[])
        .set_legend(Graph(1.0), Graph(1.0), &[], &[])
        .set_x_label("Time", &[])
        .set_y_label("Pool Conc", &[]);
    // fg.set_terminal(&"pngcairo", &"test2.png");
    fg.show();
    fg
}


// update graph as it proceeds, from SIR model

//fn update_graph(fg: &mut Figure, s: &[f32], i: &[f32], r: &[f32], d: &[f32], dt: f32) {
fn update_graph(fg: &mut Figure, con_a: &[f64], con_b: &[f64],  dt: f64) {    
    let x_axis: &Vec<f64> = &(1..=con_a.len() as i64).map(|x| x as f64 * dt).collect();
    fg.clear_axes();
    fg.axes2d()
        .lines(x_axis, con_a, &[Caption("Con A"), Color("blue")])
        .lines(x_axis, con_b, &[Caption("Con B"), Color("red")]);
    fg.show();
}
*/
