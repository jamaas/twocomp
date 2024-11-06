fn main() {
  use gnuplot::{Figure, Caption, Color};

  println!("Hello, world!");
    
    //let x = [0u32, 1, 2];
    let x: [u32; 3] = [0,1,2];
    let y = [3u32, 4, 5];
    let mut fg = Figure::new();
    fg.axes2d() .lines(&x, &y, &[Caption("A line"), Color("black")]);
    fg.show();
}
