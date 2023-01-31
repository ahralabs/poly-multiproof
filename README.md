# poly-multiproof

Polynomial commitment schemes with fast openings on multiple points. 
Methods and runtime:

The two methods here are inspired by [BDFG21](https://eprint.iacr.org/2020/081.pdf).


### Examples

An example of using pmp for a grid data availability scheme with 1d erasure encoding is in `examples/data_availability_grid.rs`. To run it with a nice timer, do
```bash
cargo run --example data_availability_grid --release --features print-trace
```

which will print out something like this example for a 256x256 grid
```
Start:   create pmp
End:     create pmp ........................................330.880ms
Start:   create grid
··Start:   erasure encoding columns
··End:     erasure encoding columns ........................32.143ms
··Start:   computing polynomials from evals
··End:     computing polynomials from evals ................22.562ms
··Start:   computing commitments
··End:     computing commitments ...........................707.127ms
End:     create grid .......................................771.776ms
Start:   opening to grid
End:     opening to grid ...................................244.446ms
Start:   verifying grid
End:     verifying grid ....................................223.557ms
```

There are nice constants in the top of the file to play with.

### Benchmarks

To run benchmarks with asm optimizations on x86 machines, run
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm
```
or to run with the goal of plotting, run
```bash
RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm --plotting-backend disabled -- --quick --quiet &> bench_out.txt
```
The logs in `bench_out.txt` can then be parsed and plotted in `Plot Benches.ipynb`. 
Using `--quick` is nice since there are many many inputs benchmarked and it will still take an hour or so to run with `--quick`.
