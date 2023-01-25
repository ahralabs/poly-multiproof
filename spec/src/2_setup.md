# 2. Methods

Here we detail the `Setup`

### Chunking

Each curve defines a given `CHUNK_SIZE`, which is one less than the size of the curve's scalar field to ensure each chunk will fit into the scalar.
> Ex: BLS12-381, this amount is 32 - 1 = 31 

### Setup
Setup defines the shared parameters for all other methods
```rust
Setup(grid_width: usize, point_sets: Vec<Vec<Scalar>>) -> Parameters
```

- `grid_width` is the exact number of scalars we will commit/open to at a time. For most curve choices, this should be a power of 2.
- `point_sets` lists different sets of points where we can open a grid to. 
  - This is done so we can do ahead-of-time computation that only depends on _which_ set of points we are opening at, and not on the actual contents of the data at those indicies.
  > Ex: `point_sets = [[0, 4], [1, 2]]` means I can open any rows of the grid at `[0, 4]` or `[1, 2]`. We can only open at both `[0, 4]` at the same time, and not `0` or `4` individually.

`Setup` does the following
1. Samples a random $x \in \mathbb{F}, g_1 \in \mathbb{G}_1, g_2 \in \mathbb{G}_2$
2. Computes $(g_1, g_1^x, g_1^{x^2}, \ldots g_1^{x^\mathtt{grid\_width - 1}})$
2. Computes $(g_2, g_2^x, g_2^{x^2}, \ldots g_2^{x^\mathtt{grid\_width - 1}})$
4. For the $k$-th set of points in `point_sets`
    - Enumerates the points at the given indicies from the FFT domain. Call them $z_0, \ldots, z_k$.
    - Constructs the Lagrange polynomials 
    $$
    l_{k,i}(X) = \frac{\prod_{j \neq i} (X - z_j)}{\prod_{j \neq i} (z_i - z_j)}
    $$
      - If any of the denominators are zero due to a repeated point, errors
    - Constructs the vanishing polynomial
    $$
    Z_k(X) = \prod_j (X - z_j)
    $$
    - Computes $\theta_k = [Z_k(x)]_2$


