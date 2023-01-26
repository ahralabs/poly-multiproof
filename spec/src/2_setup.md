# 2. Methods

Here we detail the `Setup`

### Setup
Setup defines the shared parameters for all other methods
```rust
Setup(max_coeffs: usize, point_sets: Vec<Vec<Scalar>>) -> Parameters
```

- `max_coeffs` is the maximum number of coefficiens we will commit/open to at a time. For most curve choices, this should be a power of 2.
- `point_sets` lists different sets of points where we can open a grid to. 
  - This is done so we can do ahead-of-time computation that only depends on _which_ set of points we are opening at, and not on the actual contents of the data at those indicies.
  > Ex: `point_sets = [[0, 4], [1, 2]]` means I can open any rows of the grid at `[0, 4]` or `[1, 2]`. We can only open at both `[0, 4]` at the same time, and not `0` or `4` individually.

`Setup` does the following
1. Samples a random $x \in \mathbb{F}, g_1 \in \mathbb{G}_1, g_2 \in \mathbb{G}_2$
2. Computes $(g_1, g_1^x, g_1^{x^2}, \ldots g_1^{x^\mathtt{max\_coeffs - 1}})$
2. Computes $(g_2, g_2^x, g_2^{x^2}, \ldots g_2^{x^\mathtt{max\_coeffs - 1}})$ [^1]
4. For the $i$-th set of points in `point_sets`
    - Enumerates the points at the given indicies from the FFT domain. Call them $z_1, \ldots, z_k$.
    - Constructs the unique degree $k-1$ Lagrange polynomials for each $z_j$
    $$
    l_{i,j}(X) = \frac{\prod_{l \neq j} (X - z_l)}{\prod_{l \neq j} (z_j - z_l)}
    $$
      - If the denominator is zero due to a repeated point, setup errors
    - Constructs the vanishing polynomial
    $$
    Z_i(X) = \prod_j (X - z_j)
    $$
    - Computes $\theta_i = [Z_i(x)]_2$

[^1]: For method 1, we only need as many powers of $g_2$ as points we open to. For method 2, we only need $g_2, g_2^x$.

### Commit
Commiting is exactly the same as in KZG10

```rust
Commit(p: Parameters, poly: Polynomial) -> Commitment
```

Takes the polynomial $f$ and computes $[f(x)]_1$.
