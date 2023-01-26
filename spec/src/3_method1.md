# 3. Method 1

This method biases towards a faster opening time, and requires roughly half the compute of `FastVerify`.
See [A. Optimizations](./a_optimizations) for more details.

For the following, the notation of BDFG21 is followed, namely
* $[k] = \{ 1, \ldots, k\}$
* For a set of points $(x_1, \ldots, x_k)$ and a polynomial $f_i$, $r_i$ is the unique degree $k-1$ polynomial such that $\forall j \; r_i(x_j) = f_i(x_j)$ 

### M1Open

```rust
M1Open(transcript: Transcript, 
       evals: Vec<Vec<Scalar>>, 
       polys: Vec<Polynomial>, 
       point_set_index: usize) -> Proof
```
Let the polynomials be $(f_1, \ldots, f_t)$ and the points be $(x_1, \ldots x_k)$

1. Check that `evals` has length $t$ and each element of evals has length $k$
2. Transcribes the following to the Merlin transcript
    - each of the evals row by row with the utf-8 message `open evals`
    - each of the points in order with the utf-8 message `open points`
3. Constructs $\gamma \in \mathbb{F}$ by
    - Reading `SCALAR_SIZE` bytes from `transcript` with utf-8 message `open gamma`
    - Intepreting the bytes as an integer in big endian modded by the modulus of the scalar field
4. Computing
$$
f_\mathtt{sum} = \sum_{i \in [t]} \gamma^{i-1} f_i
$$
5. Polynomial divide to get $q = f_\mathtt{sum} / Z_\mathtt{point\_set\_index}$, discarding the remainer
6. Compute $[q(x)]_1$ and serialize in compressed form.

### M1Verify

```rust
M1Verify(transcript: Transcript, 
         commits: Vec<Commitment>, 
         evals: Vec<Vec<Scalar>>,
         proof: Proof,
         point_set_index: usize) -> bool
```

Let the evals be $((y_{1, 1}, \ldots y_{1, k}), \ldots (y_{t, 1}, \ldots y_{t, k}))$, commits be $(c_1, \ldots, c_t)$, and the points be $(z_1, \ldots z_k)$.
Commit $c_i$ must be for the polynomial that evaluates to $(y_{i, 1}, \ldots, y_{i, k})$ at points $(z_1, \ldots, z_k)$.

This method
1. Checks that each element has the correct length
    - Commits has length $t$
    - Each element of evals has length $k$
1. Transcribes the points/evals the same as in the opening
2. Reads $\gamma$ same is in the opening
3. For each $j \in [k]$, computes
$$
  a_j = \sum_{i \in [t]} \gamma^{j-1} y_{i,j}
$$
$a_j$ represents the value of $\sum_{i \in [t]} \gamma^{i-1} r_i (z_j)$
4. Use the $a_j$ to interpolate $\phi = \sum_{i \in [t]} \gamma^{i-1} r_i$ using the previously computed lagrange polynomials $l_{i, j}$ as
$$
\phi = \sum_{i \in [t]} \gamma^{i-1} r_i = \sum_{j \in [k]} a_j l_{\texttt{point\_set\_index}, j}
$$
5. Compute $\alpha = [\phi(x)]_1$
6. Compute $\beta = \sum_{i \in [t]} \gamma^{i-1} c_i$
7. Return `true` if 
$$
e(\beta - \alpha, g_2) = e(\texttt{proof}, \theta_{\texttt{point\_set\_index}})
$$
