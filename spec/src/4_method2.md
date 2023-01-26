# 4. Method 2

This method biases towards a faster verifying time, and is roughly twice as slow to open on.
The API is exactly the same as in method 1.
See [A. Optimizations](./a_optimizations) for more details.

### M2Open

```rust
M2Open(transcript: Transcript, 
       evals: Vec<Vec<Scalar>>, 
       polys: Vec<Polynomial>, 
       point_set_index: usize) -> Proof
```
This method 
1. Transcribes the following (all points serialized compressed if possible)
    - each of the evals row by row with the utf-8 message `open evals`
    - each of the points in order with the utf-8 message `open points`
2. Constructs $\gamma \in \mathbb{F}$ by
    - Reading `SCALAR_SIZE` bytes from `transcript` with utf-8 message `open gamma`
    - Intepreting the bytes as an integer in big endian modded by the modulus of the scalar field
3. Computing, for each polynomial ${f_1, \ldots, f_t}$
$$
f_\mathtt{sum} = \sum_{i \in [t]} \gamma^{i-1} f_i
$$
4. Polynomial divide to get $h = f_\mathtt{sum} / Z_\mathtt{point\_set\_index}$, with remainder $v$.
5. Compute $W_1 = [q(x)]_1$ 
6. Serialize $W_1$ in compressed form and transcribe with the message `open W1`.
7. Construct $z$ in the same way as we did $\gamma$, reading with message `open z`.
8. Compute $s = \sum_{i \in [t]} \gamma^{i - 1} r_i(z)$ by computing $(Z_\texttt{point\_set\_index} \cdot v)(z)$. [^1]
9. Comupute $f_z = -s + \sum_{i \in [t]} f_i$
10. Compute $l = f_z - h Z_\mathtt{point\_set\_index}$
11. Compute $b = l / (X - z)$
12. Compute $W_2 = [b(x)]_1$
13. Return $(W_1, W_2)$ and serialize in compressed form

### M2Verify

```rust
M2Verify(transcript: Transcript, 
         commits: Vec<Commitment>, 
         evals: Vec<Vec<Scalar>>,
         (w_1, w_2): Proof,
         point_set_index: usize) -> bool
```

Let the evals be $((y_{1, 1}, \ldots y_{1, k}), \ldots (y_{t, 1}, \ldots y_{t, k}))$, commits be $(c_1, \ldots, c_t)$, and the points be $(z_1, \ldots z_k)$.
Commit $c_i$ must be for the polynomial that evaluates to $(y_{i, 1}, \ldots, y_{i, k})$ at points $(z_1, \ldots, z_k)$.

This method
1. Transcribes the points/evals the same as in the opening
2. Reads $\gamma$ same is in the opening
3. Transcribes $W_1$
4. Reads $z$ same is in opening
5. Lagrange interpolate $\phi = \sum_{i \in [t]} \gamma^{i-1} r_i$ using the given evaluations, same as in method 1.
6. Compute $\alpha = g_1^{\phi(z)}$
7. Compute $\beta = \sum_{i \in [t]} \gamma^{i-1} c_i$
8. Compute $F = \beta - \alpha - W_1^{Z_{\texttt{proof\_set\_index}}(z)}$
9. Compute $g_2^xg_2^{-z} = g_2^{x-z}$
10. Return `true` if 
$$
e(F, g_2) = e(W_1, g_2^{x-z})
$$

[^1]: For why this works, see the [A. Optimizations](./a_optimizations.md)
