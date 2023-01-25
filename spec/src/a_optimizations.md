# Method 1
This method is the fastest method for opening, and slightly slower for verification than Method 2. As written in the paper, verification is impractically slow for any appreciable number of polynomials/points. So here, we make an assumption to make the computation reasonable: each polynomial is opened at all the same points, that is (using the notation from the paper) $S_i = S_j = T \;\forall i,j \in [k]$

### Method 1, Opening
Let $T = \{ z_1, z_2, \ldots z_t\}$. For brevity, assume sums over $i$ are done as $\sum_{i \in [k]}$

The prover has to compute, $[h(x)]_1$, which means they must compute the coefficients of $h$ and then do $d+1$  $\mathbb{G}_1$ scalar multiplications and some addition. $h$ is defined as
$$
h(X) := \sum_i \gamma^{i-1} \frac{f_i(X) - r_i(X)}{Z_{S_i}(X)} = \frac{\sum_i \gamma^{i-1} \left( f_i(X) - r_i(X) \right)}{Z_T(X)}
$$
Here, $r_i$ is the unique degree $|T|$ polynomial such that $r_i(z_j) = f_i(z_j)$. Note that polynomial fraction divides cleanly, because $f_i(z_j) - r_i(z_j) = 0$ for all $j \in [t]$,  meaning it has factors of $(X-z_j)$, and so does $Z_T(X)$.  

The $r_i$ term is just the unique minimal degree polynomial you need to subtract from $f_i$ to allow $Z_T$ to cleanly divide it.  But, this also means that $r_i$ is simply the _remainder_ of dividing $f_i$ and $Z_T$.  And the same applies for $\sum_i \gamma^{i-1} f_i$: $\sum_i \gamma^{i-1} r_i$ is just the remainder you get when you divide the sum of $f_i$s by $Z_T$!  So all they need to do in order to compute $h$ is to 
1. Compute $\sum_i \gamma^{i-1} f_i$
2. Use a polynomial division algorithm to divide by $Z_T$
3. Take the result and throw away the remainder

This makes opening very fast, because no interpolation is needed.

### Method 1, Verification

For verification, the verifier needs to do a little more legwork than the prover. Let the commitment to the $i$-th polynomial be $c_i$

First, the verifier needs to compute 
$$
Z_i := [Z_{T \setminus S_i}(x)]_2 = [Z_{\emptyset}(x)]_2 = [1]_2 = g_2
$$
That's simple enough, now they have to compute
$$
\begin{aligned}
F :&=\prod_i e(\gamma^{i-1} (c_i [r_i(x)]_1), Z_i) \\
  &=e \left(\sum\gamma^{i-1} (c_i - [r_i(x)]_1), g_2\right) \\
  &=e \left(\sum\gamma^{i-1} c_i - \left[ \sum_i \gamma^{i-1} r_i(x) \right]_1, g_2\right)
\end{aligned}
$$
Constraining the points lets them compute a single pairing instead of $k$! Now let's look at each side of pairing. 
To compute the LHS, they need to compute $\sum_i \gamma^{i-1} r_i$. Sadly, there's no nice way around this here, it must be interpolated, but they can interpolate $\sum_i \gamma^{i-1} r_i$ once rather than interpolate each $r_i$.
By splitting the sum here, they can just do $2d$ $\mathbb{G}_1$ scalar multiplications, one in $\sum\gamma^{i-1} c_i$ and one in $\left[ \sum_i \gamma^{i-1} r_i(x) \right]_1$, instead of doing $(k+1)d$ scalar multiplications naively.
Then all they're left with is the single pairing.
Next they have to compute $e(W, [Z_T(x)]_2)$. Computing $Z_T$ is pretty straightforward, and all they are left to do is the $d$ $\mathbb{G}_2$ scalar multiplications.

This leaves us with the following amount of operations (other operations impact runtime fairly little)

| Operation                | Open quantity | Verify quantity     |
| ------------------------ | ------------- | ------------------- |
| Polynomial Interpolation | 0             | 1                   |
| G1 Scalar Multiplication | d             | $2d$                |
| G2 Scalar Multipcication | 0             | $d$ (can be cached) |
| Pairing                  | 0             | 2                   | 


# Method 2
Method 2 is an equally secure way to generate blocked openings, which biases computation slightly more to the opener than the verifer. They use the same assumptions as in Method 1.

### Method2, Opening

The prover computes $h(X) = f(X)/Z_T(X)$ for
$$
f(X) := \sum_i \gamma^{i-1} Z_{T \setminus S_i}(X) (f_i(X) - r_i(X)) = \sum_i \gamma^{i-1} (f_i(X) - r_i(X))
$$
since $T \setminus S_i = \emptyset$ so $Z_{T \setminus S_i} = 1$. 
This looks very familiar from method 1! Here they just compute $Z_T$ then compute the quotient $f / Z_T$ and for now, ignore the remainder. Then they compute $W := [h(x)]_1$.

This method has additional challenge other than $\gamma$, called $z$ which must take into account $W$. After seeing $W$, the verifier sends a $z \in \mathbb{F}$. In reality this should be abstracted away using a Fiat-Shamir transform.

The prover then computes $L(X) := f_z(X) - Z_T(z)h(X)$ for
$$
\begin{align}
f_z(X) :&= \sum_i \gamma^{i-1} Z_{T\S_i}(z) (f_i(X) - r_i(z)) \\
	   &= \sum_i \gamma^{i-1} f_i(X) - \sum_i \gamma^{i-1} r_i(z)
\end{align}
$$

Here we notice that we have already computed $\sum_i \gamma^{i-1} r_i(X)$ when we computed the remainder of $f / Z_T$. All they have to do is take this remainder, evaluate it at $z$, and then subtract it from $\sum_i \gamma^{i-1} f_i(X)$. 
They then evaluate $Z_T(z)$, compute  $L(X) := f_z(X) - Z_T(z) h(X)$, then compute the polynomial $\frac{L(X)}{x-z}$ and evaluate $W' := \left[\frac{L(X)}{x-z}\right]_1$, which involves $d$ $\mathbb{G}_1$ scalar multiplications. The proof then consists of $(W, W')$. 

### Method 2, Verification

The verifier receives $W, W'$, the commitments $c_i$, the evaluations of the $f_i$ at each $z_j \in T$ and computes
$$
\begin{align}
F :&= -Z_T(z) W + \sum_i \gamma^{i-1} Z_{T\setminus S_i}(z) (c_i - [r_i(z)]_1) \\
 &= -Z_T(z) W + \sum_i \gamma^{i-1} (c_i - [r_i(z)]_1) \\
 &= -Z_T(z) W + \sum_i \gamma^{i-1} c_i - \left[\sum_i \gamma^{i-1} r_i(z)\right]_1) \\
\end{align}
$$
To do this, they do lagrange interpolation of  $\sum_i \gamma^{i-1} r_i(X)$, evaluate it at $z$, then do the single $\mathbb{G}_1$ scalar multiplication. Computing $\sum_i \gamma^{i-1} c_i$ is $d$ $\mathbb{G_1}$ scalar multiplications, and one more for computing $-Z_T(z) W$.
Next the verifier checks
$e(F, g_2) = e(W', [x-z]_2)$
Which involves 2 pairings, and 2 $\mathbb{G_2}$ scalar multiplications, yielding

| Operation                | Open quantity | Verify quantity |
| ------------------------ | ------------- | --------------- |
| Polynomial Interpolation | 0             | $1$             |
| G1 Scalar Multiplication | $2d$          | $2 + d$         |
| G2 Scalar Multipcication | 0             | $2$             |
| Pairing                  | 0             | $2$             |

Comparing the two methods, we get

| Operation                | Method 1 Open | Method 2 Open | Method 1 Verify     | Method 2 Verify |
| ------------------------ | ------------- | ------------- | ------------------- | --------------- |
| Polynomial Interpolation | 0             | 0             | 1                   | $1$             |
| G1 Scalar Multiplication | d             | $2d$          | $2d$                | $2 + d$         |
| G2 Scalar Multipcication | 0             | 0             | $d$ (can be cached) | $2$             |
| Pairing                  | 0             | 0             | 2                   | $2$             |
