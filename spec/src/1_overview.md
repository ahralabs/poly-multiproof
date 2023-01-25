# Overview

The Polynomial Multiproof (PMP) scheme is a polynomial commitment scheme that allows for efficiently creating/verifying opening proofs for multiple polynomials at multiple points. 
Notably, opening speed is mostly _not dependent on the number of polynomials/points_ opened to.
Verification scales with the number of points/polynomials, but remains fast.

This specification outlines implementation requirements for the polynomial multiproof (PMP) scheme.
The  of information between two parties, a prover and a verifier, which allows a prover to commit to data, and efficiently generate small proofs of sections of that data for the verifier to check. 
This builds on previous methods such as KZG10[^1], and is heavily inspired by BDFG21[^2]. 
It should be seen as choosing a special case of BDFG21 which allows for significant optimizations that make the protocol more viable for use in applications like Data Availability.

The protocol depends on the selection of a pairing based curve. 
Attributes of an ideal curve for this protocol are, in order of importance,
1. Well-audited, acceptable security
2. Fast G1 scalar multiplication
3. Large scalar field size
4. Small compressed G1 scalar size

BLS12-381 satisfies the first two critera well, and while a larger scalar field would be nice and faster G1 scalar multiplication is possible, other existing curves are not as well audited. Hence BLS12-381 is a reasonable choice to start with.

Data must be fit into a curve's scalars, hence having a larger scalar allows for more data. 
The easiest way to fit data into a curve's scalar is to simply cut off one byte, such that the data fits within the modulus.

The scheme given in BDFG21 is interactive, so Merlin transcripts are used to do the Fiat-Shamir transform.

Here we present two methods, simply referred to as `M1` and `M2`.
The main difference is that `M1` is faster for opening and `M2` is faster for verifying.

Exactly why and how is outlined in [Optimizations](./a_optimizations.md)

1. `Setup` which sets up a structured reference string for the curve and does some 
2. `Commit` which can commit to an interpolated polynomial
3. `M1Open/M2Open` which computes a single opening proof that a set of polynomials are equal to specific values at a set of points
4. `M1Verify/M2Verify` which verifies an opening proof against commitments and evaluations

[^1]: [KZG10 Paper](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf)

[^2]: [BDFG21 Paper](https://eprint.iacr.org/2020/081.pdf)
