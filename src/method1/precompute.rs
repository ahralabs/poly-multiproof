use ark_poly::univariate::DensePolynomial;
use merlin::Transcript;
use std::usize;

use ark_ec::pairing::Pairing;

use super::{vanishing_polynomial, Error, Proof};
use crate::lagrange::LagrangeInterpContext;
use crate::Commitment;

#[derive(Clone, Debug)]
pub struct Setup<E: Pairing> {
    pub inner: super::Setup<E>,
    point_sets: Vec<Vec<E::ScalarField>>,
    vanishing_polys: Vec<DensePolynomial<E::ScalarField>>,
    g2_zeros: Vec<E::G2>,
    lagrange_ctxs: Vec<LagrangeInterpContext<E::ScalarField>>,
}

impl<E: Pairing> Setup<E> {
    pub fn new(
        inner: super::Setup<E>,
        point_sets: Vec<Vec<E::ScalarField>>,
    ) -> Result<Setup<E>, Error> {
        let vanishing_polys: Vec<_> = point_sets
            .iter()
            .map(|ps| vanishing_polynomial(ps))
            .collect();
        let g2_zeros = vanishing_polys
            .iter()
            .map(|p| crate::curve_msm::<E::G2>(&inner.powers_of_g2, p))
            .collect::<Result<Vec<_>, Error>>()?;
        let lagrange_ctxs = point_sets
            .iter()
            .map(|ps| LagrangeInterpContext::new_from_points(ps.as_ref()))
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(Setup {
            inner,
            point_sets,
            vanishing_polys,
            g2_zeros,
            lagrange_ctxs,
        })
    }

    pub fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Proof<E>, Error> {
        self.inner.open_with_vanishing_poly(
            transcript,
            evals,
            polys,
            &self.point_sets[point_set_index],
            &self.vanishing_polys[point_set_index],
        )
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        self.inner.verify_with_lag_ctx_g2_zeros(
            transcript,
            commits,
            &self.point_sets[point_set_index],
            evals,
            proof,
            &self.lagrange_ctxs[point_set_index],
            &self.g2_zeros[point_set_index],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::Setup;
    use crate::{method1, test_rng};
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use merlin::Transcript;

    #[test]
    fn test_basic_open_works() {
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let inner = method1::Setup::new(256, 32, &mut test_rng());
        let s = Setup::<Bls12_381>::new(inner, vec![points.clone()]).expect("Failed to construct");
        let polys = (0..20)
            .map(|_| DensePolynomial::<Fr>::rand(50, &mut test_rng()))
            .collect::<Vec<_>>();
        let evals: Vec<Vec<_>> = polys
            .iter()
            .map(|p| points.iter().map(|x| p.evaluate(x)).collect())
            .collect();
        let coeffs = polys.iter().map(|p| p.coeffs.clone()).collect::<Vec<_>>();
        let commits = coeffs
            .iter()
            .map(|p| s.inner.commit(p).expect("Commit failed"))
            .collect::<Vec<_>>();
        let mut transcript = Transcript::new(b"testing");
        let open = s
            .open(&mut transcript, &evals, &coeffs, 0)
            .expect("Open failed");
        let mut transcript = Transcript::new(b"testing");
        assert_eq!(
            Ok(true),
            s.verify(&mut transcript, &commits, 0, &evals, &open)
        );
    }
}
