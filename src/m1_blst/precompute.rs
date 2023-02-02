use ark_bls12_381::{Bls12_381, Fr, G2Projective as G2};
use ark_poly::univariate::DensePolynomial;
use merlin::Transcript;
use std::usize;

use super::{fast_msm, vanishing_polynomial, Error, Proof};
use crate::lagrange::LagrangeInterpContext;
use crate::traits::{Committer, PolyMultiProof, PolyMultiProofNoPrecomp};
use crate::Commitment;

pub struct M1Precomp {
    pub inner: super::M1NoPrecomp,
    point_sets: Vec<Vec<Fr>>,
    vanishing_polys: Vec<DensePolynomial<Fr>>,
    g2_zeros: Vec<G2>,
    lagrange_ctxs: Vec<LagrangeInterpContext<Fr>>,
}

impl M1Precomp {
    pub fn from_inner(inner: super::M1NoPrecomp, point_sets: Vec<Vec<Fr>>) -> Result<Self, Error> {
        let vanishing_polys: Vec<_> = point_sets
            .iter()
            .map(|ps| vanishing_polynomial(ps))
            .collect();
        let g2_zeros = vanishing_polys
            .iter()
            .map(|p| fast_msm::prep_scalars(&p))
            .map(|p| fast_msm::g2_msm(&inner.prepped_g2s, &p, inner.powers_of_g2.len()))
            .collect::<Result<Vec<_>, Error>>()?;
        let lagrange_ctxs = point_sets
            .iter()
            .map(|ps| LagrangeInterpContext::new_from_points(ps.as_ref()))
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(M1Precomp {
            inner,
            point_sets,
            vanishing_polys,
            g2_zeros,
            lagrange_ctxs,
        })
    }
}

impl Committer<Bls12_381> for M1Precomp {
    fn commit(&self, poly: impl AsRef<[Fr]>) -> Result<Commitment<Bls12_381>, Error> {
        self.inner.commit(poly)
    }
}

impl PolyMultiProof<Bls12_381> for M1Precomp {
    type Proof = Proof;

    fn new(
        max_coeffs: usize,
        point_sets: Vec<Vec<Fr>>,
        r: &mut impl ark_std::rand::RngCore,
    ) -> Result<Self, Error> {
        let inner = super::M1NoPrecomp::new(
            max_coeffs,
            point_sets
                .iter()
                .map(|set| set.len())
                .max_by(Ord::cmp)
                .ok_or(Error::NoPointsGiven)?
                .into(),
            r,
        )?;
        Self::from_inner(inner, point_sets)
    }

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[Fr]>],
        polys: &[impl AsRef<[Fr]>],
        point_set_index: usize,
    ) -> Result<Self::Proof, Error> {
        self.inner.open_with_vanishing_poly(
            transcript,
            evals,
            polys,
            &self.point_sets[point_set_index],
            &self.vanishing_polys[point_set_index],
        )
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<Bls12_381>],
        point_set_index: usize,
        evals: &[impl AsRef<[Fr]>],
        proof: &Self::Proof,
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
    use super::M1Precomp;
    use crate::{
        test_rng,
        traits::{Committer, PolyMultiProof},
    };
    use ark_bls12_381::Fr;
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use merlin::Transcript;

    #[test]
    fn test_basic_open_works() {
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let s = M1Precomp::new(256, vec![points.clone()], &mut test_rng())
            .expect("Failed to construct");
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
            .map(|p| s.commit(p).expect("Commit failed"))
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
