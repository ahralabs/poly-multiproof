use ark_poly::univariate::DensePolynomial;
use ark_std::UniformRand;
use merlin::Transcript;
use std::usize;

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::rand::RngCore;

use super::{
    gen_curve_powers, gen_powers, linear_combination, poly_div_q_r, vanishing_polynomial, Error,
};
use crate::lagrange::LagrangeInterpContext;
use crate::{get_challenge, get_field_size, transcribe_points_and_evals, Commitment};

#[derive(Clone, Debug)]
pub struct Setup<E: Pairing> {
    pub powers_of_g1: Vec<E::G1Affine>,
    pub powers_of_g2: Vec<E::G2Affine>,
    point_sets: Vec<Vec<E::ScalarField>>,
    vanishing_polys: Vec<DensePolynomial<E::ScalarField>>,
    lagrange_ctxs: Vec<LagrangeInterpContext<E::ScalarField>>,
}

#[derive(Debug)]
pub struct Proof<E: Pairing>(E::G1Affine);

impl<E: Pairing> Setup<E> {
    pub fn new(
        max_degree: usize,
        max_pts: usize,
        point_sets: Vec<Vec<E::ScalarField>>,
        rng: &mut impl RngCore,
    ) -> Result<Setup<E>, Error> {
        let num_scalars = max_degree + 1;

        let x = E::ScalarField::rand(rng);
        let x_powers = gen_powers(x, num_scalars);

        let powers_of_g1 = gen_curve_powers::<E::G1>(x_powers.as_ref(), rng);
        let powers_of_g2 = gen_curve_powers::<E::G2>(x_powers[..max_pts + 1].as_ref(), rng);
        Self::new_with_powers(powers_of_g1, powers_of_g2, point_sets)
    }
    pub fn new_with_powers(
        powers_of_g1: Vec<E::G1Affine>,
        powers_of_g2: Vec<E::G2Affine>,
        point_sets: Vec<Vec<E::ScalarField>>,
    ) -> Result<Setup<E>, Error> {

        let vanishing_polys = point_sets
            .iter()
            .map(|ps| vanishing_polynomial(ps))
            .collect();
        let lagrange_ctxs = point_sets
            .iter()
            .map(|ps| LagrangeInterpContext::new_from_points(ps.as_ref()))
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(Setup {
            powers_of_g1,
            powers_of_g2,
            point_sets,
            vanishing_polys,
            lagrange_ctxs,
        })
    }

    pub fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error> {
        let res = crate::curve_msm::<E::G1>(&self.powers_of_g1, poly.as_ref())?;
        Ok(Commitment(res.into_affine()))
    }

    pub fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Proof<E>, Error> {
        // Commit the evals and the points to the transcript
        let field_size_bytes = get_field_size::<E::ScalarField>();
        let points = &self.point_sets[point_set_index];
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;

        // Read the challenge
        let gamma = get_challenge::<E::ScalarField>(transcript, b"open gamma", field_size_bytes);
        // Make the gamma powers
        let gammas = gen_powers::<E::ScalarField>(gamma, self.powers_of_g1.len());
        // Take a linear combo of gammas with the polynomials
        let fsum = linear_combination::<E::ScalarField>(polys, &gammas)
            .ok_or(Error::NoPolynomialsGiven)?;

        // Compute Z_s
        let z_s = &self.vanishing_polys[point_set_index];
        // Polynomial divide, the remained would contain the gamma * ri_s,
        // The result is the correct quotient
        let (q, _) = poly_div_q_r(DensePolynomial { coeffs: fsum }.into(), z_s.into())?;
        // Open to the resulting polynomial
        Ok(Proof(
            crate::curve_msm::<E::G1>(&self.powers_of_g1, &q)?.into_affine(),
        ))
    }

    pub fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let zeros = &self.vanishing_polys[point_set_index];
        let zeros = crate::curve_msm::<E::G2>(&self.powers_of_g2, zeros.coeffs.as_slice())?;

        let field_size_bytes = get_field_size::<E::ScalarField>();
        let points = &self.point_sets[point_set_index];
        transcribe_points_and_evals(transcript, points, evals, field_size_bytes)?;
        let gamma = get_challenge(transcript, b"open gamma", field_size_bytes);
        // Aggregate the r_is and then do a single msm of just the ri's and gammas
        let gammas = gen_powers(gamma, evals.len());

        // Get the gamma^i r_i polynomials with lagrange interp. This does both the lagrange interp
        // and the gamma mul in one step so we can just lagrange interp once.
        let gamma_ris = self.lagrange_ctxs[point_set_index]
            .lagrange_interp_linear_combo(evals, &gammas)?
            .coeffs;
        let gamma_ris_pt = crate::curve_msm::<E::G1>(&self.powers_of_g1, gamma_ris.as_ref())?;

        // Then do a single msm of the gammas and commitments
        let cms = commits.iter().map(|i| i.0).collect::<Vec<_>>();
        let gamma_cm_pt = crate::curve_msm::<E::G1>(&cms, gammas.as_ref())?;

        let g2 = self.powers_of_g2[0];

        Ok(E::pairing(gamma_cm_pt - gamma_ris_pt, g2) == E::pairing(proof.0, zeros))
    }
}

#[cfg(test)]
mod tests {
    use super::Setup;
    use crate::test_rng;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use merlin::Transcript;

    #[test]
    fn test_basic_open_works() {
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let s = Setup::<Bls12_381>::new(256, 32, vec![points.clone()], &mut test_rng())
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
