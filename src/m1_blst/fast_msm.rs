use std::marker::PhantomData;
use ark_ff::BigInt;
use ark_serialize::CanonicalSerialize;
use blst::{blst_fp, blst_fp2, blst_p1, blst_p2, p1_affines, p2_affines};

use crate::Error;

fn convert_g1(p: &ark_bls12_381::G1Projective) -> blst_p1 {
    let x = blst_fp { l: p.x.0 .0 };
    let y = blst_fp { l: p.y.0 .0 };
    let z = blst_fp { l: p.z.0 .0 };
    blst_p1 { x, y, z }
}

fn convert_g2(p: &ark_bls12_381::G2Projective) -> blst_p2 {
    let x = blst_fp2 {
        fp: [blst_fp { l: p.x.c0.0 .0 }, blst_fp { l: p.x.c1.0 .0 }],
    };
    let y = blst_fp2 {
        fp: [blst_fp { l: p.y.c0.0 .0 }, blst_fp { l: p.y.c1.0 .0 }],
    };
    let z = blst_fp2 {
        fp: [blst_fp { l: p.z.c0.0 .0 }, blst_fp { l: p.z.c1.0 .0 }],
    };
    blst_p2 { x, y, z }
}

fn convert_g1_slice(points: &[ark_bls12_381::G1Projective]) -> Vec<blst_p1> {
    points.iter().map(convert_g1).collect()
}

fn convert_g2_slice(points: &[ark_bls12_381::G2Projective]) -> Vec<blst_p2> {
    points.iter().map(convert_g2).collect()
}

pub(crate) fn prep_g1s(points: &[ark_bls12_381::G1Projective]) -> p1_affines {
    p1_affines::from(&convert_g1_slice(points))
}

pub(crate) fn prep_g2s(points: &[ark_bls12_381::G2Projective]) -> p2_affines {
    p2_affines::from(&convert_g2_slice(points))
}

pub(crate) fn prep_scalars(scalars: &[ark_bls12_381::Fr]) -> Vec<u8> {
    let mut res = vec![0u8; 32 * scalars.len()];
    for (i, s) in scalars.iter().enumerate() {
        // This _must_ be little endian bytes for this to work
        s.serialize_compressed(&mut res[i * 32..(i + 1) * 32])
            .unwrap();
    }
    res
}

pub(crate) fn g1_msm(
    g1s: &p1_affines,
    scalars: &[u8],
    g1s_len: usize,
) -> Result<ark_bls12_381::G1Projective, Error> {
    if g1s_len < scalars.len()/32 {
        return Err(Error::PolynomialTooLarge {
            n_coeffs: scalars.len()/32,
            expected_max: g1s_len,
        });
    }
    let res_p1 = g1s.mult(scalars, 255);
    Ok(ark_bls12_381::G1Projective {
        x: ark_ff::Fp(BigInt(res_p1.x.l), PhantomData),
        y: ark_ff::Fp(BigInt(res_p1.y.l), PhantomData),
        z: ark_ff::Fp(BigInt(res_p1.z.l), PhantomData),
    })
}

pub(crate) fn g2_msm(
    g2s: &p2_affines,
    scalars: &[u8],
    g2s_len: usize,
) -> Result<ark_bls12_381::G2Projective, Error> {
    if g2s_len < scalars.len()/32 {
        return Err(Error::PolynomialTooLarge {
            n_coeffs: scalars.len()/32,
            expected_max: g2s_len,
        });
    }
    let res_p2 = g2s.mult(scalars, 255);
    Ok(ark_bls12_381::G2Projective {
        x: ark_ff::QuadExtField {
            c0: ark_ff::Fp(BigInt(res_p2.x.fp[0].l), PhantomData),
            c1: ark_ff::Fp(BigInt(res_p2.x.fp[1].l), PhantomData),
        },
        y: ark_ff::QuadExtField {
            c0: ark_ff::Fp(BigInt(res_p2.y.fp[0].l), PhantomData),
            c1: ark_ff::Fp(BigInt(res_p2.y.fp[1].l), PhantomData),
        },
        z: ark_ff::QuadExtField {
            c0: ark_ff::Fp(BigInt(res_p2.z.fp[0].l), PhantomData),
            c1: ark_ff::Fp(BigInt(res_p2.z.fp[1].l), PhantomData),
        },
    })
}

#[cfg(test)]
mod tests {
    use ark_ec::CurveGroup;
    use ark_ff::UniformRand;
    use rand::thread_rng;

    use crate::curve_msm;

    use super::*;

    #[test]
    fn test_msm_works() {
        let g1s = (0..512)
            .map(|_| ark_bls12_381::G1Projective::rand(&mut thread_rng()))
            .collect::<Vec<_>>();
        let g2s = (0..512)
            .map(|_| ark_bls12_381::G2Projective::rand(&mut thread_rng()))
            .collect::<Vec<_>>();
        let scalars = (0..512)
            .map(|_| ark_bls12_381::Fr::rand(&mut thread_rng()))
            .collect::<Vec<_>>();

        let pg1 = prep_g1s(&g1s);
        let pg2 = prep_g2s(&g2s);
        let pfr = prep_scalars(&scalars);

        let res1 = g1_msm(&pg1, &pfr, g1s.len()).unwrap();
        let res2 = g2_msm(&pg2, &pfr, g2s.len()).unwrap();

        let g1s_affine = g1s.iter().map(|p| p.into_affine()).collect::<Vec<_>>();
        let g2s_affine = g2s.iter().map(|p| p.into_affine()).collect::<Vec<_>>();

        let alt_res1 = curve_msm::<ark_bls12_381::G1Projective>(&g1s_affine, &scalars).unwrap();
        let alt_res2 = curve_msm::<ark_bls12_381::G2Projective>(&g2s_affine, &scalars).unwrap();

        assert_eq!(res1, alt_res1);
        assert_eq!(res2, alt_res2);
    }
}
