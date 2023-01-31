//! An example usage of PMP for a data availability grid
//! This packs bytes into scalars by chunking into groups of 31 bytes and 0-padding.
//! Then it puts them into a grid sized 256x256

use ark_bls12_381::Fr;
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_ff::{PrimeField, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{end_timer, start_timer};
use merlin::Transcript;
use poly_multiproof::{
    method1::precompute::M1Precomp,
    traits::{Committer, PolyMultiProof},
    Commitment,
};
use rand::{thread_rng, RngCore};
use rayon::prelude::*;

//******************************
//* Play with these constants! *
//******************************

// The width of the grid
const GRID_WIDTH: usize = 256;
// The height of the grid (before erasure encoding)
const GRID_HEIGHT: usize = 4096;
// The number of pieces to break the grid into horizontally.
// The smaller this is, the more time PMP setup will take, but the faster opening will be.
const N_CHUNKS_W: usize = 4;
// The number of pieces to break the grid into vertically
// The bigger this is, the faster verification will be
const N_CHUNKS_H: usize = 1024;

// Can leave these alone
const CHUNK_W: usize = GRID_WIDTH / N_CHUNKS_W;
const CHUNK_H: usize = GRID_HEIGHT / N_CHUNKS_H;

struct Grid<E: Pairing> {
    pub evals: Vec<Vec<E::ScalarField>>,
    pub polys: Vec<Vec<E::ScalarField>>,
    pub commits: Vec<Commitment<E>>,
}

impl<E: Pairing> Grid<E> {
    fn from_data(data: Vec<u8>, c: &(impl Committer<E> + Sync)) -> Self {
        let pt_size = E::ScalarField::zero().serialized_size(Compress::Yes) - 1;
        let points: Vec<_> = data
            .chunks(pt_size)
            .map(|chunk| E::ScalarField::from_be_bytes_mod_order(chunk))
            .collect();

        let mut rows: Vec<_> = points
            .chunks(GRID_WIDTH)
            .map(|row| {
                let mut row: Vec<_> = row.to_vec();
                if row.len() != GRID_WIDTH {
                    println!("Padding end of row with {} elems to fit grid", row.len());
                    row.resize(GRID_WIDTH, E::ScalarField::zero());
                }
                row
            })
            .collect();

        let domain_h = GeneralEvaluationDomain::<E::ScalarField>::new(rows.len()).unwrap();
        if domain_h.size() != rows.len() {
            println!(
                "Padding {} rows to fit the best evaluation domain of size {}",
                rows.len(),
                domain_h.size()
            );
            rows.resize(domain_h.size(), vec![E::ScalarField::zero(); GRID_WIDTH]);
        }
        let domain_2h =
            GeneralEvaluationDomain::<E::ScalarField>::new(2 * domain_h.size()).unwrap();
        assert_eq!(domain_h.size(), rows.len());

        let mut interp_rows = vec![vec![E::ScalarField::zero(); GRID_WIDTH]; 2 * rows.len()];

        let erasure_t = start_timer!(|| "erasure encoding columns");
        for j in 0..GRID_WIDTH {
            let mut col = Vec::with_capacity(rows.len());
            for i in 0..rows.len() {
                col.push(rows[i][j]);
            }
            domain_h.ifft_in_place(&mut col);
            domain_2h.fft_in_place(&mut col);
            for i in 0..col.len() {
                interp_rows[i][j] = col[i];
            }
        }
        end_timer!(erasure_t);

        let domain_w = GeneralEvaluationDomain::<E::ScalarField>::new(GRID_WIDTH).unwrap();

        let poly_t = start_timer!(|| "computing polynomials from evals");
        let polys: Vec<_> = interp_rows.iter().map(|row| domain_w.ifft(&row)).collect();
        end_timer!(poly_t);

        let commit_t = start_timer!(|| "computing commitments");
        let mut commits: Vec<_> = polys
            .par_iter()
            .step_by(2)
            .map(|row| c.commit(row).expect("Commit failed").0.into_group())
            .collect();

        assert_eq!(commits.len(), domain_h.size());
        domain_h.ifft_in_place(&mut commits);
        domain_2h.fft_in_place(&mut commits);
        end_timer!(commit_t);

        Self {
            commits: commits
                .iter()
                .map(|c| Commitment(c.into_affine()))
                .collect(),
            evals: interp_rows,
            polys,
        }
    }
}

fn main() {
    let data_len = 31 * GRID_HEIGHT * GRID_WIDTH;
    let mut data = vec![0; data_len];
    rand::thread_rng().fill_bytes(&mut data);
    let domain = Radix2EvaluationDomain::<Fr>::new(GRID_WIDTH)
        .expect("Failed to make grid width eval domain");
    // Open in chunks of width/4
    let mut point_sets = vec![Vec::new(); N_CHUNKS_W];
    for (i, p) in domain.elements().enumerate() {
        let point_set = i / CHUNK_W;
        point_sets[point_set].push(p);
    }

    let pmp_t = start_timer!(|| "create pmp");
    let pmp = M1Precomp::<ark_bls12_381::Bls12_381>::new(
        GRID_WIDTH,
        point_sets.clone(),
        &mut thread_rng(),
    )
    .expect("Failed to make pmp");
    end_timer!(pmp_t);

    let grid_t = start_timer!(|| "create grid");
    let grid = Grid::from_data(data, &pmp);
    end_timer!(grid_t);

    let coords: Vec<_> = (0..N_CHUNKS_H)
        .flat_map(|i| (0..N_CHUNKS_W).map(move |j| (i, j)))
        .collect();

    let open_t = start_timer!(|| "opening to grid");
    let opens: Vec<_> = coords
        .par_iter()
        .map(|(i, j)| {
            // j is the point set index
            let start_i = i * CHUNK_H;
            let start_j = j * CHUNK_W;
            let end_i = start_i + CHUNK_H;
            let end_j = start_j + CHUNK_W;

            let polys = &grid.polys[start_i..end_i];
            let evals = &grid.evals[start_i..end_i]
                .iter()
                .map(|row| &row[start_j..end_j])
                .collect::<Vec<_>>();

            let open = pmp
                .open(&mut Transcript::new(b"example open"), evals, polys, *j)
                .expect("Failed to open");
            (*i, *j, open)
        })
        .collect();
    end_timer!(open_t);

    let veri_t = start_timer!(|| "verifying grid");
    opens.par_iter().for_each(|(i, j, proof)| {
        let start_i = *i * CHUNK_H;
        let start_j = *j * CHUNK_W;
        let end_i = start_i + CHUNK_H;
        let end_j = start_j + CHUNK_W;

        let evals = &grid.evals[start_i..end_i]
            .iter()
            .map(|row| &row[start_j..end_j])
            .collect::<Vec<_>>();
        let res = pmp
            .verify(
                &mut Transcript::new(b"example open"),
                &grid.commits[start_i..end_i],
                *j,
                evals,
                proof,
            )
            .expect(format!("Verify errored at {:>3}, {:>3}", i, j).as_str());
        if !res {
            println!("Verify failed at {:>3}, {:>3}", i, j);
        }
    });
    end_timer!(veri_t);
}
