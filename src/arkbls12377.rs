use crate::impl_ark_rcb15;
use crate::impl_core_arith;
use ark_bls12_377::g1::Parameters;
use ark_bls12_377::Fq;
use ark_bls12_377::G1Projective;
use ark_ec::{ProjectiveCurve, SWModelParameters};
use ark_std::Zero;

impl_core_arith!(Fq);
impl_ark_rcb15!(G1Projective, Fq, Parameters);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impl_ark_rcb15_tests;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    const REPEAT: usize = 5;
    impl_ark_rcb15_tests!(G1Projective);
}
