use crate::impl_ark_rcb15;
use crate::impl_core_arith;
use crate::impl_naive_msm;
use ark_bn254::g1::Parameters;
use ark_bn254::Fq;
use ark_bn254::Fr;
use ark_bn254::G1Projective;
use ark_ec::{ProjectiveCurve, SWModelParameters};
use ark_ff::PrimeField;
use ark_std::Zero;

impl_core_arith!(Fq);
impl_ark_rcb15!(G1Projective, Fq, Parameters);
impl_naive_msm!(G1Projective, Fr);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impl_ark_rcb15_tests;
    use ark_ec::msm::VariableBaseMSM;
    use ark_ff::PrimeField;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    const REPEAT: usize = 5;

    impl_ark_rcb15_tests!(G1Projective);
}
