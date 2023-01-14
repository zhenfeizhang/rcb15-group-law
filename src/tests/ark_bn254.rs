use crate::arkworks::add;
use crate::arkworks::double;
use crate::arkworks::homogeneous_form_to_affine;
use crate::arkworks::mul;
use crate::arkworks::naive_msm;
use crate::impl_ark_rcb15_tests;
use ark_bn254::g1::Parameters;
use ark_ec::msm::VariableBaseMSM;
use ark_ec::short_weierstrass_jacobian::GroupProjective;
use ark_ec::ModelParameters;
use ark_ec::ProjectiveCurve;
use ark_ff::PrimeField;
use ark_std::test_rng;
use ark_std::UniformRand;

const REPEAT: usize = 5;

impl_ark_rcb15_tests!(Parameters);
