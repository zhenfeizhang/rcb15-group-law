use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_ec::{short_weierstrass_jacobian::GroupProjective, SWModelParameters};
use ark_ff::PrimeField;

pub fn homogeneous_form_to_affine<P: SWModelParameters>(x: &GroupProjective<P>) -> GroupAffine<P> {
    GroupAffine::<P>::new(x.x / x.z, x.y / x.z, false)
}

pub fn add<P: SWModelParameters>(
    p1: &GroupProjective<P>,
    p2: &GroupProjective<P>,
) -> GroupProjective<P> {
    let b3 = P::BaseField::from(3u64) * P::COEFF_B;
    let (x3, y3, z3) = core_add::<P>(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, b3);

    GroupProjective::<P>::new(x3, y3, z3)
}

pub fn double<P: SWModelParameters>(p: &GroupProjective<P>) -> GroupProjective<P> {
    let b3 = P::BaseField::from(3u64) * P::COEFF_B;
    let (x3, y3, z3) = core_double::<P>(p.x, p.y, p.z, b3);

    GroupProjective::<P>::new(x3, y3, z3)
}

/// Naive double-then-add method for group multiplications.
pub fn mul<P: SWModelParameters>(
    base: &GroupProjective<P>,
    scalar: &P::ScalarField,
) -> GroupProjective<P> {
    let mut res = None;
    for b in ark_ff::BitIteratorBE::without_leading_zeros(scalar.into_repr()) {
        if res.is_some() {
            res = Some(double(&res.unwrap()));
        }
        if b {
            if res.is_some() {
                res = Some(add(&res.unwrap(), base));
            } else {
                res = Some(*base);
            }
        }
    }
    res.unwrap()
}

/// Naive msm that does the sum of product without any optimizations.
pub fn naive_msm<P: SWModelParameters>(
    points: &[GroupProjective<P>],
    scalars: &[P::ScalarField],
) -> GroupProjective<P> {
    let mut res = mul(&points[0], &scalars[0]);
    for (p, s) in points.iter().zip(scalars.iter()).skip(1) {
        let tmp = mul(p, s);
        res = add(&res, &tmp)
    }
    res
}

fn core_add<P: SWModelParameters>(
    x1: P::BaseField,
    y1: P::BaseField,
    z1: P::BaseField,
    x2: P::BaseField,
    y2: P::BaseField,
    z2: P::BaseField,
    b3: P::BaseField,
) -> (P::BaseField, P::BaseField, P::BaseField) {
    // Algorithm 7 of eprint:2015-1060
    // Source code from A.3
    let t0_1 = x1 * x2; // mul #1
    let t1_2 = y1 * y2; // mul #2
    let t2_3 = z1 * z2; // mul #3

    let t3_4 = x1 + y1;
    let t4_5 = x2 + y2;
    let t3_6 = t3_4 * t4_5; // mul #4

    let t4_7 = t0_1 + t1_2;
    let t3_8 = t3_6 - t4_7;
    let t4_9 = y1 + z1;

    let x3_10 = y2 + z2;
    let t4_11 = t4_9 * x3_10; // mul #5
    let x3_12 = t1_2 + t2_3;

    let t4_13 = t4_11 - x3_12;
    let x3_14 = x1 + z1;
    let y3_15 = x2 + z2;

    let x3_16 = x3_14 * y3_15; // mul #6
    let y3_17 = t0_1 + t2_3;
    let y3_18 = x3_16 - y3_17;

    let x3_19 = t0_1 + t0_1;
    let t0_20 = x3_19 + t0_1;
    // b3 is a constant -- so this multiplication is cheap
    let t2_21 = b3 * t2_3;

    let z3_22 = t1_2 + t2_21;
    let t1_23 = t1_2 - t2_21;
    // b3 is a constant -- so this multiplication is cheap
    let y3_24 = b3 * y3_18;

    let x3_25 = t4_13 * y3_24; // mul #7
    let t2_26 = t3_8 * t1_23; // mul #8
    let x3_27 = t2_26 - x3_25;

    let y3_28 = y3_24 * t0_20; // mul #9
    let t1_29 = t1_23 * z3_22; // mul #10
    let y3_30 = t1_29 + y3_28;

    let t0_31 = t0_20 * t3_8; // mul #11
    let z3_32 = z3_22 * t4_13; // mul #12
    let z3_33 = z3_32 + t0_31;

    (x3_27, y3_30, z3_33)
}

fn core_double<P: SWModelParameters>(
    x: P::BaseField,
    y: P::BaseField,
    z: P::BaseField,
    b3: P::BaseField,
) -> (P::BaseField, P::BaseField, P::BaseField) {
    // Algorithm 7 of eprint:2015-1060
    // Source code from A.3
    let t0_1 = y * y;
    let z3_2 = t0_1 + t0_1;
    let z3_3 = z3_2 + z3_2;

    let z3_4 = z3_3 + z3_3;
    let t1_5 = y * z;
    let t2_6 = z * z;

    let t2_7 = b3 * t2_6;
    let x3_8 = t2_7 * z3_4;
    let y3_9 = t0_1 + t2_7;

    let z3_10 = t1_5 * z3_4;
    let t1_11 = t2_7 + t2_7;
    let t2_12 = t1_11 + t2_7;

    let t0_13 = t0_1 - t2_12;
    let y3_14 = t0_13 * y3_9;
    let y3_15 = x3_8 + y3_14;

    let t1_16 = x * y;
    let x3_17 = t0_13 * t1_16;
    let x3_18 = x3_17 + x3_17;

    (x3_18, y3_15, z3_10)
}
