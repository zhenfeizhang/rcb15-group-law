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
    let t0 = x1 * x2; // mul #1
    let t1 = y1 * y2; // mul #2
    let t2 = z1 * z2; // mul #3

    let t3 = x1 + y1;
    let t4 = x2 + y2;
    let t3 = t3 * t4; // mul #4

    let t4 = t0 + t1;
    let t3 = t3 - t4;
    let t4 = y1 + z1;

    let x3 = y2 + z2;
    let t4 = t4 * x3; // mul #5
    let x3 = t1 + t2;

    let t4 = t4 - x3;
    let x3 = x1 + z1;
    let y3 = x2 + z2;

    let x3 = x3 * y3; // mul #6
    let y3 = t0 + t2;
    let y3 = x3 - y3;

    let x3 = t0 + t0;
    let t0 = x3 + t0;
    // b3 is a constant -- so this multiplication is cheap
    let t2 = b3 * t2;

    let z3 = t1 + t2;
    let t1 = t1 - t2;

    // b3 is a constant -- so this multiplication is cheap
    let y3 = b3 * y3;

    let x3 = t4 * y3; // mul #7
    let t2 = t3 * t1; // mul #8
    let x3 = t2 - x3;

    let y3 = y3 * t0; // mul #9
    let t1 = t1 * z3; // mul #10
    let y3 = t1 + y3;

    let t0 = t0 * t3; // mul #11
    let z3 = z3 * t4; // mul #12
    let z3 = z3 + t0;

    (x3, y3, z3)
}

fn core_double<P: SWModelParameters>(
    x: P::BaseField,
    y: P::BaseField,
    z: P::BaseField,
    b3: P::BaseField,
) -> (P::BaseField, P::BaseField, P::BaseField) {
    // Algorithm 7 of eprint:2015-1060
    // Source code from A.3
    let t0 = y * y;
    let z3 = t0 + t0;
    let z3 = z3 + z3;

    let z3 = z3 + z3;
    let t1 = y * z;
    let t2 = z * z;

    let t2 = b3 * t2;
    let x3 = t2 * z3;
    let y3 = t0 + t2;

    let z3 = t1 * z3;
    let t1 = t2 + t2;
    let t2 = t1 + t2;

    let t0 = t0 - t2;
    let y3 = t0 * y3;
    let y3 = x3 + y3;

    let t1 = x * y;
    let x3 = t0 * t1;
    let x3 = x3 + x3;

    (x3, y3, z3)
}
