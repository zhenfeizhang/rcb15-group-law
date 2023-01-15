use halo2curves::group::ff::Field;
use halo2curves::group::ff::PrimeField;
use halo2curves::CurveExt;

pub fn add<C: CurveExt>(p1: &C, p2: &C) -> C {
    let b3 = C::Base::from(3) * C::b();
    let (x1, y1, z1) = p1.jacobian_coordinates();
    let (x2, y2, z2) = p2.jacobian_coordinates();
    let (x3, y3, z3) = core_add::<C>(x1, y1, z1, x2, y2, z2, b3);

    // we need to path halo2curve in order to skip this step
    let (x3, y3, z3) = homogeneous_to_jacobian::<C>(x3, y3, z3);
    C::new_jacobian(x3, y3, z3).unwrap()
}

pub fn double<C: CurveExt>(p: &C) -> C {
    let b3 = C::Base::from(3) * C::b();
    let (x1, y1, z1) = p.jacobian_coordinates();
    let (x3, y3, z3) = core_double::<C>(x1, y1, z1, b3);

    // we need to path halo2curve in order to skip this step
    let (x3, y3, z3) = homogeneous_to_jacobian::<C>(x3, y3, z3);
    C::new_jacobian(x3, y3, z3).unwrap()
}

pub fn mul<C: CurveExt>(base: &C, scalar: &C::ScalarExt) -> C {
    let mut res = None;
    for b in scalar
        .to_repr()
        .as_ref()
        .iter()
        .rev()
        .flat_map(|byte| (0..8).rev().map(move |i| (byte >> i) & 1u8))
    {
        if res.is_some() {
            res = Some(double(&res.unwrap()));
        }
        if b == 1 {
            if res.is_some() {
                res = Some(add(&res.unwrap(), base));
            } else {
                res = Some(*base);
            }
        }
    }
    res.unwrap()
}

#[inline]
pub fn homogeneous_to_jacobian<C: CurveExt>(
    x: C::Base,
    y: C::Base,
    z: C::Base,
) -> (C::Base, C::Base, C::Base) {
    let z = z.invert().unwrap();
    (x * z, y * z, C::Base::one())
}

pub fn homogeneous_form_to_affine<C: CurveExt>(p: &C) -> C::Affine {
    let (x, y, z) = p.jacobian_coordinates();
    let z = z.invert().unwrap();
    C::new_jacobian(x * z, y * z, C::Base::one())
        .unwrap()
        .to_affine()
}

pub fn naive_msm<C: CurveExt>(points: &[C], scalars: &[C::ScalarExt]) -> C {
    let mut res = mul(&points[0], &scalars[0]);
    for (p, s) in points.iter().zip(scalars.iter()).skip(1) {
        let tmp = mul(p, s);
        res = add(&res, &tmp)
    }
    res
}

fn core_add<C: CurveExt>(
    x1: C::Base,
    y1: C::Base,
    z1: C::Base,
    x2: C::Base,
    y2: C::Base,
    z2: C::Base,
    b3: C::Base,
) -> (C::Base, C::Base, C::Base) {
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

    // {
    //     println!("x, y, z in add {:?} {:?} {:?}", x3, y3, z3);
    //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - Fq::from(3u64) *z3 * z3 *z3;
    //     println!("on curve: {:?}", on_curve);
    // }

    (x3_27, y3_30, z3_33)
}

fn core_double<C: CurveExt>(
    x: C::Base,
    y: C::Base,
    z: C::Base,
    b3: C::Base,
) -> (C::Base, C::Base, C::Base) {
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

    // {
    //     println!("x, y, z in double {:?} {:?} {:?}", x3, y3, z3);
    //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - Fq::from(3u64) *z3 * z3 *z3;
    //     println!("on curve: {:?}\n", on_curve);
    // }
    (x3_18, y3_15, z3_10)
}
