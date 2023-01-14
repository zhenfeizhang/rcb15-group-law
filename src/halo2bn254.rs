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

    // {
    //     println!("x, y, z in add {:?} {:?} {:?}", x3, y3, z3);
    //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - Fq::from(3u64) *z3 * z3 *z3;
    //     println!("on curve: {:?}", on_curve);
    // }

    (x3, y3, z3)
}

fn core_double<C: CurveExt>(
    x: C::Base,
    y: C::Base,
    z: C::Base,
    b3: C::Base,
) -> (C::Base, C::Base, C::Base) {
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

    // {
    //     println!("x, y, z in double {:?} {:?} {:?}", x3, y3, z3);
    //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - Fq::from(3u64) *z3 * z3 *z3;
    //     println!("on curve: {:?}\n", on_curve);
    // }
    (x3, y3, z3)
}
