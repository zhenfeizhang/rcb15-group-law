use halo2curves::bn256::Fq;
use halo2curves::bn256::Fr;
use halo2curves::bn256::G1Affine;
use halo2curves::bn256::G1;
use halo2curves::group::ff::Field;
use halo2curves::group::ff::PrimeField;
use halo2curves::CurveExt;

use crate::impl_core_arith;

impl_core_arith!(Fq);

pub fn add(p1: &G1, p2: &G1) -> G1 {
    let x1 = p1.x;
    let y1 = p1.y;
    let z1 = p1.z;
    let x2 = p2.x;
    let y2 = p2.y;
    let z2 = p2.z;

    let b3 = Fq::from(3) * G1::b();

    let (x3, y3, z3) = core_add(x1, y1, z1, x2, y2, z2, b3);

    G1 {
        x: x3,
        y: y3,
        z: z3,
    }
}

pub fn double(p: &G1) -> G1 {
    let x1 = p.x;
    let y1 = p.y;
    let z1 = p.z;

    let b3 = Fq::from(3) * G1::b();

    let (x3, y3, z3) = core_double(x1, y1, z1, b3);

    G1 {
        x: x3,
        y: y3,
        z: z3,
    }
}

pub fn mul(base: &G1, scalar: &Fr) -> G1 {
    let mut res = None;
    for b in scalar
        .to_repr()
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

pub fn homogeneous_form_to_affine(x: &G1) -> G1Affine {
    let z = x.z.invert().unwrap();

    G1Affine {
        x: x.x * z,
        y: x.y * z,
    }
}

pub fn naive_msm(points: &[G1], scalars: &[Fr]) -> G1 {
    let mut res = mul(&points[0], &scalars[0]);
    for (p, s) in points.iter().zip(scalars.iter()).skip(1) {
        let tmp = mul(p, s);
        res = add(&res, &tmp)
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use halo2curves::bn256::G1;
    use halo2curves::group::{Curve, Group};
    use halo2curves::CurveAffine;

    const REPEAT: usize = 5;
    #[test]
    fn test_add() {
        let mut rng = test_rng();
        for _ in 0..REPEAT {
            let x = G1::random(&mut rng);
            let y = G1::random(&mut rng);

            // test additions
            let z = x + y;
            let res = add(&x, &y);
            let res_affine = homogeneous_form_to_affine(&res);
            assert_eq!(res_affine.is_on_curve().unwrap_u8(), 1);
            #[cfg(debug_assertions)]
            {
                println!("halo2: {:?}", z.to_affine());
                println!("rcb15: {:?}", homogeneous_form_to_affine(&res));
            }

            assert_eq!(z.to_affine(), homogeneous_form_to_affine(&res));
        }
    }

    #[test]
    fn test_double() {
        let mut rng = test_rng();

        for _ in 0..REPEAT {
            let x = G1::random(&mut rng);
            let res = double(&x);
            // test doubling via addition formula
            let res2 = add(&x, &x);
            let res3 = x + x;
            assert_eq!(res, res2);

            #[cfg(debug_assertions)]
            {
                println!("halo2: {:?}", x.to_affine());
                println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
                println!("rcb15 add: {:?}\n", homogeneous_form_to_affine(&res2));
            }

            assert_eq!(res3.to_affine(), homogeneous_form_to_affine(&res));
            assert_eq!(res3.to_affine(), homogeneous_form_to_affine(&res2));
        }
    }

    #[test]
    fn test_msm() {
        let mut rng = test_rng();

        for i in 1..REPEAT {
            let dim = 1 << i;
            let bases: Vec<_> = (0..dim).map(|_| G1::random(&mut rng)).collect();
            let bases_affine: Vec<_> = bases.iter().map(|x| x.to_affine()).collect();
            let scalars: Vec<_> = (0..dim).map(|_| Fr::random(&mut rng)).collect();

            let res = naive_msm(&bases, &scalars);
            let mut res2 = G1::identity();
            multiexp_serial(&scalars, &bases_affine, &mut res2);

            #[cfg(debug_assertions)]
            {
                println!("arkworks: {:?}", res2.to_affine());
                println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
            }

            assert_eq!(res2.to_affine(), homogeneous_form_to_affine(&res));
        }
    }

    fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
        let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

        let c = if bases.len() < 4 {
            1
        } else if bases.len() < 32 {
            3
        } else {
            (f64::from(bases.len() as u32)).ln().ceil() as usize
        };

        fn get_at<F: PrimeField>(segment: usize, c: usize, bytes: &F::Repr) -> usize {
            let skip_bits = segment * c;
            let skip_bytes = skip_bits / 8;

            if skip_bytes >= 32 {
                return 0;
            }

            let mut v = [0; 8];
            for (v, o) in v.iter_mut().zip(bytes.as_ref()[skip_bytes..].iter()) {
                *v = *o;
            }

            let mut tmp = u64::from_le_bytes(v);
            tmp >>= skip_bits - (skip_bytes * 8);
            tmp = tmp % (1 << c);

            tmp as usize
        }

        let segments = (256 / c) + 1;

        for current_segment in (0..segments).rev() {
            for _ in 0..c {
                *acc = acc.double();
            }

            #[derive(Clone, Copy)]
            enum Bucket<C: CurveAffine> {
                None,
                Affine(C),
                Projective(C::Curve),
            }

            impl<C: CurveAffine> Bucket<C> {
                fn add_assign(&mut self, other: &C) {
                    *self = match *self {
                        Bucket::None => Bucket::Affine(*other),
                        Bucket::Affine(a) => Bucket::Projective(a + *other),
                        Bucket::Projective(mut a) => {
                            a += *other;
                            Bucket::Projective(a)
                        }
                    }
                }

                fn add(self, mut other: C::Curve) -> C::Curve {
                    match self {
                        Bucket::None => other,
                        Bucket::Affine(a) => {
                            other += a;
                            other
                        }
                        Bucket::Projective(a) => other + &a,
                    }
                }
            }

            let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; (1 << c) - 1];

            for (coeff, base) in coeffs.iter().zip(bases.iter()) {
                let coeff = get_at::<C::Scalar>(current_segment, c, coeff);
                if coeff != 0 {
                    buckets[coeff - 1].add_assign(base);
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = C::Curve::identity();
            for exp in buckets.into_iter().rev() {
                running_sum = exp.add(running_sum);
                *acc = *acc + &running_sum;
            }
        }
    }
}
