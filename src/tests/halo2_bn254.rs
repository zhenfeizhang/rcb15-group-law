use std::ops::Mul;

use ark_std::test_rng;
use halo2curves::bn256::Fr;
use halo2curves::bn256::G1;
use halo2curves::group::ff::Field;
use halo2curves::group::ff::PrimeField;
use halo2curves::group::{Curve, Group};
use halo2curves::CurveAffine;

use crate::halo2::add;
use crate::halo2::double;
use crate::halo2::homogeneous_form_to_affine;
use crate::halo2::mul;
use crate::halo2::naive_msm;

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
fn test_mul() {
    let mut rng = test_rng();

    for _ in 0..REPEAT {
        let base = G1::random(&mut rng);
        let scalar = Fr::random(&mut rng);

        let res = mul(&base, &scalar);
        let res2 = base.mul(scalar);

        #[cfg(debug_assertions)]
        {
            println!("halo2: {:?}", res2.to_affine());
            println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
        }

        assert_eq!(res2.to_affine(), homogeneous_form_to_affine(&res));
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
            println!("halo2: {:?}", res2.to_affine());
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
