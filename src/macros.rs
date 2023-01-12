#[macro_export]
macro_rules! impl_core_arith {
    ($base: ident) => {
        fn core_add(
            x1: $base,
            y1: $base,
            z1: $base,
            x2: $base,
            y2: $base,
            z2: $base,
            b3: $base,
        ) -> ($base, $base, $base) {
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
            //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - $base::from(3u64) *z3 * z3 *z3;
            //     println!("on curve: {:?}", on_curve);
            // }

            (x3, y3, z3)
        }

        fn core_double(x: $base, y: $base, z: $base, b3: $base) -> ($base, $base, $base) {
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
            //     let on_curve = y3 * y3 * z3 - x3 * x3 * x3 - $base::from(3u64) *z3 * z3 *z3;
            //     println!("on curve: {:?}\n", on_curve);
            // }
            (x3, y3, z3)
        }
    };
}

#[macro_export]
macro_rules! impl_naive_msm {
    ($proj: ident, $scalar: ident) => {
        pub fn naive_msm(points: &[$proj], scalars: &[$scalar]) -> $proj {
            let mut res = mul(&points[0], scalars[0].into_repr().as_ref());
            for (p, s) in points.iter().zip(scalars.iter()).skip(1) {
                let tmp = mul(p, s.into_repr().as_ref());
                res = add(&res, &tmp)
            }
            res
        }
    };
}

#[macro_export]
macro_rules! impl_ark_rcb15 {
    ($proj: ident, $base: ident, $param: ident) => {
        pub fn add(p1: &$proj, p2: &$proj) -> $proj {
            let x1 = p1.x;
            let y1 = p1.y;
            let z1 = p1.z;
            let x2 = p2.x;
            let y2 = p2.y;
            let z2 = p2.z;

            let b3 = $base::from(3) * $param::COEFF_B;

            let (x3, y3, z3) = core_add(x1, y1, z1, x2, y2, z2, b3);

            let mut res = $proj::zero();
            res.x = x3;
            res.y = y3;
            res.z = z3;

            res
        }

        pub fn double(p: &$proj) -> $proj {
            let x = p.x;
            let y = p.y;
            let z = p.z;

            let b3 = $base::from(3) * $param::COEFF_B;

            let (x3, y3, z3) = core_double(x, y, z, b3);

            let mut res = $proj::zero();
            res.x = x3;
            res.y = y3;
            res.z = z3;

            res
        }

        pub fn homogeneous_form_to_affine(x: &$proj) -> <$proj as ProjectiveCurve>::Affine {
            <$proj as ProjectiveCurve>::Affine::new(x.x / x.z, x.y / x.z, false)
        }

        pub fn mul(base: &$proj, scalar: &[u64]) -> $proj {
            let mut res = None;
            for b in ark_ff::BitIteratorBE::without_leading_zeros(scalar) {
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
    };
}

#[macro_export]
macro_rules! impl_ark_rcb15_tests {
    ($proj: ident) => {
        #[test]
        fn test_add() {
            let mut rng = test_rng();

            // random additions
            for _ in 0..REPEAT {
                let x = $proj::rand(&mut rng);
                let y = $proj::rand(&mut rng);

                // test additions
                let z = x + y;
                let res = add(&x, &y);

                let res_affine = homogeneous_form_to_affine(&res);
                assert!(res_affine.is_on_curve());
                assert!(res_affine.is_in_correct_subgroup_assuming_on_curve());

                #[cfg(debug_assertions)]
                {
                    println!("arkworks: {:?}", z.into_affine());
                    println!("rcb15: {:?}", homogeneous_form_to_affine(&res));
                }

                assert_eq!(
                    z.into_affine(),
                    homogeneous_form_to_affine(&res),
                    "random add failed"
                );
            }
        }

        #[test]
        fn test_double() {
            let mut rng = test_rng();

            for _ in 0..REPEAT {
                let x = $proj::rand(&mut rng);
                let res = double(&x);
                // test doubling via addition formula
                let res2 = add(&x, &x);
                let res3 = x + x;
                assert_eq!(res, res2);

                #[cfg(debug_assertions)]
                {
                    println!("arkworks: {:?}", x.into_affine());
                    println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
                    println!("rcb15 add: {:?}\n", homogeneous_form_to_affine(&res2));
                }

                assert_eq!(res3.into_affine(), homogeneous_form_to_affine(&res));
                assert_eq!(res3.into_affine(), homogeneous_form_to_affine(&res2));
            }
        }

        #[test]
        fn test_mul() {
            let mut rng = test_rng();

            for _ in 0..REPEAT {
                let base = $proj::rand(&mut rng);
                let scalar = <$proj as ProjectiveCurve>::ScalarField::rand(&mut rng).into_repr();

                let res = mul(&base, &scalar.as_ref());
                let res2 = base.mul(scalar);

                #[cfg(debug_assertions)]
                {
                    println!("arkworks: {:?}", res2.into_affine());
                    println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
                }

                assert_eq!(res2.into_affine(), homogeneous_form_to_affine(&res));
            }
        }

        #[test]
        fn test_msm() {
            let mut rng = test_rng();

            for i in 1..REPEAT {
                let dim = 1 << i;
                let bases: Vec<_> = (0..dim).map(|_| $proj::rand(&mut rng)).collect();
                let bases_affine: Vec<_> = bases.iter().map(|x| x.into_affine()).collect();
                let scalars: Vec<_> = (0..dim)
                    .map(|_| <$proj as ProjectiveCurve>::ScalarField::rand(&mut rng))
                    .collect();
                let scalars_repr: Vec<_> = scalars.iter().map(|x| x.into_repr()).collect();

                let res = naive_msm(&bases, &scalars);
                let res2 = VariableBaseMSM::multi_scalar_mul(&bases_affine, &scalars_repr);

                #[cfg(debug_assertions)]
                {
                    println!("arkworks: {:?}", res2.into_affine());
                    println!("rcb15 double: {:?}", homogeneous_form_to_affine(&res));
                }

                assert_eq!(res2.into_affine(), homogeneous_form_to_affine(&res));
            }
        }
    };
}
