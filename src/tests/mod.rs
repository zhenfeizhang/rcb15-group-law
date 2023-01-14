// mod ark_bls12_377;
mod ark_bls12_381;
mod ark_bn254;
mod halo2_bn254;
mod halo2_pallas;

#[macro_export]
macro_rules! impl_ark_rcb15_tests {
    ($proj: ident) => {
        #[test]
        fn test_add() {
            let mut rng = test_rng();

            // random additions
            for _ in 0..REPEAT {
                let x = GroupProjective::<$proj>::rand(&mut rng);
                let y = GroupProjective::<$proj>::rand(&mut rng);

                // test additions
                let z = x + y;
                let res = add::<$proj>(&x, &y);

                println!(
                    "{}",
                    res.y * res.y * res.z
                        - res.x * res.x * res.x
                        - <$proj as ModelParameters>::BaseField::from(3u64) * res.z * res.z * res.z
                );

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
                let x = GroupProjective::<$proj>::rand(&mut rng);
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
                let base = GroupProjective::<$proj>::rand(&mut rng);
                let scalar = <$proj as ModelParameters>::ScalarField::rand(&mut rng);

                let res = mul(&base, &scalar);
                let res2 = base.mul(scalar.into_repr());

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
                let bases: Vec<_> = (0..dim)
                    .map(|_| GroupProjective::<$proj>::rand(&mut rng))
                    .collect();
                let bases_affine: Vec<_> = bases.iter().map(|x| x.into_affine()).collect();
                let scalars: Vec<_> = (0..dim)
                    .map(|_| <$proj as ModelParameters>::ScalarField::rand(&mut rng))
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
