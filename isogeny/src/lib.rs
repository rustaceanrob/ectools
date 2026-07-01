pub mod csidh {
    use curve::{MontgomeryCurve, MontgomeryPoint, Scalar};
    use field::{Csidh512FieldOrder, FieldElement};

    pub const PRIMES: [u64; 74] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
        191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
    ];

    type Fp = Csidh512FieldOrder;
    type Fe = FieldElement<Fp>;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, std::hash::Hash)]
    pub struct Curve {
        pub a: Fe,
    }

    impl MontgomeryCurve<Fp> for Curve {
        fn a(&self) -> Fe {
            self.a
        }
        fn b(&self) -> Fe {
            Fe::one()
        }
    }

    pub fn base_curve() -> Curve {
        Curve { a: Fe::zero() }
    }

    fn twist(e: &Curve) -> Curve {
        Curve {
            a: Fe::zero() - e.a,
        }
    }

    fn cofactor_for(ell: u64) -> Scalar<Fp> {
        let mut limbs = [0u64; 8];
        limbs[0] = 4;
        for &l in PRIMES.iter() {
            if l != ell {
                mul_u64_in_place(&mut limbs, l);
            }
        }
        Scalar::from_limbs(limbs)
    }

    fn mul_u64_in_place(limbs: &mut [u64; 8], factor: u64) {
        let mut carry: u64 = 0;
        for l in limbs.iter_mut() {
            let (product, next_carry) = l.carrying_mul_add(factor, 0, carry);
            *l = product;
            carry = next_carry;
        }
        debug_assert_eq!(carry, 0, "cofactor overflowed 512 bits");
    }

    fn random_field<R: FnMut() -> u64>(next_u64: &mut R) -> Fe {
        loop {
            let mut limbs = [0u64; 8];
            for l in limbs.iter_mut() {
                *l = next_u64();
            }
            if let Some(fe) = Fe::from_limbs_checked(limbs) {
                return fe;
            }
        }
    }

    /// CSIDH-reference codomain formula for an odd-degree ell isogeny on a
    /// Montgomery curve. With (X_i : Z_i) = [i]·kernel for i = 1..(ell-1)/2,
    /// T = Π(X_i + Z_i), Ted = Π(X_i - Z_i), and Edwards projective
    /// (d_x : d_z) = (A − 2 : A + 2), the codomain update is
    /// d_x' = (A − 2)^ell · Ted^8, d_z' = (A + 2)^ell · T^8, and
    /// (A' : C') = (2·(d_z' + d_x') : d_z' − d_x').
    fn codomain(curve: &Curve, kernel: MontgomeryPoint<Fp>, ell: u64) -> (Fe, Fe) {
        let s = ((ell - 1) / 2) as usize;
        if s == 0 {
            return (curve.a, Fe::one());
        }
        let a24_num = curve.a + Fe::two();
        let a24_den = Fe::from_u64(4);

        let mut t_prod = Fe::one();
        let mut ted_prod = Fe::one();
        let mut prev = MontgomeryPoint::infinity();
        let mut curr = kernel;

        for i in 1..=s {
            t_prod *= curr.x() + curr.z();
            ted_prod *= curr.x() - curr.z();
            if i < s {
                let next = if i == 1 {
                    x_double_proj(&kernel, a24_num, a24_den)
                } else {
                    curr.x_add(&kernel, &prev)
                };
                prev = curr;
                curr = next;
            }
        }

        let ted_2 = ted_prod * ted_prod;
        let ted_4 = ted_2 * ted_2;
        let ted_8 = ted_4 * ted_4;
        let t_2 = t_prod * t_prod;
        let t_4 = t_2 * t_2;
        let t_8 = t_4 * t_4;

        let a_plus_pow = a24_num.pow_u64(ell);
        let a_minus_pow = (curve.a - Fe::two()).pow_u64(ell);

        let d_x = a_minus_pow * ted_8;
        let d_z = a_plus_pow * t_8;

        (Fe::two() * (d_z + d_x), d_z - d_x)
    }

    fn x_double_proj(p: &MontgomeryPoint<Fp>, a24_num: Fe, a24_den: Fe) -> MontgomeryPoint<Fp> {
        let a = p.x() + p.z();
        let aa = a * a;
        let b = p.x() - p.z();
        let bb = b * b;
        let e = aa - bb;
        let x = a24_den * aa * bb;
        let z = e * (a24_den * bb + a24_num * e);
        MontgomeryPoint::from_projective_unchecked(x, z)
    }

    pub fn action<R: FnMut() -> u64>(
        start: &Curve,
        private_key: &[i8; PRIMES.len()],
        next_u64: &mut R,
    ) -> Curve {
        let mut curve = *start;
        let mut e = *private_key;

        while e.iter().any(|&x| x != 0) {
            for i in 0..PRIMES.len() {
                if e[i] == 0 {
                    continue;
                }
                let ell = PRIMES[i];
                let sign: i8 = if e[i] > 0 { 1 } else { -1 };
                let work = if sign > 0 { curve } else { twist(&curve) };
                let cofactor = cofactor_for(ell);

                let torsion = loop {
                    let x = random_field(next_u64);
                    if !work.is_on_curve(x) {
                        continue;
                    }
                    let p = MontgomeryPoint::from_affine(x);
                    let q = work.mult(cofactor, p);
                    if !q.is_infinity() {
                        break q;
                    }
                };

                let (a_new, c_new) = codomain(&work, torsion, ell);
                let new_work = Curve {
                    a: a_new.normalize_by(c_new),
                };
                curve = if sign > 0 { new_work } else { twist(&new_work) };
                e[i] -= sign;
            }
        }
        curve
    }
}

#[cfg(test)]
mod csidh_tests {
    use super::csidh::*;
    use curve::MontgomeryCurve;
    use field::{Csidh512FieldOrder, FieldElement};

    #[test]
    fn base_curve_is_e0() {
        let e0 = base_curve();
        assert_eq!(e0.a, FieldElement::<Csidh512FieldOrder>::zero());
    }

    #[test]
    fn base_curve_j_invariant_is_1728() {
        assert_eq!(
            base_curve().j_invariant(),
            FieldElement::<Csidh512FieldOrder>::from_u64(1728)
        );
    }

    #[test]
    fn action_with_zero_key_is_identity() {
        let e0 = base_curve();
        let key = [0i8; PRIMES.len()];
        let result = action(&e0, &key, &mut || panic!("no RNG needed for zero key"));
        assert_eq!(result.a, e0.a);
    }

    #[test]
    fn primes_count_is_74() {
        assert_eq!(PRIMES.len(), 74);
        assert_eq!(PRIMES[0], 3);
        assert_eq!(PRIMES[PRIMES.len() - 1], 587);
    }

    #[test]
    #[ignore = "slow — run with cargo test --release -- --ignored"]
    fn action_commutes() {
        use rand::{RngExt, SeedableRng, rngs::StdRng};

        let mut rng = StdRng::seed_from_u64(0xC51D_DEAD_BEEF_0001);
        let e0 = base_curve();
        let mut a_key = [0i8; PRIMES.len()];
        a_key[0] = 1;
        let mut b_key = [0i8; PRIMES.len()];
        b_key[1] = 1;

        let a_pub = action(&e0, &a_key, &mut || rng.random::<u64>());
        let b_pub = action(&e0, &b_key, &mut || rng.random::<u64>());
        let a_shared = action(&b_pub, &a_key, &mut || rng.random::<u64>());
        let b_shared = action(&a_pub, &b_key, &mut || rng.random::<u64>());

        assert_eq!(a_shared.j_invariant(), b_shared.j_invariant());
    }
}
