use curve::MontgomeryCurve;
use isogeny::csidh;
use rand::{RngExt, rng};

const KEY_BOUND: usize = 74;

fn random_private_key<R: RngExt>(rng: &mut R) -> [i8; csidh::PRIMES.len()] {
    let mut key = [0i8; csidh::PRIMES.len()];
    for slot in key.iter_mut().take(KEY_BOUND) {
        *slot = (rng.random::<u8>() % 11) as i8 - 5;
    }
    key
}

fn main() {
    let mut rng = rng();
    let a_priv = random_private_key(&mut rng);
    let b_priv = random_private_key(&mut rng);
    println!("Alice's private:      {:?}", &a_priv[..KEY_BOUND]);
    println!("Bob's private:        {:?}", &b_priv[..KEY_BOUND]);

    let start = csidh::base_curve();

    println!("\nDeriving public curves...");
    let t = std::time::Instant::now();
    let a_pub = csidh::action(&start, &a_priv, &mut || rng.random::<u64>());
    let b_pub = csidh::action(&start, &b_priv, &mut || rng.random::<u64>());
    println!("done in {:?}", t.elapsed());

    println!("\nDeriving shared secrets...");
    let t = std::time::Instant::now();
    let a_shared = csidh::action(&b_pub, &a_priv, &mut || rng.random::<u64>());
    let b_shared = csidh::action(&a_pub, &b_priv, &mut || rng.random::<u64>());
    println!("done in {:?}", t.elapsed());

    let j_a = a_shared.j_invariant();
    let j_b = b_shared.j_invariant();
    assert_eq!(j_a, j_b, "shared j-invariants disagreed");
    println!("\nShared j-invariant: {:?}", j_a);
    println!("CSIDH shared secret established.");
}
