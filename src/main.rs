extern crate num_bigint;
extern crate random_fast_rng;
mod fastmod;

use num_bigint::BigUint;
use std::thread;
use std::sync::Arc;
use std::str::FromStr;
use random_fast_rng::{Random,local_rng};
use fastmod::*;

//checks against the first 100 primes
//returns true if none of the first 100 primes divide this number
fn first_pass_primality_check(n: &BigUint, prime_list: &Vec<(u32,Vec<u32>)>) -> bool {
    //don't bother checking divisibility by 2, we just don't generate even numbers at all
    for (p,pre) in prime_list.iter() {
        if fastmod_bignum_u32(n,*p,pre) == 0 {
            return false;
        }
    }
    return true;
}

//Miller-Rabin primality check
//use 48 random bases to get a > 99.999999% probability of true primeness
fn miller_rabin (n: &BigUint) -> bool {
    //first factor n - 1 to q*2^k
    let mut k = 0;
    let negative_one: BigUint = n - 1u32;
    for d in negative_one.to_u32_digits() {
        k += d.trailing_zeros();
        if d != 0 {
            break;
        }
    }
    let q: BigUint = &negative_one >> (k as usize);
    //generate our precomputation for modulo by n
    //setting an appropriate limit is important: we need to make sure that we have a big
    //enough vector to deal with the square of a number which is < n. In other words the limit
    //has to be 2*length(n)
    let limit = n.bits()/16;
    let precomputed = precompute_bignum_bignum(limit,n);
    //start doing miller-rabin with 10 pseudo random numbers
    //start with the base 2, which is a relatively good base
    let a = 2u32;
    let mut cur_power = BigUint::from(a).modpow(&q,n); //use built-in modpow for this step; it will be much faster for a nontrivial power
    if cur_power != BigUint::from(1u32) {
        for _ in 0..k {
            cur_power = fastmod_bignum_bignum(&(&cur_power*&cur_power),n,&precomputed);
            if cur_power == negative_one {
                break; //2^(2^l)q = -1 for some l, so this base is not a witness.
            }
        }
        //if we got here, then 2 is a witness, so n is definitely not prime
        return false;
    }
    //try again with base 3, another relatively good base
    let a = 3u32;
    let mut cur_power = BigUint::from(a).modpow(&q,n); //use built-in modpow for this step; it will be much faster for a nontrivial power
    if cur_power != BigUint::from(1u32) {
        for _ in 0..k {
            cur_power = fastmod_bignum_bignum(&(&cur_power*&cur_power),n,&precomputed);
            if cur_power == negative_one {
                break; //3^(2^l)q = -1 for some l, so this base is not a witness.
            }
        }
        //if we got here, then 3 is a witness, so n is definitely not prime
        return false;
    }
    //now switch to using "small" random numbers
    'try_a: for _ in 2..50 {
        let a = local_rng().get_u16() | 2; //set a bit to guarantee that a > 1. Just use a u16 for this, since small bases are usually just as good as big for Miller-Rabin
        let mut cur_power = BigUint::from(a).modpow(&q,n); //use built-in modpow for this step; it will be much faster for a nontrivial power
        if cur_power == BigUint::from(1u32) {
            continue 'try_a; //a^q = 1, so this base is not a witness
        }
        for _ in 0..k {
            cur_power = fastmod_bignum_bignum(&(&cur_power*&cur_power),n,&precomputed);
            if cur_power == negative_one {
                continue 'try_a; //a^(2^l)q = -1 for some l, so this base is not a witness.
            }
        }
        //if we got here, then a is a witness, so n is not prime
        return false;
    }
    //we were unable to find a witness, so n is probably prime
    return true;
}

//generate several "candidate" images from the supplied bitmap. The seperate utilities codify.rs
//and aks.rs can be used to convert a .bmp image to argument and check for true primality,
//rspectively
fn main() {
    //args
    let a: Vec<String> = std::env::args().collect();
    let val_inner: BigUint = BigUint::from_bytes_be(&a[1].as_bytes()) << 24;
    let num_threads = usize::from_str(a[2].as_str()).unwrap();
    let limit = val_inner.bits()/32 + 1;

    println!("performing search through {} bits with {} threads", val_inner.bits(), num_threads);
    
    //first: use the seive of erasthones to generate the first 100 primes and their
    //precomputations
    let mut primes_inner: Vec<(u32,Vec<u32>)> = Vec::new();
    let mut seive: Vec<u32> = (3..104732).collect();
    for i in 0..104729 {
        let p = seive[i];
        if p != 0 {
            let mut this_mult: usize = (p*2 - 3) as usize;
            while this_mult < 104728 {
                seive[this_mult] = 0;
                this_mult += p as usize;
            }
            primes_inner.push((p,precompute_bignum_u32(limit,p)));
        }
    }

    //put our thread variables in an Arc
    let val_global = Arc::new(val_inner);
    let primes_global = Arc::new(primes_inner);
    //spawn worker threads
    //fancy iterator version, because I'm a douchebag who hates readable code
    let child_threads: Vec<_> = (0..num_threads).map(|_| { 
        let val = Arc::clone(&val_global);
        let primes = Arc::clone(&primes_global);
        thread::spawn( move || {
            let mut try_count = 0;
            let mut rng = local_rng();
            loop {
                let mut random_part_vec: [u8;4] = [0,0,0,0];
                rng.fill_bytes(&mut random_part_vec);
                random_part_vec[0] = 0;
                random_part_vec[3] |= 1; //set the low bit so the number is never even
                for i in 1..4 {
                    random_part_vec[i] &= 0x7f; // make it valid ASCII
                    if random_part_vec[i] < 0x20 { //get rid of non-printable characters
                        random_part_vec[i] |= 0x20;
                    } else if random_part_vec[i] == 0x7f { 
                        random_part_vec[i] = 0x7d; //make sure the last bit stays set if it was set before
                    }
                }
                let fixed_random_part = u32::from_be_bytes(random_part_vec);
                let this_guess = (&val as &BigUint) + fixed_random_part; 
                try_count += 1;
                if first_pass_primality_check(&this_guess,&primes) {
                    if miller_rabin(&this_guess) {
                        let this_guess_bytes = this_guess.to_bytes_be();
                        println!("new guess:");
                        println!("\"{}\"",std::str::from_utf8(&this_guess_bytes).unwrap());
                        for b in this_guess_bytes.iter() {
                            print!("{:0>8b} ",b);
                        }
                        println!("after {} tries",try_count);
                        std::process::exit(0);
                    }
                }
            }
        })
    }).collect();
    for c in child_threads {
        c.join();
    }
}
