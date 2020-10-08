use num_bigint::BigUint;

mod fastmod {
    pub fn precompute_bignum_bignum(limit: usize, k: &BigUint) -> Vec<BigUint> {
        let mut last_pow_pow_32 = BigUint::from(1u32);
        let mut ret: Vec<BigUint> = Vec::with_capacity(limit);    
        for _ in 0..limit+1 {    
            ret.push(last_pow_pow_32.clone());    
            last_pow_pow_32 = (last_pow_pow_32 << 32) % k;    
        }    
        ret    
    }

    pub fn precompute_bignum_u32(limit: usize, k: u32) -> Vec<u32> {
        let mut last_pow_pow_32: u32 = 1;
        let mut ret: Vec<u32> = Vec::with_capacity(limit);
        for _ in 0..limit+1 {
            ret.push(last_pow_pow_32);
            last_pow_pow_32 = ((u64::from(last_pow_pow_32) << 32) % u64::from(k)) as u32;
        }
        ret
    }
    
    //assumes you have enough digits precomputed to do the thing here
    pub fn fastmod_bignum_bignum(n: &BigUint, k: &BigUint, precomputed: &Vec<BigUint>) -> BigUint {
        let mut sum = BigUint::from(0u32);
        let digits = n.to_u32_digits();
        for (i,d) in digits.iter().enumerate() {
            sum += &precomputed[i]*d;
        }
        sum % k
    }
    
    pub fn fastmod_bignum_u32(n: &BigUint, k: u32, precomputed: &Vec<u32>) -> u32 {
        let mut sum: u64 = 0;
        let digits = n.to_u32_digits();
        for (i,d) in digits.iter().enumerate() {
            sum += u64::from(precomputed[i])*u64::from(*d);
        }
        (sum % u64::from(k)) as u32
    }
}
