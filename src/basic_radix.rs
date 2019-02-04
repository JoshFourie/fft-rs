// Forked from Lib-SNARK.

// basic_radix2_domain_aux.tcc 45-79
pub fn _basic_serial_radix2_FFT<FieldT>(mut a: Vec<FieldT>, omega: FieldT) 
where
    FieldT: From<usize>
    + num::One
    + std::ops::MulAssign
    + Copy,

{
    let n: usize = a.len();
    let logn: f64 = (n as f64).log2();
    
    // impl check for n != (1u << logn) 48.
    
    // 51-56
    /* swapping in place for fft */
    for k in 0..n 
    {
        let rk = subroutines::bit_reverse(k, logn as usize);
        if k < rk { a[k] = FieldT::from(rk) }
    }

    // invariant: m = 2^{s-1} 
    // 58-79.
    let m: usize = 1;
    for _ in 0..=logn as usize
    {
        // 'w_m is 2^s-th root of unity now'
        // DOUBLE CHECK THAT THIS RAISES THE POWER CORRECTLY
        let w_m: FieldT =  num::pow( omega, n/(2*m) );

        // asm volatile  ("/* pre-inner */"); 64
        // there is a rust macro for this.
        for k in 0..n
        {
            let w: FieldT = FieldT::one();
            for k in 0..m
            {
                let t: FieldT = w * a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] += t;
                w *=w_m; 
            }   
        }
        // asm volatile ("/* post-inner */"); 76
        m *= 2;   
    }
}

mod subroutines {

    // libff/libff/common.tcc 60-69
    pub fn bit_reverse(mut n: usize, l: usize) -> usize
    {
        let mut r=0;
        for _ in 0..l
        {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    #[cfg(test)]
    mod tests {

        #[test]
        // expected that the bit_reverse fn() wil return a mirrored bit repr.
        fn bit_reverse() 
        {
            assert_eq!(4, super::bit_reverse(1, 3));
            assert_eq!(5, super::bit_reverse(5, 3));
        }

    }
}