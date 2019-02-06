// Forked from Lib-SNARK.

use std::ops::{ AddAssign, MulAssign, Sub };
use num::{ One, Zero };

pub struct Field<FieldT>{ omega: FieldT }

pub trait FieldTraits: 
    One
    + Zero
    + From<usize>
    + Sub<Self, Output=Self>
    + MulAssign
    + AddAssign
    + PartialEq
    + Copy
{ 
    fn s<S>() -> S; 
    fn inverse<U>(self) ->  U;
}

// basic_radix2_domain_aux.tcc 45-79
// note that: 'it's the caller's responsibility to multiply by 1/N'.
pub fn _basic_serial_radix2_FFT<FieldT>(mut a: Vec<FieldT>, omega: FieldT) 
where
    FieldT: FieldTraits
{
    let mut n: usize = a.len();
    let mut logn: f64 = (n as f64).log2();
    
    // impl check for n != (1u << logn) 48.
    
    // 51-56
    /* swapping in place for fft */
    for k in 0..n 
    {
        let mut rk = subroutines::bit_reverse(k, logn as usize);
        if k < rk { a[k] = FieldT::from(rk) }
    }

    // invariant: m = 2^{s-1} 
    // 58-79.
    let mut m: usize = 1;
    for _ in 0..=logn as usize
    {
        // 'w_m is 2^s-th root of unity now'
        // DOUBLE CHECK THAT THIS RAISES THE POWER CORRECTLY
        let mut w_m: FieldT =  num::pow( omega, n/(2*m) );

        // asm volatile  ("/* pre-inner */"); 64
        // there is a rust macro for this.
        for k in 0..n
        {
            let mut w: FieldT = FieldT::one();
            for j in 0..m
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

// there are then some multi-core options.Copy

// basic_radix2_domain_aux.tcc 182-236
pub fn _basic_radix2_evaluate_all_lagrange_polynomials<FieldT>(m: usize, t: FieldT) -> Vec<FieldT>
where
    FieldT: FieldTraits
{
    if m==1 { return vec![ FieldT::one() ] };
    // if m != 1u check 190.

    let omega: FieldT = subroutines::get_root_of_unity::<FieldT>(m as f64);
    let mut u = vec![FieldT::zero(); m];

    // NOTE: If t equals one of the roots of unity in S={omega^{0},...,omega^{m-1}}
    // then output 1 at the right place, and 0 elsewhere.

    if num::pow(t, m) == FieldT::one()
    {
        let mut omega_i = FieldT::one();
        for i in 0..m
        {
            // "t = omega^i"
            if omega_i == t 
            { 
                u[i] = FieldT::one(); 
                return u
            } else { omega_i *= omega; }
        };
    };

    // NOTE:   Otherwise, if t does not equal any of the roots of unity in S,
    // then compute each L_{i,S}(t) as Z_{S}(t) * v_i / (t-\omega^i)
    // where:
    // - Z_{S}(t) = \prod_{j} (t-\omega^j) = (t^m-1), and
    // - v_{i} = 1 / \prod_{j \neq i} (\omega^i-\omega^j).
    // Below we use the fact that v_{0} = 1/m and v_{i+1} = \omega * v_{i}.
    
    let Z: FieldT = num::pow(t, m)-FieldT::one();
    let mut l: FieldT = Z * FieldT::from(m).inverse(); 
    let mut r: FieldT = FieldT::one();
    for i in 0..m 
    {
        u[i] = l * (t-r).inverse();
        l *= omega;
        r *= omega;
    }
    u
}

mod subroutines {
    use super::FieldTraits;

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

    // basic_radix2_domain_aux.tcc 171-180.
    pub fn _multiply_by_coset<FieldT>( a: Vec<FieldT>, g: FieldT)
    where
        FieldT: FieldTraits,
    {
        let mut u: FieldT = g;
        for mut i in a.into_iter()
        {
            i *= u;
            u *= g;
        }
    }

    // libff/algebra/fields/field_utils.tcc 38-51;
    pub fn get_root_of_unity<FieldT>(n: f64) -> FieldT
    where
        FieldT: FieldTraits
    {
        let logn = n.log2();
        // if (n!=1u<<logn) && if lohn > FieldT::s
        let omega: FieldT = FieldT::root_of_unity();
        for i in [FieldT::s::<f64>()..logn].iter().rev()
        {
            omega *= omega;
        };
        omega
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