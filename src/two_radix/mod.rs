use num::{ Complex, Num, Float};
use std::ops::{ Div, Mul };

pub mod leaf;
pub mod node;
pub mod tree;
/*************************************************** Structures ****************************************/

#[derive(Debug, PartialEq, Clone)] pub struct PointWise<F> { pub seq: Vec<Complex<F>> }

/*********************************************** Implementations **************************************/

impl<F> From<Vec<F>> for PointWise<F> where F: Clone + Num
{
    fn from(target: Vec<F>) -> Self 
    {
        PointWise{ seq: target.into_iter().map(|num| Complex::from(num)).collect::<Vec<Complex<F>>>() } 
    }
}
/*************************************************** Functions ****************************************/

pub fn butterfly(lhs: Complex<f64>, rhs: Complex<f64>, small_n: usize, big_n: usize) -> (Complex<f64>, Complex<f64>)
{
    let ohm=rhs*twiddle(small_n, big_n);
    (lhs + ohm, lhs - ohm)
}

pub fn twiddle(small_n: usize, big_n: usize) -> Complex<f64>
{
    let pi = std::f64::consts::PI;        
    return (
        Complex::from(
        [-2.0, pi, small_n as f64, 1.0.div( big_n as f64)]
            .iter()
            .product::<f64>()
        ) * Complex::i()
    ).exp()
}

pub fn danielson_lanczos_pattern<F>(target: &mut Vec<F>, bits: usize)
{
    let bit_reverse = |mut num: usize, bits: usize| -> usize 
    {
        let mut rev=0;
        for _ in 0..bits
        {
            rev = (rev << 1) | (num & 1);
            num >>= 1;
        } 
        rev
    };
    for i in 0..target.len().div(2) { target.swap( i, bit_reverse(i, bits) ) };
}

/*************************************************** Convenience ****************************************/

pub trait RoundTo<T> { fn round_to(self, x: T) -> Self; }

impl<T> RoundTo<T> for Complex<T> 
where
    T: Mul + Clone + Num + Float
{ 
    fn round_to(self, x: T) -> Self 
    { 
        Complex::new(
            (self.re * x).round()/x,
            (self.im * x).round()/x, 
        )
    } 
}

#[cfg(test)]
mod tests
{
    use super::*;
    use num::Complex;
    use std::ops::Div;
    use assert_approx_eq::assert_approx_eq; 

    pub fn de_moivre(element_deg: f64, iter_deg: f64, deg_bound: f64) -> Complex<f64> 
    { 
        let pi = std::f64::consts::PI;            
        return (
            Complex::from(
            [-2_f64, pi, element_deg, iter_deg, 1_f64.div( deg_bound )]
                .iter()
                .product::<f64>()
            ) * Complex::i()
        ).exp()
    
    }

    pub fn single_dft(coeffs: &Vec<Complex<f64>>, element_deg:f64, max_deg_divisor: f64) -> Complex<f64>
    {
        let mut sigma: Vec<Complex<f64>> = Vec::new();
        for n in 0..coeffs.len()
        {
            sigma.push( de_moivre( element_deg, n as f64, (coeffs.len() as f64).div(max_deg_divisor) ))
        }
        sigma.iter()
            .zip(coeffs)
            .map( |(x, ohm)| x * ohm )
            .sum()
    }

    pub fn seq_dft(coeffs: &Vec<Complex<f64>>) -> Vec<Complex<f64>> 
    {
        coeffs.iter()
            .enumerate()
            .map(|(k, _)| single_dft( coeffs, k as f64, 1_f64 ))
            .collect::<Vec<_>>()
    }

    #[test]
    fn test_twiddle()
    {
        let exp_x: Complex<f64> = Complex::i();
        let test_x: Complex<f64> = twiddle(6, 8);
        assert_approx_eq!(exp_x.re, test_x.re);
        assert_approx_eq!(exp_x.im, test_x.im);

        let exp_y: Complex<f64> = Complex::new(-1_f64.div(2.0.sqrt()), 1_f64.div( 2.0.sqrt() ));
        let test_y: Complex<f64> = twiddle(5, 8);
        assert_approx_eq!(exp_y.re, test_y.re);
        assert_approx_eq!(exp_y.im, test_y.im);

        let exp_z: Complex<f64> = Complex::new(-1_f64.div(2.0.sqrt()), -1_f64.div( 2.0.sqrt() ));
        let test_z: Complex<f64> = twiddle(3, 8);
        assert_approx_eq!(exp_z.re, test_z.re);
        assert_approx_eq!(exp_z.im, test_z.im);

        let dummy: Complex<f64> = Complex::new(-5_f64.div(2.0.sqrt()), 1_f64.div( 2.0.sqrt() ));
        let test_dummy: Complex<f64> = twiddle(3, 8);
        assert!(dummy != test_dummy)        
    }

    #[test]
    fn test_butterfly_assumption()
    {
        let test_one = |x, y|
        {
            let twiddle_x=twiddle(x,8);
            let twiddle_y=twiddle(y,8);
            assert_approx_eq!(twiddle_x.re, -twiddle_y.re);
            assert_approx_eq!(twiddle_x.im, -twiddle_y.im);
        };

        test_one(4, 0);
        test_one(7, 3);

        let test_two = |x, y|
        {
            let (test_lhs, test_rhs) = butterfly(Complex::from(2.0), Complex::from(1.0), x, 8);
            let exp_lhs = Complex::from(2.0) + Complex::from(1.0) * twiddle(x, 8);
            let exp_rhs = Complex::from(2.0) + Complex::from(1.0) * twiddle(y, 8);

            assert_approx_eq!(test_lhs.re, exp_lhs.re);
            assert_approx_eq!(test_lhs.im, exp_lhs.im);
            assert_approx_eq!(test_rhs.re, exp_rhs.re);
            assert_approx_eq!(test_rhs.im, exp_rhs.im);
        };

        test_two(1, 5);
        test_two(3, 7);
    }

    #[test]
    fn test_dl_pattern_and_bit_rev() 
    {
        let mut target_8_input = vec![0, 1, 2, 3, 4, 5, 6, 7];  
        let expected_8_input = vec![0, 4, 2, 6, 1, 5, 3, 7];
        let dummy_8_input = vec![0, 4, 2, 6, 1, 5, 3, 8];

        // check that positions not vals are swapped.
        let mut target_4_input = vec![0, 134, 275, 319];
        let expected_4_input = vec![0, 275, 134, 319];
        let dummy_4_input = vec![3, 1, 2, 0];

        danielson_lanczos_pattern(&mut target_4_input, 2); 
        danielson_lanczos_pattern(&mut target_8_input, 3);

        assert_eq!( &target_8_input , &expected_8_input );
        assert_eq!( &target_4_input , &expected_4_input );
        assert!( &target_8_input != &dummy_8_input );
        assert!( &target_4_input != &dummy_4_input );
    }
}