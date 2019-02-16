use std::ops::{Add, Mul, Neg};
use zksnark::{CoefficientPoly, field::Field};

/************************** TRAIT INTERFACE ***********************/

pub trait PointRep<P>: Copy + Add + Neg + Mul {
    fn degree(self) -> P;
    fn y(self) -> P;  
}

pub trait PointWiseRep: Clone + Add + Neg + Mul { }

/************************** STRUCTS ***********************/

#[derive(Debug, PartialEq, Copy, Clone)] pub struct Points<P> { degree: P, y: P }

#[derive(Debug, PartialEq, Clone)] pub struct PointWise<P> { points: Vec<Points<P>> }

/************************** IMPLEMENTATIONS ***********************/

impl<P> PointRep<P> for Points<P> 
where 
    P: Add<P, Output=P> + Neg<Output=P> + Mul<P, Output=P> + Copy
{
    fn degree(self) -> P {
        self.degree
    }
    fn y(self) -> P {
        self.y
    }
}

impl<P> From<(P, P)> for Points<P> 
{ 
    fn from( (degree, y): (P, P) ) -> Self { 
        Self { degree , y } 
    } 
}

impl<P> Add<Self> for Points<P> 
where 
    P: Add<Output=P>
{ 
    type Output=Self;
    fn add(self, rhs: Self) -> Self { 
        Self::from( (self.degree, self.y + rhs.y) ) 
    } 
}

impl<P> Neg for Points<P> 
where
    P: Neg<Output=P>
{ 
    type Output=Self;
    fn neg(self) -> Self { 
        Self::from( (self.degree, -self.y) ) 
    } 
}

impl<P> Mul<Self> for Points<P> 
where
    P: Mul<Output=P>
{
    type Output=Self;
    fn mul(self, rhs: Self) -> Self { 
        Self::from( (self.degree, self.y * rhs.y) ) 
    } 
}

impl<P> From<Vec<(P, P)>> for PointWise<P> {
    fn from(object: Vec<(P, P)>) -> Self {
        Self {
            points: object
                .into_iter()
                .map( |(degree,y)| Points::from( (degree, y ) ) )
                .collect::<Vec<_>>()
        }   
    }
}

impl<P> Add<Self> for PointWise<P> 
where 
    P: Add<P, Output=P>
{
    type Output = Self;    
    fn add(self, rhs: Self) -> Self {
        Self::from(
            self.points
                .into_iter()
                .zip(rhs.points.into_iter())
                .map( |(a, b)| (a.degree, a.y + b.y) )
                .collect::<Vec<_>>()
        )
    }
}       

impl<P> Neg for PointWise<P> 
where 
    P: Neg<Output=P>
{
    type Output=Self;
    fn neg(self) -> Self {
        Self::from( 
            self.points
                .into_iter()
                .map( |i| (-i.degree, -i.y ) )
                .collect::<Vec<_>>()
        )
    }
}

impl<P> Mul<Self> for PointWise<P> 
where 
    P: Mul<P, Output=P>
    + Neg<Output=P>
    + Clone
{
    type Output=Self;
    fn mul(self, rhs: Self) -> Self {
        Self::from (
            self.points
                .into_iter()
                .zip(rhs.points.into_iter())
                .map( |(a, b)| (a.degree, a.y * b.y) )
                .collect::<Vec<_>>()
        )
    }
}

pub mod fft_subroutines {

    use std::ops::*;
    use zksnark::{*, field::Field};
    use num::{Complex, Integer};

    pub fn filter<T>(object: Vec<T>) -> (Vec<T>, Vec<T>)
    {
        let n = object.len()/2;
        let mut Ax: Vec<T> = Vec::with_capacity(n);
        let mut Bx: Vec<T> = Vec::with_capacity(n);
        for (i, val) in object.into_iter().enumerate() {
            match i.is_even() {
                true => Ax.push(val),
                false => Bx.push(val),
            }
        }
        (Ax, Bx)
    }

    pub fn check_pow<T>(n: T) -> bool
    where
        T: Field 
        + BitAnd<T, Output=T>
        + Sub<T, Output=T>,
        <T as BitAnd>::Output: PartialEq<T>
        + Copy,
    {
        n&(n-T::one()) == T::zero()
    }

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

    pub fn de_moivre(element_deg: f64, iter_deg: f64, deg_bound: f64) -> Complex<f64> 
    { 
        let pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;        
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


}


#[cfg(test)]
mod test {
    use crate::poly_rep::{*, fft_subroutines};
    use zksnark::{field::z251::Z251};
    use num::{Num, Complex, Float};

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

    #[test]
    fn bit_reverse() 
    {
        assert_eq!(4, fft_subroutines::bit_reverse(1, 3));
        assert_eq!(5, fft_subroutines::bit_reverse(5, 3));
    }

    #[test]
    fn single_dft_test() {
        let coeffs = vec![ Complex::from(1_f64), Complex::from(2_f64) - Complex::i(), -Complex::i(), -Complex::from(1_f64) + 2_f64 * Complex::i() ];
        let expected_output = Complex::from(-2_f64) - 2.0_f64*Complex::i();
        let dummy_output = Complex::from(-210_f64) - 210.0_f64*Complex::i();

        let fft_vec = fft_subroutines::single_dft(&coeffs, 1_f64, 1_f64).round_to(10_f64);
        assert_eq!(&fft_vec,  &expected_output);
        assert!(&fft_vec !=  &dummy_output);
    }

    #[test]
    fn seq_dft_test() {
        let coeffs = vec![Complex::from(1_f64), Complex::from(2_f64) - Complex::i(), -Complex::i(), -Complex::from(1_f64) + 2_f64 * Complex::i()];
        let expected_output = vec![Complex::from(2_f64), Complex::from(-2_f64) - 2_f64*Complex::i(), -2_f64*Complex::i(), Complex::from(4_f64) + 4_f64*Complex::i()];
        let dummy_output = vec![Complex::from(1_f64), Complex::from(-1_f64) - 2_f64*Complex::i(), -3_f64*Complex::i(), Complex::from(1_f64) + 4_f64*Complex::i()];
        let fft_vec = fft_subroutines::seq_dft(&coeffs).iter().map(|c| c.round_to(10_f64)).collect::<Vec<_>>(); 

        assert_eq!(&fft_vec,  &expected_output);
        assert!(&fft_vec !=  &dummy_output);
    }

    #[test]
    fn check_pow_test() {
        assert!(fft_subroutines::check_pow( Z251::from(8) ));
        assert!(!fft_subroutines::check_pow( Z251::from(5) ));
    }

    #[test]
    fn filter_test() {
        let (Ax, Bx) = fft_subroutines::filter( vec![1, 2, 3, 4] );
        
        assert_eq!(Ax, vec![1, 3]);
        assert_eq!(Bx, vec![2, 4]); 
        assert!(Ax != vec![2, 4]);
        assert!(Bx != vec![1, 3]);
    }
    
    #[test]
    fn pointwise_addition_test() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (0,2), (1,3), (2,18), (3,59) ]);

        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)
    }

    #[test]
    fn pointwise_subtraction_test() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,-5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,-1), (1,-3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (0,0), (1,-3), (2,8), (3,59) ]);
        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)   
    }

    #[test]
    fn pointwise_multiplication_test() {
        let Ax = PointWise::from( vec![ (-3,17), (-2,-3), (-1,1), (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (-3,20), (-2,-3), (-1,2), (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (-3,340), (-2,9), (-1,2), (0,1), (1,0), (2,65), (3,814) ]);
        assert_eq!(Ax.clone() * Bx.clone(), Cx.clone());
        assert!(Ax.clone() * Ax != Cx);
    }
}