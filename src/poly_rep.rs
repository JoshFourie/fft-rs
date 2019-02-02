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
    use zksnark::field::Field;
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

    pub fn unity_roots<T>(object: Vec<T>) -> Vec<Complex<f64>>
    {
        let de_moivre = |k, n| -> Complex<f64> {
            Complex::exp( &([
                Complex::from(2.0), 
                Complex::from(std::f64::consts::PI), 
                Complex::i(),
                Complex::from(k),
            ].iter()
                .product::<Complex<f64>>()
                .div( Complex::from(n as f64) ))
            )
            // ENFORCE LEMMA CHECKS
        };
        let n = object.len();
        object.into_iter()
            .enumerate()
            .map( |(k, _)| de_moivre(k as f64, n) )
            .collect::<Vec<Complex<f64>>>()           
    }
}

#[cfg(test)]
mod test {
    use crate::poly_rep::{*, fft_subroutines};
    use zksnark::{field::z251::Z251};
    use num::Complex;

    #[test]
    fn unity_roots_fft () {
        assert_eq!(
            fft_subroutines::unity_roots(vec![0, 1, 2, 3])
                .into_iter()
                .map( |k| (k.re as isize, k.im as isize ))
                .collect::<Vec<(isize, isize)>>(),
            vec![
                Complex::<isize>::from(1), 
                Complex::i(),
                Complex::from(-1), 
                -Complex::i()
            ].into_iter()
                .map(|k| (k.re as isize, k.im as isize ))
                .collect::<Vec<(isize, isize)>>()
        );
        assert_eq!(
            fft_subroutines::unity_roots(vec![0, 1, 2, 3, 4])
                .into_iter()
                .map( |k| (k.re as isize, k.im as isize ))
                .collect::<Vec<(isize, isize)>>(),
            vec![
                Complex::from(1_f32), 
                Complex::new(0.3090_f32, 0.95106_f32),
                Complex::new(-0.8090_f32, 0.5878_f32), 
                Complex::new(-0.8090_f32, -0.5878_f32),
                Complex::new(0.3090_f32, -0.95106_f32)
            ].into_iter()
                .map(|k| (k.re as isize, k.im as isize ))
                .collect::<Vec<(isize, isize)>>()
        );
    }

    #[test]
    fn check_pow_fft() {
        assert!(fft_subroutines::check_pow( Z251::from(8) ));
        assert!(!fft_subroutines::check_pow( Z251::from(5) ));
    }

    #[test]
    fn filter_fft() {
        let (Ax, Bx) = fft_subroutines::filter( vec![1, 2, 3, 4] );
        assert_eq!(Ax, vec![1, 3]);
        assert_eq!(Bx, vec![2, 4]); 
        assert!(Ax != vec![2, 4]);
        assert!(Bx != vec![1, 3]);
    }
    
    #[test]
    fn pointwise_addition() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (0,2), (1,3), (2,18), (3,59) ]);
        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)
    }

    #[test]
    fn pointwise_subtraction() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,-5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,-1), (1,-3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (0,0), (1,-3), (2,8), (3,59) ]);
        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)   
    }

    #[test]
    fn pointwise_multiplication() {
        let Ax = PointWise::from( vec![ (-3,17), (-2,-3), (-1,1), (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (-3,20), (-2,-3), (-1,2), (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (-3,340), (-2,9), (-1,2), (0,1), (1,0), (2,65), (3,814) ]);
        assert_eq!(Ax.clone() * Bx.clone(), Cx.clone());
        assert!(Ax.clone() * Ax != Cx);
    }
}