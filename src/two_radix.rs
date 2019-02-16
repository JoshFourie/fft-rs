use num::{ Complex, Num };
use std::ops::{ Div };

/*************************************************** Structures ****************************************/

#[derive(Debug, PartialEq, Clone)] pub struct PointWise<F> { seq: Vec<Complex<F>> }

#[derive(Debug, PartialEq)] pub struct DecimationNode<F> { pair: [Complex<F>; 2] }

#[derive(Debug, PartialEq)] pub struct DecimationTree<F> 
{ 
    root: Vec<F>, 
    nodes: Vec<DecimationNode<F>>, 
}

/*********************************************** Implementations **************************************/

impl<F> From<Vec<F>> for PointWise<F> where F: Clone + Num
{
    fn from(target: Vec<F>) -> Self 
    {
        PointWise{ seq: target.into_iter().map(|num| Complex::from(num)).collect::<Vec<Complex<F>>>() } 
    }
}

impl<F> DecimationNode<F>
{
    pub fn new(lhs: Complex<F>, rhs: Complex<F>) -> Self { Self{ pair: [lhs, rhs] } }

    pub fn concat(&self, rhs: &Self) // -> Self
    {
        
    }
}

impl<F: Copy + Clone> DecimationTree<F>
where
    for <'f> f64: From<&'f F>,
    for <'f> Complex<F>: From<&'f F>
{
    pub fn new_with_root(root: Vec<F>) -> Self { Self { root, nodes: Vec::new() } }

    pub fn init(mut self, bits: usize ) 
    {
        let mut root=self.root;
        danielson_lanczos_pattern(&mut root, bits);
        for (i, lhs) in root.iter().enumerate().step_by(2) 
        {
            let rhs=Complex::from( &root[i+1] );
            let node=DecimationNode::new( Complex::from(lhs), rhs );
            self.nodes.push(node);
        }             
    }
}

/*************************************************** Functions ****************************************/

fn butterfly(x0: f64, x1: f64, n: usize, N: usize) -> (Complex<f64>, Complex<f64>)
{
    let ohm=x1*twiddle(n, N);
    let f0= x0 + ohm;
    let f1= x0 - ohm;
    (f0, f1)
}

fn twiddle(n: usize, N: usize) -> Complex<f64>
{
    let pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;        
    return (
        Complex::from(
        [-2_f64, pi, n as f64, 1_f64.div( N as f64)]
            .iter()
            .product::<f64>()
        ) * Complex::i()
    ).exp()
}

fn danielson_lanczos_pattern<F>(target: &mut Vec<F>, bits: usize)
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

#[cfg(test)]
mod tests 
{
    use super::*;
    
    fn de_moivre(element_deg: f64, iter_deg: f64, deg_bound: f64) -> Complex<f64> 
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

    fn single_dft(coeffs: &Vec<Complex<f64>>, element_deg:f64, max_deg_divisor: f64) -> Complex<f64>
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

    fn seq_dft(coeffs: &Vec<Complex<f64>>) -> Vec<Complex<f64>> 
    {
        coeffs.iter()
            .enumerate()
            .map(|(k, _)| single_dft( coeffs, k as f64, 1_f64 ))
            .collect::<Vec<_>>()
    }

    #[test]
    fn test_dl_pattern_and_bit_rev() 
    {
        let mut target_8_input = vec![0, 1, 2, 3, 4, 5, 6, 7];  
        let expected_8_input = vec![0, 4, 2, 6, 1, 5, 3, 7];
        let dummy_8_input = vec![0, 4, 2, 6, 1, 5, 3, 8];

        // check that positions are swapped.
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