use num::{ Complex, Num, ToPrimitive };
use std::ops::{ Div };

/*************************************************** Structures ****************************************/

#[derive(Debug, PartialEq, Clone)] pub struct PointWise<F> { seq: Vec<Complex<F>> }

#[derive(Debug, PartialEq)] 
pub struct DecimationNode
{ 
    element: Complex<f64>,
    index: usize, 
    stage: usize,
    twiddle: bool, 
}

#[derive(Debug, PartialEq)]
pub struct DecimationLeaf
{ 
    lhs: Vec<DecimationNode>, 
    rhs: Vec<DecimationNode>, 
    stage: usize 
}

#[derive(Debug, PartialEq)] 
pub struct DecimationTree
{ 
    root: Vec<f64>, 
    leaves: Vec<DecimationLeaf>, 
}

/*********************************************** Implementations **************************************/

impl<F> From<Vec<F>> for PointWise<F> where F: Clone + Num
{
    fn from(target: Vec<F>) -> Self 
    {
        PointWise{ seq: target.into_iter().map(|num| Complex::from(num)).collect::<Vec<Complex<F>>>() } 
    }
}

impl DecimationLeaf
{
    pub fn new(lhs: Vec<DecimationNode>, rhs: Vec<DecimationNode>) -> Self { Self{ lhs, rhs, stage: 1 } }

    pub fn new_empty() -> Self { Self { lhs: Vec::new(), rhs: Vec::new(), stage: 1 } }

    // takes two leaves and builds them into a single parent leaf.
    // concatenating assigns the stage to the nodes, and indexes them correctly.
    pub fn generate_parent(mut self, mut other: Self) -> Self
    {
        assert!(self.stage == other.stage);
        assert!(self.lhs.len() == self.rhs.len());
        assert!(self.lhs.len() == other.lhs.len());
        assert!(self.rhs.len() == other.rhs.len());

        self.lhs.append(&mut self.rhs);
        other.lhs.append(&mut self.rhs);

        let mut parent=Self::new_empty();
        parent.stage=self.stage+1;
        parent.lhs=self.lhs.into_iter()
            .enumerate()
            .map(|(i, mut node)|
            {
                node.stage=parent.stage;
                node.twiddle=false;
                node.index=i;
                node
            }).collect::<Vec<_>>();
        parent.rhs=other.lhs.into_iter()
            .enumerate()
            .map(|(i, mut node)| 
            {
                node.stage=parent.stage;
                node.twiddle=false;
                node.index=i;
                node
            }).collect::<Vec<_>>();        
        parent
    }
}

impl DecimationNode
{
    pub fn new(element: Complex<f64>) -> Self { Self{ element, index: 0, stage: 1, twiddle: false } }

    // takes two decimation nodes and returns two that have been correctly twiddled.
    pub fn map_butterflies(mut self, mut other: Self) -> (Self, Self) 
    {
        assert!(self.stage == other.stage);
        assert!(!self.twiddle, other.twiddle);

        let big_n=(2 as usize).pow(self.stage as u32);
        let small_n=self.index;
        let (lhs, rhs)=butterfly( 
            self.element.to_f64()
                .unwrap(), 
            other.element.to_f64()
                .unwrap(), 
            small_n, 
            big_n
        );
        self.element=lhs;
        other.element=rhs;
        (self, other)
    }
}

impl DecimationTree
{
    pub fn new_with_root(root: Vec<f64>) -> Self { Self { root, leaves: Vec::new() } }

    pub fn init(mut self, bits: usize ) -> Self
    {
        danielson_lanczos_pattern(&mut self.root, bits);
        self.leaves=self.root
            .iter()
            .enumerate()
            .step_by(2)
            .map(|(i, lhs)| 
            {
                let mut leaf=DecimationLeaf::new( Vec::new(), Vec::new() );
                let rhs=self.root[i+1];
                leaf.lhs.push( DecimationNode::new( Complex::from(lhs) ) );
                leaf.rhs.push( DecimationNode::new( Complex::from(rhs) ) );
                leaf
            }).collect::<Vec<_>>();
        self
    }

    pub fn process_tree(mut self) -> Vec<Complex<f64>>
    {
        
    }
} 

/*************************************************** Functions ****************************************/

fn butterfly(lhs: f64, rhs: f64, n: usize, N: usize) -> (Complex<f64>, Complex<f64>)
{
    let ohm=rhs*twiddle(n, N);
    (lhs + ohm, lhs - ohm)
}

fn twiddle(n: usize, N: usize) -> Complex<f64>
{
    let pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;        
    return (
        Complex::from(
        [-2.0, pi, n as f64, 1.0.div( N as f64)]
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

    #[test]
    fn test_decimation_tree()
    {
        let root: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let tree=DecimationTree::new_with_root(root).init(3);
       
        // todo: change to assert_eq!
        println!("{:?}", tree);
    }

    #[test]
    fn test_decimation_node_concatenation()
    {
        let mut root: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        danielson_lanczos_pattern(&mut root, 3);
        let tree=DecimationTree::new_with_root(root).init(3);
        let concat=&tree.nodes[2].concatenate(&tree.nodes[3]);
        println!( "{:?}",  concat);
        println!( "{:?}", &tree.nodes);   
    }
}