use crate::two_radix::{danielson_lanczos_pattern, node::DecimationNode, leaf::DecimationLeaf};
use num::Complex;
use itertools::Itertools;

#[derive(Debug, PartialEq, Clone)] 
pub struct DecimationTree
{ 
    pub root: Vec<f64>, 
    pub leaves: Vec<DecimationLeaf>, 
}

impl DecimationTree
{
    pub fn new_with_root(root: Vec<f64>) -> Self { Self { root, leaves: Vec::new() } }

    pub fn new_empty() -> Self { Self{ root: Vec::new(), leaves: Vec::new() } }

    pub fn init_base(mut self, bits: usize ) -> Self
    {
        let stage = 1;
        danielson_lanczos_pattern(&mut self.root, bits);
        self.leaves=self.root
            .iter()
            .tuples()
            .map(|(lhs,rhs)| 
            {
                DecimationLeaf::new(
                    vec![DecimationNode::new( Complex::from(lhs), stage, 0, false )],
                    vec![DecimationNode::new( Complex::from(rhs), stage, 0, true )],
                    stage
                ).transform()
            }).collect::<Vec<_>>();
        self
    }

    pub fn grow(self) -> Self
    {
        let mut tree = Self::new_with_root(self.root);
        for (lhs, rhs) in self.leaves
            .into_iter()
            .tuples()
        { 
            let parent=DecimationLeaf::generate_parent(
                lhs, rhs
            ).transform();
            tree.leaves.push(parent);
        }
        tree
    }

    pub fn process(mut self) -> Self
    {
        let discrete_log = {    
            let mut x: u32=0;
            while 2_usize.pow(x) != self.root.len() 
            {
                x += 1;
            }
            x 
        };
        for _ in 0..(discrete_log-1)
        { 
            self=self.grow();
        }
        self
    }

    pub fn extract_seq(self) -> Vec<Complex<f64>>
    {
        assert_eq!(self.leaves.len(), 1);
        let mut seq: Vec<Complex<f64>> = Vec::new();        
        for leaf in self.leaves
            .into_iter() 
        {
            seq.append(
                &mut leaf.lhs
                    .into_iter()
                    .chain(leaf.rhs.into_iter())
                    .map(|node|
                    {
                        node.element
                    }).collect::<Vec<_>>()
            )
        }
        seq
    }
} 

#[cfg(test)]
mod tests
{
    extern crate test;

    use super::*;
    use num::Zero;
    use rustfft::FFTplanner;
    use assert_approx_eq::assert_approx_eq;
    
    #[test]
    fn test_discrete_log() 
    {
        let discrete_log = |target| -> u32
        {   
            let mut x: u32 = 0;
            while 2_usize.pow(x) != target 
            {
                x += 1;
            }
            x 
        };
        assert_eq!(discrete_log(8), 3);
        assert!(discrete_log(8) != 2);
    }

    #[test]
    fn test_process_mapping()
    {
        /* Manual Process */
        let a0 = DecimationNode::new(Complex::from(0.0), 1, 0, false);
        let a4 = DecimationNode::new(Complex::from(4.0), 1, 0, true);
        let a2 = DecimationNode::new(Complex::from(2.0), 1, 0, false);
        let a6 = DecimationNode::new(Complex::from(6.0), 1, 0, true);
        let a1 = DecimationNode::new(Complex::from(1.0), 1, 0, false);
        let a5 = DecimationNode::new(Complex::from(5.0), 1, 0, true);
        let a3 = DecimationNode::new(Complex::from(3.0), 1, 0, false);
        let a7 = DecimationNode::new(Complex::from(7.0), 1, 0, true);

        let leaf_stg1_even_lhs = DecimationLeaf::new(vec![a0], vec![a4], 1).transform();
        let leaf_stg1_even_rhs = DecimationLeaf::new(vec![a2], vec![a6], 1).transform();
        let leaf_stg1_odd_lhs = DecimationLeaf::new(vec![a1], vec![a5], 1).transform();
        let leaf_stg1_odd_rhs = DecimationLeaf::new(vec![a3], vec![a7], 1).transform();  

        let leaf_stg2_lhs = DecimationLeaf::generate_parent(
            leaf_stg1_even_lhs, leaf_stg1_even_rhs
        ).transform();
        let leaf_stg2_rhs = DecimationLeaf::generate_parent(
            leaf_stg1_odd_lhs, leaf_stg1_odd_rhs
        ).transform();        
        let leaf_stg3 = DecimationLeaf::generate_parent(
            leaf_stg2_lhs, leaf_stg2_rhs
        ).transform();

        let root: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let tree = DecimationTree::new_with_root(root).init_base(3).process().leaves[0].clone();
        assert_eq!(leaf_stg3, tree);
    }

    #[test]
    fn test_process_tree()
    {
        let mut input:  Vec<Complex<f64>> = vec![Complex::zero(); 8];
        let mut output: Vec<Complex<f64>> = vec![Complex::zero(); 8];
        input[1] = Complex::new(1.0, 0.0);

        let mut planner = FFTplanner::new(false);
        let fft = planner.plan_fft(8);
        fft.process(&mut input, &mut output);

        let root: Vec<f64> = vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let tree=DecimationTree::new_with_root(root).init_base(3).process().extract_seq();
        for (exp, test) in output.into_iter()
            .zip( tree.into_iter() )
        {
            assert_approx_eq!(exp.re, test.re);
            assert_approx_eq!(exp.im, test.im);
        }

        let mut root = vec![10.0; 4096];
        root[1] = 1.0;
        root[3] = 4.0;

        let mut input:  Vec<Complex<f64>> = vec![Complex::from(10.0); 4096];
        let mut output: Vec<Complex<f64>> = vec![Complex::zero(); 4096];
        input[1] = Complex::new(1.0, 0.0);
        input[3] = Complex::new(4.0, 0.0);

        let mut planner = FFTplanner::new(false);
        let fft = planner.plan_fft(4096);
        fft.process(&mut input, &mut output);

        let tree=DecimationTree::new_with_root(root).init_base(12).process().extract_seq();
        for (exp, test) in output.into_iter()
            .zip( tree.into_iter() )
        {
            assert_approx_eq!(exp.re, test.re);
            assert_approx_eq!(exp.im, test.im);
        }
    }

    #[bench]
    fn stress_bench(b: &mut test::Bencher )
    {
        b.iter(||
        {
            let mut root = vec![10.0; 4096];
            root[1] = 1.0;
            root[3] = 4.0;
            DecimationTree::new_with_root(root).init_base(12).process().extract_seq();
        })
    }

    #[bench]
    fn comparison_bench(b: &mut test::Bencher )
    {
        b.iter(||
        {
            let mut root = vec![Complex::from(10.0); 4096];
            root[1] = Complex::from(1.0);
            root[3] = Complex::from(4.0);
            let mut planner = FFTplanner::new(false);
            let fft = planner.plan_fft(4096);
            fft.process(&mut root.clone(), &mut root);
        })
    }
}