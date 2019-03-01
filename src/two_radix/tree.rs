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
            .enumerate()
            .map(|(i, (lhs,rhs))| 
            {
                let mut leaf=DecimationLeaf::new( Vec::new(), Vec::new(), stage );
                leaf.lhs.push( DecimationNode::new( Complex::from(lhs), stage, i, false ) );
                leaf.rhs.push( DecimationNode::new( Complex::from(rhs), stage, i, true ) );
                leaf
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
                lhs.transform(), rhs.transform()
            );
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
    use super::*;
    use num::Zero;
    use rustfft::FFTplanner;


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
        assert_eq!(tree, output);
    }
}