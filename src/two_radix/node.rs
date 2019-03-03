use num::Complex;
use crate::two_radix::butterfly;

#[derive(Debug, PartialEq)] 
pub enum NodeError<'a>
{
    Stage(usize, usize),
    Index(usize, usize),
    Twiddle(&'a mut DecimationNode),
    Len(usize, usize),
}

impl<'a> NodeError<'a>
{
    pub fn flip_node(node: &mut DecimationNode) {
        if node.twiddle==true { node.twiddle=false } 
        else if node.twiddle==false { node.twiddle=true }
    } 
}

#[derive(Debug, PartialEq, Copy, Clone)] 
pub struct DecimationNode
{ 
    pub element: Complex<f64>,
    pub index: usize, 
    pub stage: usize,
    pub twiddle: bool, 
}

impl DecimationNode
{
    pub fn new(element: Complex<f64>, stage: usize, index: usize, twiddle: bool) -> Self 
    {
        Self{ element, index, stage, twiddle} 
    }

    // takes two decimation nodes and returns two that have been correctly twiddled.
    pub fn map_butterflies(mut self, mut other: Self) -> (Self, Self)
    {
        let closure = |mut s: Self, mut o: Self| -> (Self, Self)
        {
            let big_n=(2 as usize).pow(s.stage as u32);
            let small_n=s.index;
            let (lhs, rhs)=butterfly( 
                s.element, 
                o.element, 
                small_n, 
                big_n
            );
            s.element=lhs;
            o.element=rhs;
            (s, o)    
        };
        match Self::butterfly_check(&mut self, &mut other)
        {
            Ok(_) => { },
            Err(err) => match err 
            {
                NodeError::Stage(x, y) => panic!("self.stage: {} != other.stage: {}", x, y),
                NodeError::Index(x, y) => panic!("self.index: {} != other.index: {}", x, y),
                NodeError::Len(x, y) => panic!("self.len: {} != other.len: {}", x, y),
                NodeError::Twiddle(lhs) => {
                    NodeError::flip_node(lhs);
                    match Self::butterfly_check(&mut self, &mut other)
                    {
                        Err(e) => match e
                        {
                            NodeError::Twiddle(rhs) => NodeError::flip_node(rhs),
                            _ => panic!("Unreachable Error: second pass through butterfly")
                        },
                        Ok(_) => { },
                    }
                    return closure(self, other) 
                }
            }
        }
        closure(self, other)    
    }

    pub fn butterfly_check<'a>(&'a mut self, other: &'a mut Self) -> Result<(), NodeError>
    {        
        match self.stage == other.stage
        {
            true => {},
            false => return Err(NodeError::Stage(self.stage, other.stage)),
        }
        match self.index == other.index
        {
            true => {},
            false => return Err(NodeError::Index(self.index, other.index))
        }
        match self.twiddle 
        {
            true => return Err(NodeError::Twiddle(self)),
            false => { },
        }
        match other.twiddle 
        {
            true => { },
            false => return Err(NodeError::Twiddle(self)),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests
{
    use super::*;
    use itertools::Itertools;
    use crate::two_radix::{*, tests::seq_dft, leaf::DecimationLeaf, node::DecimationNode };
    use rustfft::FFTplanner;
    use assert_approx_eq::assert_approx_eq;

    fn assert_approx_vec(x: std::vec::IntoIter<Complex<f64>>, y: std::vec::IntoIter<Complex<f64>>)
    {
        for (exp, test) in x.zip(y)
        {
            assert_approx_eq!(exp.re, test.re);
            assert_approx_eq!(exp.im, test.im);
        }
    }

    #[test]
    fn test_two_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0), 1, 0, false);
        let a1 = DecimationNode::new(Complex::from(1.0), 1, 0, true);
        let leaf = DecimationLeaf::new(vec![a0], vec![a1], 3).transform();
        let test = butterfly(a0.element, a1.element, 0, 8);
        let leaf_lhs = leaf.lhs[0].element;
        let leaf_rhs = leaf.rhs[0].element;
        assert_eq!((leaf_lhs, leaf_rhs), test);

        let (x, y) = test;
        let check = seq_dft(&vec![0.0, 1.0].into_iter().map(|n| Complex::<f64>::from(n)).collect_vec());
        assert_approx_vec(check.into_iter(), vec![x, y].into_iter());
    }

    #[test]
    fn test_four_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0), 1, 0, false);
        let a1 = DecimationNode::new(Complex::from(1.0), 1, 0, true);
        let a2 = DecimationNode::new(Complex::from(2.0), 1, 0, false);
        let a3 = DecimationNode::new(Complex::from(3.0), 1, 0, true); 
        
        /* Stage One */
        let (stage_one_a0, stage_one_a2) = butterfly(a0.element, a2.element, 0, 2);
        let (stage_one_a1, stage_one_a3) = butterfly(a1.element, a3.element, 0, 2);

        /* Stage Two */ 
        let (stage_two_a0, stage_two_a1) = butterfly(stage_one_a0, stage_one_a1, 0, 4);
        let (stage_two_a2, stage_two_a3) = butterfly(stage_one_a2, stage_one_a3, 1, 4);

        let stage_one_temp_leaf_lhs = DecimationLeaf::new(
            vec![a0],
            vec![a2],
            1
        ).transform();
        let stage_one_temp_leaf_rhs = DecimationLeaf::new(
            vec![a1],
            vec![a3],
            1
        ).transform();

        let stage_two_temp_leaf = DecimationLeaf::generate_parent(
            stage_one_temp_leaf_lhs,
            stage_one_temp_leaf_rhs
        ).transform();

        let mut test = stage_two_temp_leaf.lhs.into_iter()
            .chain(stage_two_temp_leaf.rhs.into_iter())
            .map(|node| node.element.round_to(1.0))
            .collect::<Vec<_>>();

        danielson_lanczos_pattern(&mut test, 2);
        let exp = vec![stage_two_a0, stage_two_a1, stage_two_a2, stage_two_a3].into_iter().map(|x|x.round_to(1.0)).collect::<Vec<_>>();
        let mut check = seq_dft( 
            &vec![0.0, 1.0, 2.0, 3.0].into_iter()
                .map(|n| Complex::from(n).round_to(1.0))
                .collect::<Vec<_>>() 
        ).into_iter()
            .map(|x|x.round_to(1.0))
            .collect::<Vec<_>>();
        danielson_lanczos_pattern(&mut check, 2);
        assert_approx_vec(test.clone().into_iter(), exp.into_iter());
        assert_approx_vec(test.into_iter(), check.into_iter());
    }

    #[test]
    fn test_eight_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0), 1, 0, false);
        let a4 = DecimationNode::new(Complex::from(4.0), 1, 0, true);
        let a2 = DecimationNode::new(Complex::from(2.0), 1, 0, false);
        let a6 = DecimationNode::new(Complex::from(6.0), 1, 0, true);
        let a1 = DecimationNode::new(Complex::from(1.0), 1, 0, false);
        let a5 = DecimationNode::new(Complex::from(5.0), 1, 0, true);
        let a3 = DecimationNode::new(Complex::from(3.0), 1, 0, false);
        let a7 = DecimationNode::new(Complex::from(7.0), 1, 0, true);

        /* Stage One */
        let (stg1_a0, stg1_a4)=butterfly(a0.element, a4.element, 0, 2);
        let (stg1_a2, stg1_a6)=butterfly(a2.element, a6.element, 0, 2);
        let (stg1_a1, stg1_a5)=butterfly(a1.element, a5.element, 0, 2);
        let (stg1_a3, stg1_a7)=butterfly(a3.element, a7.element, 0, 2); 

        let leaf_stg1_even_lhs = DecimationLeaf::new(vec![a0], vec![a4], 1).transform();
        let leaf_stg1_even_rhs = DecimationLeaf::new(vec![a2], vec![a6], 1).transform();
        let leaf_stg1_odd_lhs = DecimationLeaf::new(vec![a1], vec![a5], 1).transform();
        let leaf_stg1_odd_rhs = DecimationLeaf::new(vec![a3], vec![a7], 1).transform();  

        // Expected:
        // a0 is equivalent in both & a5 is equivalent  
        assert_eq!(&stg1_a0, &leaf_stg1_even_lhs.lhs[0].element);
        assert_eq!(&stg1_a5, &leaf_stg1_odd_lhs.rhs[0].element);

        // twiddle only the rhs branch.
        assert_eq!(&true, &leaf_stg1_odd_lhs.rhs[0].twiddle);
        assert_eq!(&false, &leaf_stg1_even_lhs.lhs[0].twiddle);
        assert_eq!(&false, &leaf_stg1_odd_lhs.lhs[0].twiddle);
        assert_eq!(&true, &leaf_stg1_even_lhs.rhs[0].twiddle);
        assert!(&leaf_stg1_odd_rhs.rhs[0].twiddle != &leaf_stg1_odd_rhs.lhs[0].twiddle);
        assert!(&leaf_stg1_odd_lhs.rhs[0].twiddle != &leaf_stg1_odd_lhs.lhs[0].twiddle);
        assert!(&leaf_stg1_even_rhs.rhs[0].twiddle != &leaf_stg1_even_rhs.lhs[0].twiddle);
        assert!(&leaf_stg1_even_lhs.rhs[0].twiddle != &leaf_stg1_even_lhs.lhs[0].twiddle);

        // all stages are set to 1 so that we check N=2 when we called transform.
        for (w, (x, (y, z))) in leaf_stg1_even_lhs.clone()
            .into_iter()
            .zip( leaf_stg1_odd_lhs.clone().into_iter() 
            .zip( leaf_stg1_even_rhs.clone()
                .into_iter()
                .zip( leaf_stg1_odd_rhs.clone().into_iter() )))
        {
            assert!( w.stage==1 && x.stage==1 && y.stage==1 && z.stage==1)
        }

        /* Stage Two */ 
        let (stg2_a0, stg2_a2)=butterfly(stg1_a0, stg1_a2, 0, 4);
        let (stg2_a4, stg2_a6)=butterfly(stg1_a4, stg1_a6, 1, 4);
        let (stg2_a1, stg2_a3)=butterfly(stg1_a1, stg1_a3, 0, 4);
        let (stg2_a5, stg2_a7)=butterfly(stg1_a5, stg1_a7, 1, 4);

        let leaf_stg2_lhs = DecimationLeaf::generate_parent(
            leaf_stg1_even_lhs, leaf_stg1_even_rhs
        ).transform();
        let leaf_stg2_rhs = DecimationLeaf::generate_parent(
            leaf_stg1_odd_lhs, leaf_stg1_odd_rhs
        ).transform();

        // Expected:
        // indexes are identical and correct when we zip
        let stg2_test_vec = vec![
            (&stg2_a0, &stg2_a2), (&stg2_a4, &stg2_a6), (&stg2_a1, &stg2_a3), (&stg2_a5, &stg2_a7)
        ];
        for (i, ((lhs, rhs), (lhs_num, rhs_num))) in leaf_stg2_lhs.lhs.iter()
            .zip( leaf_stg2_lhs.rhs.iter() )
            .zip( stg2_test_vec.into_iter() )
            .enumerate()
        {
            assert_eq!(i, lhs.index);
            assert_eq!(lhs.index, rhs.index);
            assert!(i+1 != lhs.index);

            assert_eq!(rhs.stage, 2);
            assert_eq!(lhs.stage, 2);
            
            assert_eq!(&lhs.element, lhs_num);
            assert_eq!(&rhs.element, rhs_num);

            assert_eq!(rhs.twiddle, true);
            assert_eq!(lhs.twiddle, false);
        }

        /* Stage Three */
        let (stg3_a0, stg3_a1)=butterfly(stg2_a0, stg2_a1, 0, 8);
        let (stg3_a4, stg3_a5)=butterfly(stg2_a4, stg2_a5, 1, 8);
        let (stg3_a2, stg3_a3)=butterfly(stg2_a2, stg2_a3, 2, 8);
        let (stg3_a6, stg3_a7)=butterfly(stg2_a6, stg2_a7, 3, 8);

        /* Manual DecimationLeaf */        
        let leaf_stg3 = DecimationLeaf::generate_parent(
            leaf_stg2_lhs, leaf_stg2_rhs
        ).transform();

        // Expected:
        for (i, (lhs, rhs)) in leaf_stg3.lhs.iter()
            .zip(leaf_stg3.rhs.iter())
            .enumerate()
        {            
            assert_eq!(i, lhs.index);
            assert_eq!(lhs.index, rhs.index);
            assert!(i+1 != lhs.index);

            assert_eq!(rhs.stage, 3);
            assert_eq!(lhs.stage, 3);

            assert_eq!(rhs.twiddle, true);
            assert_eq!(lhs.twiddle, false);
        }
        /* Test */
        let test = leaf_stg3.lhs.into_iter()
            .chain(leaf_stg3.rhs.into_iter())
            .map(|node| node.element)
            .collect::<Vec<_>>();

        // We can't assert_eq! as the numbers are trivially different around 10 sig. fig. 
        let mut exp = vec![stg3_a0, stg3_a1, stg3_a2, stg3_a3, stg3_a4, stg3_a5, stg3_a6, stg3_a7];
        danielson_lanczos_pattern(&mut exp, 3);
        
        let mut input: Vec<Complex<f64>> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
            .into_iter()
            .map_into::<Complex<f64>>()
            .collect_vec();
        let mut output: Vec<Complex<f64>> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
            .into_iter()
            .map_into::<Complex<f64>>()
            .collect_vec();
        let mut planner = FFTplanner::new(false);
        let fft = planner.plan_fft(8);
        fft.process(&mut input, &mut output);

        let check = seq_dft( 
            &vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0].into_iter()
                .map(|n| Complex::from(n))
                .collect::<Vec<_>>() 
        );

        assert_approx_vec(output.into_iter(), test.clone().into_iter());
        assert_approx_vec(check.into_iter(), test.clone().into_iter());
        assert_approx_vec(exp.into_iter(), test.into_iter());
    }
}