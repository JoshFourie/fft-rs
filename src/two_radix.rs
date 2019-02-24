use num::{ Complex, Num, ToPrimitive };
use itertools::Itertools;
use std::ops::{ Div };

/*************************************************** Structures ****************************************/

#[derive(Debug, PartialEq, Clone)] pub struct PointWise<F> { seq: Vec<Complex<F>> }

#[derive(Debug, PartialEq, Copy, Clone)] 
pub struct DecimationNode
{ 
    element: Complex<f64>,
    index: usize, 
    stage: usize,
    twiddle: bool, 
}

#[derive(Debug, PartialEq, Clone)]
pub struct DecimationLeaf
{ 
    lhs: Vec<DecimationNode>, 
    rhs: Vec<DecimationNode>, 
    stage: usize 
}

#[derive(Debug, PartialEq, Clone)] 
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
    fn new(lhs: Vec<DecimationNode>, rhs: Vec<DecimationNode>, stage: usize) -> Self { Self{ lhs, rhs, stage } }

    fn new_empty() -> Self { Self { lhs: Vec::new(), rhs: Vec::new(), stage: 1 } }

    // takes two leaves and builds them into a single parent leaf.
    // concatenating assigns the stage to the nodes, and indexes them correctly.
    fn generate_parent(self, other: Self) -> Self
    {
        let mut parent = Self::new_empty();
        parent.stage = self.stage+1;
        for (i, (mut lhs, mut rhs)) in self.lhs
            .into_iter()
            .chain( self.rhs.into_iter() )
            .zip( other.lhs.into_iter().chain( other.rhs.into_iter() ) )
            .enumerate()
        {
            lhs.index = i;
            rhs.index = i;
            lhs.stage = parent.stage;
            rhs.stage = parent.stage;
            lhs.twiddle = false;
            rhs.twiddle = true;
            parent.lhs.push(lhs);
            parent.rhs.push(rhs);
        };
        parent        
    }

    // applies the butterfly input on the decimation nodes.
    fn transform(self) -> Self
    {
        let mut leaf = Self::new_empty();
        for (lhs, rhs) in self.lhs
            .into_iter()
            .zip( self.rhs.into_iter() )
        {
            let (lhs_node, rhs_node)=lhs.map_butterflies(rhs);
            leaf.lhs.push(lhs_node);
            leaf.rhs.push(rhs_node);
        }
        leaf
    }
}

impl DecimationNode
{
    fn new(element: Complex<f64>) -> Self { Self{ element, index: 0, stage: 1, twiddle: false } }

    // takes two decimation nodes and returns two that have been correctly twiddled.
    fn map_butterflies(mut self, mut other: Self) -> (Self, Self)
    {
        assert!(self.stage == other.stage);
        assert!(!self.twiddle, other.twiddle);

        let big_n=(2 as usize).pow(self.stage as u32);
        let small_n=self.index;
        let (lhs, rhs)=butterfly( 
            self.element, 
            other.element, 
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
                let mut leaf=DecimationLeaf::new( Vec::new(), Vec::new(), stage );
                leaf.lhs.push( DecimationNode::new( Complex::from(lhs) ) );
                leaf.rhs.push( DecimationNode::new( Complex::from(rhs) ) );
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
            let parent=lhs.transform()
                .generate_parent( rhs.transform() );
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
        assert!(self.leaves.len() == 1);
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

/*************************************************** Functions ****************************************/

fn butterfly(lhs: Complex<f64>, rhs: Complex<f64>, small_n: usize, big_n: usize) -> (Complex<f64>, Complex<f64>)
{
    let ohm=rhs*twiddle(small_n, big_n);
    (lhs + ohm, lhs - ohm)
}

fn twiddle(small_n: usize, big_n: usize) -> Complex<f64>
{
    let pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;        
    return (
        Complex::from(
        [-2.0, pi, small_n as f64, 1.0.div( big_n as f64)]
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
    use num::Float;
    use std::ops::Mul;

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
    fn test_concatenation()
    {
        let a0 = DecimationNode::new(Complex::from(0.0));
        let a4 = DecimationNode::new(Complex::from(4.0));
        let a2 = DecimationNode::new(Complex::from(2.0));
        let a6 = DecimationNode::new(Complex::from(6.0));
        let a1 = DecimationNode::new(Complex::from(1.0));
        let a5 = DecimationNode::new(Complex::from(5.0));
        let a3 = DecimationNode::new(Complex::from(3.0));
        let a7 = DecimationNode::new(Complex::from(7.0));
        let concat = DecimationLeaf::generate_parent( 
            DecimationLeaf::new(vec![a0, a4], vec![a2, a6], 2), 
            DecimationLeaf::new(vec![a1, a5], vec![a3, a7], 2) 
        );
        let exp = DecimationLeaf::new(vec![a0, a4, a2, a6], vec![a1, a5, a3, a7], 3);
        let conc_lhs = concat.lhs
            .into_iter()
            .map(|node| node.element)
            .collect::<Vec<_>>();
        let conc_rhs = concat.rhs
            .into_iter()
            .map(|node| node.element)
            .collect::<Vec<_>>();
        let exp_lhs = exp.lhs
            .into_iter()
            .map(|node| node.element)
            .collect::<Vec<_>>();
        let exp_rhs = exp.rhs
            .into_iter()
            .map(|node| node.element)
            .collect::<Vec<_>>();
        assert_eq!(exp_lhs, conc_lhs);
        assert_eq!(exp_rhs, conc_rhs);
    }

    #[test]
    fn test_two_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0));
        let a1 = DecimationNode::new(Complex::from(1.0));
        let leaf = DecimationLeaf::new(vec![a0], vec![a1], 1).transform();
        let exp = butterfly(a0.element, a1.element, 0, 2);
        let leaf_lhs = leaf.lhs[0].element;
        let leaf_rhs = leaf.rhs[0].element;
        assert_eq!((leaf_lhs, leaf_rhs), exp);
    }

    #[test]
    fn test_four_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0));
        let a1 = DecimationNode::new(Complex::from(1.0));
        let a2 = DecimationNode::new(Complex::from(2.0));
        let a3 = DecimationNode::new(Complex::from(3.0)); 
        
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
            .map(|node| node.element.round_to(10.0))
            .collect::<Vec<_>>();

        danielson_lanczos_pattern(&mut test, 2);
        let exp = vec![stage_two_a0, stage_two_a1, stage_two_a2, stage_two_a3].into_iter().map(|x|x.round_to(10.0)).collect::<Vec<_>>();
        let mut check = seq_dft( 
            &vec![0.0, 1.0, 2.0, 3.0].into_iter()
                .map(|n| Complex::from(n))
                .collect::<Vec<_>>() 
        ).into_iter()
            .map(|x|x.round_to(10.0))
            .collect::<Vec<_>>();
        danielson_lanczos_pattern(&mut check, 2);
        assert_eq!(test, exp);
        assert_eq!(exp, check);
    }

    #[test]
    fn test_eight_input_butterfly_map()
    {
        let a0 = DecimationNode::new(Complex::from(0.0));
        let a4 = DecimationNode::new(Complex::from(4.0));
        let a2 = DecimationNode::new(Complex::from(2.0));
        let a6 = DecimationNode::new(Complex::from(6.0));
        let a1 = DecimationNode::new(Complex::from(1.0));
        let a5 = DecimationNode::new(Complex::from(5.0));
        let a3 = DecimationNode::new(Complex::from(3.0));
        let a7 = DecimationNode::new(Complex::from(7.0));

        /* Stage One */
        let (stg1_a0, stg1_a4)=butterfly(a0.element, a4.element, 0, 2);
        let (stg1_a2, stg1_a6)=butterfly(a2.element, a6.element, 0, 2);
        let (stg1_a1, stg1_a5)=butterfly(a1.element, a5.element, 0, 2);
        let (stg1_a3, stg1_a7)=butterfly(a3.element, a7.element, 0, 2); 

        /* Stage Two */ 
        let (stg2_a0, stg2_a2)=butterfly(stg1_a0, stg1_a2, 0, 4);
        let (stg2_a4, stg2_a6)=butterfly(stg1_a4, stg1_a6, 1, 4);
        let (stg2_a1, stg2_a3)=butterfly(stg1_a1, stg1_a3, 0, 4);
        let (stg2_a5, stg2_a7)=butterfly(stg1_a5, stg1_a7, 1, 4);

        /* Stage Three */
        let (stg3_a0, stg3_a1)=butterfly(stg2_a0, stg2_a1, 0, 8);
        let (stg3_a4, stg3_a5)=butterfly(stg2_a4, stg2_a5, 1, 8);
        let (stg3_a2, stg3_a3)=butterfly(stg2_a2, stg2_a3, 2, 8);
        let (stg3_a6, stg3_a7)=butterfly(stg2_a6, stg2_a7, 3, 8);

        /* Manual DecimationLeaf */
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

        /* Test */
        let exp = vec![stg3_a0, stg3_a1, stg3_a2, stg3_a3, stg3_a4, stg3_a5, stg3_a6, stg3_a7].into_iter()
            .map(|x|x.round_to(10.0))
            .collect::<Vec<_>>();
        let mut check = seq_dft( 
            &vec![0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0]
                .into_iter()
                .map(|n| Complex::from(n))
                .collect::<Vec<_>>() 
            ).into_iter()
            .map(|x|x.round_to(10.0))
            .collect::<Vec<_>>();
        let mut test = leaf_stg3.lhs.into_iter()
            .chain(leaf_stg3.rhs.into_iter())
            .map(|node| node.element)
            .collect::<Vec<_>>();   
        danielson_lanczos_pattern(&mut test, 3);
        danielson_lanczos_pattern(&mut check, 3);
        assert_eq!(check, exp);
    }

    #[test]
    fn test_process_tree()
    {
        let root: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let mut expected=seq_dft(&root.clone().into_iter().map_into::<Complex<f64>>().collect::<Vec<_>>());
        let tree=DecimationTree::new_with_root(root).init_base(3);
        danielson_lanczos_pattern(&mut expected, 3);
        assert_eq!(expected, tree.process().extract_seq());
    }
}