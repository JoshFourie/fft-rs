use crate::two_radix::node::DecimationNode;

#[derive(Debug, PartialEq, Clone)]
pub struct DecimationLeaf
{ 
    pub lhs: Vec<DecimationNode>, 
    pub rhs: Vec<DecimationNode>, 
    pub stage: usize, 
}

impl IntoIterator for DecimationLeaf
{
    type Item = DecimationNode;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter
    {
        self.lhs.into_iter()
            .chain(self.rhs.into_iter())
            .collect::<Vec<_>>()
            .into_iter()
    }
} 

impl DecimationLeaf
{
    pub fn new(lhs: Vec<DecimationNode>, rhs: Vec<DecimationNode>, stage: usize) -> Self { Self{ lhs, rhs, stage } }

    pub fn new_empty() -> Self { Self { lhs: Vec::new(), rhs: Vec::new(), stage: 1 } }

    // takes two leaves and builds them into a single parent leaf.
    // concatenating assigns the stage to the nodes, and indexes them correctly.
    pub fn generate_parent(self, other: Self) -> Self
    {
        assert_eq!(self.stage, other.stage);
        assert_eq!(self.lhs.len(), other.lhs.len());
        assert_eq!(self.rhs.len(), other.rhs.len());
        assert_eq!(self.lhs.len(), self.rhs.len());

        let mut parent = Self::new_empty();
        parent.stage = self.stage+1;
        for (i, (mut lhs, mut rhs)) in self.into_iter()
            .zip( other.into_iter() )
            .enumerate()
        {
            lhs.index=i;
            rhs.index=i;
            lhs.stage+=1;
            rhs.stage+=1;
            lhs.twiddle=false;
            rhs.twiddle=true;
            parent.lhs.push(lhs);
            parent.rhs.push(rhs);
        }
        parent        
    }

    // applies the butterfly input on the decimation nodes.
    pub fn transform(self) -> Self
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

#[cfg(tests)]
mod tests
{
    use super::*;
    use crate::two_radix::node::DecimationNode;
    
    #[test]
    fn test_concatenation()
    {
        let a0 = DecimationNode::new(Complex::from(0.0), 1, 0, false);
        let a4 = DecimationNode::new(Complex::from(4.0), 1, 1, true);
        let a2 = DecimationNode::new(Complex::from(2.0), 1, 0, false);
        let a6 = DecimationNode::new(Complex::from(6.0), 1, 1, true);
        let a1 = DecimationNode::new(Complex::from(1.0), 1, 0, false);
        let a5 = DecimationNode::new(Complex::from(5.0), 1, 1, true);
        let a3 = DecimationNode::new(Complex::from(3.0), 1, 0, false);
        let a7 = DecimationNode::new(Complex::from(7.0), 1, 1, true);
        let concat = DecimationLeaf::generate_parent( 
            DecimationLeaf::new(vec![a0, a4], vec![a2, a6], 1), 
            DecimationLeaf::new(vec![a1, a5], vec![a3, a7], 1) 
        );
        let exp = DecimationLeaf::new(vec![a0, a4, a2, a6], vec![a1, a5, a3, a7], 2);
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
}