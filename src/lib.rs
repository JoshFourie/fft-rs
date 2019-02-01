use std::ops::{Add, Mul, Neg};

// TODO pointwise_with_degree for multiplication


#[derive(Debug, PartialEq, Clone)]
pub struct PointWise<P> { points: Vec<Points<P>> }

#[derive(Debug, PartialEq, Clone)]
pub struct Points<P> { degree: P, y: P }

impl<P> From<Vec<(P, P)>> for PointWise<P> {
    fn from(object: Vec<(P, P)>) -> Self {
        Self {
            points: 
                object.into_iter()
                .map( |(degree,y)| Points::from( (degree, y ) ) )
                .collect::<Vec<_>>()
        }   
    }
}

impl<P> From<(P, P)> for Points<P> { 
    fn from((degree, y): (P, P)) -> Self { Self { degree , y } }
}

impl<P> Add<Self> for PointWise<P> 
where
    P: Add<P, Output=P>,
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
    + Clone,
{
    type Output=Self;
    fn mul(self, rhs: Self) -> Self {
        let n = self
            .clone()
            .points
            .into_iter()
            .zip(rhs.points.clone().into_iter())
            .map( |(a, b)| (a.degree, a.y * b.y) );
        let n_2 = self
            .points
            .into_iter()
            .zip(rhs.points.into_iter())
            .map( |(a, b)| (-a.degree, ( -a.y * -b.y ) ) );
        Self::from( n_2.rev().chain(n).collect::<Vec<_>>() )
    }
}

#[cfg(test)]
mod test {
    use super::PointWise;
    
    #[test]
    fn pointwise_addition() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (0,2), (1,3), (2,18), (3,59) ]);
        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)
    }

    #[test]
    fn pointwise_multiplication() {
        let Ax = PointWise::from( vec![ (0,1), (1,0), (2,5), (3,22) ]);
        let Bx = PointWise::from( vec![ (0,1), (1,3), (2,13), (3,37) ]);
        let Cx = PointWise::from( vec![ (-3,340), (-2,9), (-1,2), (0,1), (1,0), (2,65), (3,814) ]);
        println!("
        {:?}
        
        {:?}", Ax * Bx, Cx);
    }
}