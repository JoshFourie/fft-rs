use std::ops::{Add, Mul, Sub};

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
                .map(|(a, b)| {
                    (a.degree, a.y + b.y)
                })
                .collect::<Vec<_>>()
        )
    }
}

#[cfg(test)]
mod test {
    use super::PointWise;
    
    #[test]
    fn pointwise_addition() {
        let Ax = PointWise::from( vec![ (0, 1), (1, 0), (2, 5), (3, 22) ] );
        let Bx = PointWise::from( vec![ (0, 1), (1, 3), (2, 13), (3, 37) ] );
        let Cx = PointWise::from( vec![ (0, 2), (1, 3), (2, 18), (3, 59) ] );
        assert_eq!(true, Ax.clone() + Bx == Cx);
        assert_eq!(false, Ax.clone() + Ax.clone() == Cx)
    }
}