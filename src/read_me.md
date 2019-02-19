# FFT-RS: DECIMATION IN TIME w/ RECURSION TREES
               
34.5

## Preliminaries

A recursion tree bit-reverses and decomposes a vector input. Note that this isn't a very useful visualisation for the actual butterfly proccess, but is useful for the particular straegy:

            a0, a1, a2, a3, a4, a5, a6, a7
                       /      \       
                LHS             RHS
            [a0,a4,a2,a6], [a1,a5,a3,a7]
              /	    	          \
        LHS     RHS                LHS     RHS
      [a0,a4], [a2,a6]           [a1,a5], [a3,a7]
       /   \    /   \             /   \    /   \
     a0    a4  a2    a6          a1   a5  a3    a7


    
The two-input butterfly equation ('Butterfly-DFT') produces two outputs, where _Twiddle_ represents the exp( 2 * pi * i * n *(1/N) ) : 

> **F(A) = xL + (xR * Twiddle)** 
> && 
>**F(B) = xL - (xR * Twiddle)**

These can be combined to provide a programmatic butterfly solution. 
* * *

## Strategy

The goal is to butterfly recursively such that if the recursive tree decomposes into N stages, then a given input is also 'twiddled' N times.

### Extrapolating from the 8-Input Instance

We can examine the process for a few factors run through an 8-input butterfly for information.

**Observation One: Mapping the variables.**

_Stage One._ 
The first application of the Butterfly-DFT maps a2 -> a6, which is equivalent to taking the vectors at their index (input[ i ], input[ i+1 ]) and running them through the butterfly DFT.

_Stage Two._ 
The next application maps a2-> a0, examining the recursion tree, we can see that zipping the two inputs lets us map a0 -> a2 with relative ease, and also lets us map a4-a6 correctly.

_Stage Three._ 
Conveniently, the same process from stage two lets us map a0 -> a1, a4 -> a5 and so on.

**Observation Two: Applying the Twiddle Factor.**

The Twiddle factor needs to adjust at each stage for the *n* and *N* variables in our equation. The *n* variable iterates upwards according to the index of the element and we can derive *N* from a simple 2^(stage) stored within the element. 

Consider the movement of a3 through the 8-input butterfly: ( stg = 1, i=0, none ) ->  ( stg=2, i=0, W(n=0, N=4) )-> ( stg=3, i=1, W(n=1, N=8) )

**Strategic Rules.** 
We can determine a set of rules from these two observations.
> 1. The variables can be appropriately mapped according to their index at a given stage within the vector.
> 2. The Twiddle factor can be determined by the index of an element within the vector and the relevant stage.

Note that both rules are dependent on the relative position of the element in the vector of our recursion tree, whether as a part of the RHS input, or the LHS input, and a persistent variable recording the 'stage' of the process. 

**Indexing.**
We are modifying each element consistently, so the recursion tree can burn the leaf nodes in a concatenation function for memory efficiency. Note that this requires a solution to the Discrete Log Problem.