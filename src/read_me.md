# FFT-RS
											
          DECIMATION IN TIME w/ RECURSION TREES

            a0, a1, a2, a3, a4, a5, a6, a7
                       /      \       
                LHS             RHS
            [a0,a4,a2,a6], [a1,a5,a3,a7]
              /	    	          \
        LHS     RHS                LHS     RHS
      [a0,a4], [a2,a6]           [a1,a5], [a3,a7]



[0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0]

[Complex { re: 0.0, im: 0.0 }, Complex { re: 4.0, im: 0.0 }]

[Complex { re: 2.0, im: 0.0 }, Complex { re: 6.0, im: 0.0 }]

[Complex { re: 1.0, im: 0.0 }, Complex { re: 5.0, im: 0.0 }]

[Complex { re: 3.0, im: 0.0 }, Complex { re: 7.0, im: 0.0 }]
