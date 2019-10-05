// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// The original Fortran code are from
//   ALGORITHM 757, COLLECTED ALGORITHMS FROM ACM.
//   THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
//   VOL. 22, NO. 3, September, 1996, P.  288--301.
//   Allan MacLeod, Dept. of Mathematics and Statistics, University of Paisley

package toms

import (
	"math"
)

// CHEVAL evaluates a Chebyshev series, using the
// Clenshaw method with Reinsch modification, as analysed
// in the paper by Oliver.
//
// INPUT PARAMETERS
//     N - INTEGER - The no. of terms in the sequence
//     A - ARRAY, dimension 0 to N - The coefficients of
//         the Chebyshev series
//     T - The value at which the series is to be
//         evaluated
// References:
//      "An error analysis of the modified Clenshaw method for
//       evaluating Chebyshev and Fourier series" J. Oliver,
//       J.I.M.A., vol. 20, 1977, pp379-391
func CHEVAL(N int, A []float64, T float64) float64 {

	const (
		ZERO = 0.0e0
		HALF = 0.5e0
		TEST = 0.6e0
		TWO  = 2.0e0
	)

	var I int
	var D1, D2, TT, U0, U1, U2, RET float64

	U1 = ZERO

	// If ABS ( T )  < 0.6 use the standard Clenshaw method
	if math.Abs(T) < TEST {
		U0 = ZERO
		TT = T + T
		for I = N; I >= 0; I-- {
			U2 = U1
			U1 = U0
			U0 = TT*U1 + A[I] - U2
		}
		RET = (U0 - U2) / TWO
	} else {
		// If ABS ( T )  > =  0.6 use the Reinsch modification
		D1 = ZERO
		// T > =  0.6 code
		if T > ZERO {
			TT = (T - HALF) - HALF
			TT = TT + TT
			for I = N; I >= 0; I-- {
				D2 = D1
				U2 = U1
				D1 = TT*U2 + A[I] + D2
				U1 = D1 + U2
			}
			RET = (D1 + D2) / TWO
		} else {
			// T < =  -0.6 code
			TT = (T + HALF) + HALF
			TT = TT + TT
			for I = N; I >= 0; I-- {
				D2 = D1
				U2 = U1
				D1 = TT*U2 + A[I] - D2
				U1 = D1 - U2
			}
			RET = (D1 - D2) / TWO
		}
	}
	return RET
}
