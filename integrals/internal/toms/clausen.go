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
	"github.com/dreading/gospecfunc/machine"
	"github.com/dreading/gospecfunc/utils"
	"math"
)

// CLAUSN calculates Clausen's integral defined by
//   ∫ 0 to x of (-ln(2*sin(t/2))) dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func CLAUSN(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		ONEHUN = 100.0e0
		TWOPI  = 2 * math.Pi
		PISQ   = 9.8696044010893586188e0
		TWOPIA = 6.28125e0
		TWOPIB = 0.19353071795864769253e-2
	)

	var INDX, NTERMS int
	var RET, T, X, XHIGH, XSMALL float64

	var ACLAUS = []float64{2.14269436376668844709e0,
		0.7233242812212579245e-1,
		0.101642475021151164e-2,
		0.3245250328531645e-4,
		0.133315187571472e-5,
		0.6213240591653e-7,
		0.313004135337e-8,
		0.16635723056e-9,
		0.919659293e-11,
		0.52400462e-12,
		0.3058040e-13,
		0.181969e-14,
		0.11004e-15,
		0.675e-17,
		0.42e-18,
		0.3e-19}

	X = XVALUE

	// Compute the machine-dependent constants.
	T = machine.D1MACH[3]
	XHIGH = ONE / T
	// Error test
	if math.Abs(X) > XHIGH {
		return ZERO
	}
	XSMALL = math.Pi * math.Sqrt(HALF*T)
	T = T / ONEHUN
	NTERMS = 15

	INDX = 1
	if X < ZERO {
		X = -X
		INDX = -1
	}

	//  Argument reduced using simulated extra precision
	if X > TWOPI {
		T = float64(int(X / TWOPI))
		X = (X - T*TWOPIA) - T*TWOPIB
	}
	if X > math.Pi {
		X = (TWOPIA - X) + TWOPIB
		INDX = -INDX
	}
	// Set result to zero if X multiple of PI
	if X == ZERO {
		return ZERO
	}
	// Code for X < XSMALL
	if X < XSMALL {
		RET = X * (ONE - math.Log(X))
	} else {
		// Code for XSMALL < =  X < =  PI
		T = (X*X)/PISQ - HALF
		T = T + T
		if T > ONE {
			T = ONE
		}
		RET = X*utils.Cheval(NTERMS, ACLAUS, T) - X*math.Log(X)
	}
	if INDX < 0 {
		RET = -RET
	}
	return RET
}
