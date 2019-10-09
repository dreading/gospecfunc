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

// ATNINT calculates the value of the inverse-tangent integral defined by
//   ∫ 0 to x ( (arctan t)/t ) dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func ATNINT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		ONEHUN = 100.0e0
		TWOBPI = 0.63661977236758134308e0
	)

	var IND, NTERMS int
	var RET, T, X, XLOW, XUPPER float64
	var ATNINA = []float64{
		1.91040361296235937512e0,
		-0.4176351437656746940e-1,
		0.275392550786367434e-2,
		-0.25051809526248881e-3,
		0.2666981285121171e-4,
		-0.311890514107001e-5,
		0.38833853132249e-6,
		-0.5057274584964e-7,
		0.681225282949e-8,
		-0.94212561654e-9,
		0.13307878816e-9,
		-0.1912678075e-10,
		0.278912620e-11,
		-0.41174820e-12,
		0.6142987e-13,
		-0.924929e-14,
		0.140387e-14,
		-0.21460e-15,
		0.3301e-16,
		-0.511e-17,
		0.79e-18,
		-0.12e-18,
		0.2e-19}

	// Compute the machine-dependent constants.
	T = machine.D1MACH[4] / ONEHUN
	NTERMS = 22
	T = machine.D1MACH[3]
	XLOW = math.Sqrt(T / (ONE + ONE))
	XUPPER = ONE / T

	IND = 1
	X = XVALUE
	if X < ZERO {
		X = -X
		IND = -1
	}
	// Code for X < =  1.0
	if X <= ONE {
		if X < XLOW {
			RET = X
		} else {
			T = X * X
			T = (T - HALF) + (T - HALF)
			RET = X * utils.Cheval(NTERMS, ATNINA, T)
		}
	} else {
		// Code for X > 1.0
		if X > XUPPER {
			RET = math.Log(X) / TWOBPI
		} else {
			T = ONE / (X * X)
			T = (T - HALF) + (T - HALF)
			RET = math.Log(X)/TWOBPI + utils.Cheval(NTERMS, ATNINA, T)/X
		}
	}
	if IND < 0 {
		RET = -RET
	}
	return RET
}
