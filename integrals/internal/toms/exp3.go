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

// EXP3 calculates the value of ∫ 0 to x (exp(-t*t*t)) dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func EXP3(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		SIXTEN = 16.0e0
		ONEHUN = 100.0e0
		FUNINF = 0.89297951156924921122e0
	)
	var NTERM1, NTERM2 int
	var RET, T, X, XLOW, XUPPER float64

	var AEXP3 = []float64{
		1.26919841422112601434e0,
		-0.24884644638414098226e0,
		0.8052622071723104125e-1,
		-0.2577273325196832934e-1,
		0.759987887307377429e-2,
		-0.203069558194040510e-2,
		0.49083458669932917e-3,
		-0.10768223914202077e-3,
		0.2155172626428984e-4,
		-0.395670513738429e-5,
		0.66992409338956e-6,
		-0.10513218080703e-6,
		0.1536258019825e-7,
		-0.209909603636e-8,
		0.26921095381e-9,
		-0.3251952422e-10,
		0.371148157e-11,
		-0.40136518e-12,
		0.4123346e-13,
		-0.403375e-14,
		0.37658e-15,
		-0.3362e-16,
		0.288e-17,
		-0.24e-18,
		0.2e-19}

	var AEXP3A = []float64{
		1.92704649550682737293e0,
		-0.3492935652048138054e-1,
		0.145033837189830093e-2,
		-0.8925336718327903e-4,
		0.705423921911838e-5,
		-0.66717274547611e-6,
		0.7242675899824e-7,
		-0.878258256056e-8,
		0.116722344278e-8,
		-0.16766312812e-9,
		0.2575501577e-10,
		-0.419578881e-11,
		0.72010412e-12,
		-0.12949055e-12,
		0.2428703e-13,
		-0.473311e-14,
		0.95531e-15,
		-0.19914e-15,
		0.4277e-16,
		-0.944e-17,
		0.214e-17,
		-0.50e-18,
		0.12e-18,
		-0.3e-19,
		0.1e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	T = machine.D1MACH[3]
	XLOW = math.Pow(FOUR*T, ONE/THREE)
	XUPPER = math.Pow(-math.Log(T), ONE/THREE)
	T = T / ONEHUN
	if X <= TWO {
		NTERM1 = 24
	} else {
		NTERM2 = 24
	}
	// Code for XVALUE < =  2
	if X <= TWO {
		if X < XLOW {
			RET = X
		} else {
			T = ((X * X * X / FOUR) - HALF) - HALF
			RET = X * utils.Cheval(NTERM1, AEXP3, T)
		}
	} else {
		// Code for XVALUE > 2
		if X > XUPPER {
			RET = FUNINF
		} else {
			T = ((SIXTEN / (X * X * X)) - HALF) - HALF
			T = utils.Cheval(NTERM2, AEXP3A, T)
			T = T * math.Exp(-X*X*X) / (THREE * X * X)
			RET = FUNINF - T
		}
	}
	return RET
}
