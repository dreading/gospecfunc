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

// GOODST calculates the function defined as
//   ∫ 0 to infinity of ( exp(-u*u)/(u+x) ) du
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func GOODST(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		SIX    = 6.0e0
		ONEHUN = 100.0e0
		GAMBY2 = 0.28860783245076643030e0
		RTPIB2 = 0.88622692545275801365e0
	)

	var NTERM1, NTERM2 int
	var FVAL, RET, T, X, XHIGH, XLOW float64

	var AGOST = []float64{
		0.63106560560398446247e0,
		0.25051737793216708827e0,
		-0.28466205979018940757e0,
		0.8761587523948623552e-1,
		0.682602267221252724e-2,
		-0.1081129544192254677e-1,
		0.169101244117152176e-2,
		0.50272984622615186e-3,
		-0.18576687204100084e-3,
		-0.428703674168474e-5,
		0.1009598903202905e-4,
		-0.86529913517382e-6,
		-0.34983874320734e-6,
		0.6483278683494e-7,
		0.757592498583e-8,
		-0.277935424362e-8,
		-0.4830235135e-10,
		0.8663221283e-10,
		-0.394339687e-11,
		-0.209529625e-11,
		0.21501759e-12,
		0.3959015e-13,
		-0.692279e-14,
		-0.54829e-15,
		0.17108e-15,
		0.376e-17,
		-0.349e-17,
		0.7e-19,
		0.6e-19}

	var AGOSTA = []float64{
		1.81775467984718758767e0,
		-0.9921146570744097467e-1,
		-0.894058645254819243e-2,
		-0.94955331277726785e-3,
		-0.10971379966759665e-3,
		-0.1346694539578590e-4,
		-0.172749274308265e-5,
		-0.22931380199498e-6,
		-0.3127844178918e-7,
		-0.436197973671e-8,
		-0.61958464743e-9,
		-0.8937991276e-10,
		-0.1306511094e-10,
		-0.193166876e-11,
		-0.28844270e-12,
		-0.4344796e-13,
		-0.659518e-14,
		-0.100801e-14,
		-0.15502e-15,
		-0.2397e-16,
		-0.373e-17,
		-0.58e-18,
		-0.9e-19,
		-0.1e-19}

	X = XVALUE

	// Error test

	if X <= ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.

	FVAL = machine.D1MACH[3]
	T = FVAL / ONEHUN
	if X <= TWO {
		NTERM1 = 28
		XLOW = FVAL
	} else {
		NTERM2 = 23
		XHIGH = TWO / FVAL
	}

	// Computation for 0 < x <= 2
	if X <= TWO {
		if X < XLOW {
			RET = -GAMBY2 - math.Log(X)
		} else {
			T = (X - HALF) - HALF
			RET = utils.Cheval(NTERM1, AGOST, T) - math.Exp(-X*X)*math.Log(X)
		}
	} else {
		// Computation for x > 2
		FVAL = RTPIB2 / X
		if X > XHIGH {
			RET = FVAL
		} else {
			T = (SIX - X) / (TWO + X)
			RET = FVAL * utils.Cheval(NTERM2, AGOSTA, T)
		}
	}
	return RET
}
