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

// STROM calculates Stromgren's integral
//   ∫ 0 to x { t^7 exp(2t)/[exp(t)-1]^3 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func STROM(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		FOUR   = 4.0e0
		SEVEN  = 7.0e0
		ONEHUN = 100.0e0
		ONE30  = 130.0e0
		ONE5LN = 0.4055e0
		F15BP4 = 0.38497433455066256959e-1
		PI4B3  = 1.29878788045336582982e2
		VALINF = 196.51956920868988261257e0
	)

	var K1, K2, NTERMS, NUMEXP int
	var EPNGLN, EPSLN, RET, RK, SUMEXP, SUM2, T, X, XHIGH, XK, XK1, XLOW0, XLOW1 float64

	var ASTROM = []float64{
		0.56556120872539155290e0,
		0.4555731969101785525e-1,
		-0.4039535875936869170e-1,
		-0.133390572021486815e-2,
		0.185862506250538030e-2,
		-0.4685555868053659e-4,
		-0.6343475643422949e-4,
		0.572548708143200e-5,
		0.159352812216822e-5,
		-0.28884328431036e-6,
		-0.2446633604801e-7,
		0.1007250382374e-7,
		-0.12482986104e-9,
		-0.26300625283e-9,
		0.2490407578e-10,
		0.485454902e-11,
		-0.105378913e-11,
		-0.3604417e-13,
		0.2992078e-13,
		-0.163971e-14,
		-0.61061e-15,
		0.9335e-16,
		0.709e-17,
		-0.291e-17,
		0.8e-19,
		0.6e-19,
		-0.1e-19}

	X = XVALUE

	//  Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	XK = machine.D1MACH[3]
	T = XK / ONEHUN
	if X <= FOUR {
		NTERMS = 26
		XLOW0 = math.Pow(ONE30*machine.D1MACH[1], ONE/(SEVEN-TWO))
		XLOW1 = TWO * XK
	} else {
		EPSLN = math.Log(machine.D1MACH[4])
		EPNGLN = math.Log(XK)
		XHIGH = SEVEN / XK
	}

	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW0 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, 5) / PI4B3
			} else {
				T = ((X / TWO) - HALF) - HALF
				RET = math.Pow(X, 5) * utils.Cheval(NTERMS, ASTROM, T) * F15BP4
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH {
			SUMEXP = ONE
		} else {
			NUMEXP = int(EPSLN/(ONE5LN-X)) + 1
			if NUMEXP > 1 {
				T = math.Exp(-X)
			} else {
				T = ONE
			}
			RK = ZERO
			for K1 = 1; K1 <= NUMEXP; K1++ {
				RK = RK + ONE
			}
			SUMEXP = ZERO
			for K1 = 1; K1 <= NUMEXP; K1++ {
				SUM2 = ONE
				XK = ONE / (RK * X)
				XK1 = ONE
				for K2 = 1; K2 <= 7; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUM2 = SUM2 * (RK + ONE) / TWO
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = SEVEN*math.Log(X) - X + math.Log(SUMEXP)
		if T < EPNGLN {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)*F15BP4
		}
	}
	return RET
}
