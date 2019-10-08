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

// LOBACH calculates the the Lobachewsky function L(x)
//   ∫ 0 to x of -ln | cos t | dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func LOBACH(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		SIX    = 6.0e0
		TEN    = 10.0e0
		ONEHUN = 100.0e0
		LOBPIA = 1115.0e0
		LOBPIB = 512.0e0
		LOBPI2 = -1.48284696397869499311e-4
		LBPB22 = -7.41423481989347496556e-5
		PI11   = 201.0e0
		PI12   = 64.0e0
		PI2    = 9.67653589793238462643e-4
		PIBY22 = 4.83826794896619231322e-4
		TCON   = 3.24227787655480868620e0
	)
	var INDPI2, INDSGN, NTERM1, NTERM2 int
	var FVAL, FVAL1, LBPB21, NPI, PI, PIBY2, LOBPI1, PIBY21, PIBY4, PI1,
		RET, T, X, XCUB, XHIGH, XLOW1, XLOW2, XLOW3, XR float64

	var ARLOB1 = []float64{
		0.34464884953481300507e0,
		0.584198357190277669e-2,
		0.19175029694600330e-3,
		0.787251606456769e-5,
		0.36507477415804e-6,
		0.1830287272680e-7,
		0.96890333005e-9,
		0.5339055444e-10,
		0.303408025e-11,
		0.17667875e-12,
		0.1049393e-13,
		0.63359e-15,
		0.3878e-16,
		0.240e-17,
		0.15e-18,
		0.1e-19}

	var ARLOB2 = []float64{
		2.03459418036132851087e0,
		0.1735185882027407681e-1,
		0.5516280426090521e-4,
		0.39781646276598e-6,
		0.369018028918e-8,
		0.3880409214e-10,
		0.44069698e-12,
		0.527674e-14,
		0.6568e-16,
		0.84e-18,
		0.1e-19}

	X = math.Abs(XVALUE)
	INDSGN = 1
	if XVALUE < ZERO {
		INDSGN = -1
	}

	// Compute the machine-dependent constants.

	XR = machine.D1MACH[3]
	XHIGH = ONE / XR

	// Error test
	if X > XHIGH {
		return ZERO
	}

	T = XR / ONEHUN
	NTERM1 = 15
	NTERM2 = 10
	XLOW1 = math.Pow(SIX*machine.D1MACH[1], TWO/SIX)
	XLOW2 = math.Sqrt(TEN * XR)
	T = TWO*TEN - TWO
	XLOW3 = math.Sqrt(T * XR)

	// Reduce argument to [0,pi]

	PI1 = PI11 / PI12
	PI = PI1 + PI2
	PIBY2 = PI / TWO
	PIBY21 = PI1 / TWO
	PIBY4 = PIBY2 / TWO
	NPI = float64(int(X / PI))
	XR = (X - NPI*PI1) - NPI*PI2

	// Reduce argument to [0,pi/2]
	INDPI2 = 0
	if XR > PIBY2 {
		INDPI2 = 1
		XR = (PI1 - XR) + PI2
	}

	// Code for argument in [0,pi/4]

	if XR <= PIBY4 {
		if XR < XLOW1 {
			FVAL = ZERO
		} else {
			XCUB = XR * XR * XR
			if XR < XLOW2 {
				FVAL = XCUB / SIX
			} else {
				T = (TCON*XR*XR - HALF) - HALF
				FVAL = XCUB * utils.Cheval(NTERM1, ARLOB1, T)
			}
		}
	} else {
		// Code for argument in [pi/4,pi/2]
		XR = (PIBY21 - XR) + PIBY22
		if XR == ZERO {
			FVAL1 = ZERO
		} else {
			if XR < XLOW3 {
				FVAL1 = XR * (ONE - math.Log(XR))
			} else {
				T = (TCON*XR*XR - HALF) - HALF
				FVAL1 = XR * (utils.Cheval(NTERM2, ARLOB2, T) - math.Log(XR))
			}
		}
		LBPB21 = LOBPIA / (LOBPIB + LOBPIB)
		FVAL = (LBPB21 - FVAL1) + LBPB22
	}
	LOBPI1 = LOBPIA / LOBPIB
	// Compute value for argument in [pi/2,pi]
	if INDPI2 == 1 {
		FVAL = (LOBPI1 - FVAL) + LOBPI2
	}
	RET = FVAL

	// Scale up for arguments > pi
	if NPI > 0 {
		RET = (FVAL + NPI*LOBPI2) + NPI*LOBPI1
	}
	if INDSGN == -1 {
		RET = -RET
	}
	return RET
}
