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

// SYNCH1 calculates the synchrotron radiation function
//   x ∫ 0 to infinity { K(5/3)(t) } dt
// where K(5/3) is a modified Bessel function of order 5/3.
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func SYNCH1(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		TWELVE = 12.0e0
		ONEHUN = 100.0e0
		CONLOW = 2.14952824153447863671e0
		PIBRT3 = 1.81379936423421785059e0
		LNRTP2 = 0.22579135264472743236e0
	)
	var NTERM1, NTERM2, NTERM3 int
	var CHEB1, CHEB2, RET, T, X, XHIGH1, XHIGH2, XLOW, XPOWTH float64

	var ASYNC1 = []float64{
		30.36468298250107627340e0,
		17.07939527740839457449e0,
		4.56013213354507288887e0,
		0.54928124673041997963e0,
		0.3729760750693011724e-1,
		0.161362430201041242e-2,
		0.4819167721203707e-4,
		0.105124252889384e-5,
		0.1746385046697e-7,
		0.22815486544e-9,
		0.240443082e-11,
		0.2086588e-13,
		0.15167e-15,
		0.94e-18}

	var ASYNC2 = []float64{
		0.44907216235326608443e0,
		0.8983536779941872179e-1,
		0.810445737721512894e-2,
		0.42617169910891619e-3,
		0.1476096312707460e-4,
		0.36286336153998e-6,
		0.666348074984e-8,
		0.9490771655e-10,
		0.107912491e-11,
		0.1002201e-13,
		0.7745e-16,
		0.51e-18}

	var ASYNCA = []float64{
		2.13293051613550009848e0,
		0.7413528649542002401e-1,
		0.869680999099641978e-2,
		0.117038262487756921e-2,
		0.16451057986191915e-3,
		0.2402010214206403e-4,
		0.358277563893885e-5,
		0.54477476269837e-6,
		0.8388028561957e-7,
		0.1306988268416e-7,
		0.205309907144e-8,
		0.32518753688e-9,
		0.5179140412e-10,
		0.830029881e-11,
		0.133527277e-11,
		0.21591498e-12,
		0.3499673e-13,
		0.569942e-14,
		0.92906e-15,
		0.15222e-15,
		0.2491e-16,
		0.411e-17,
		0.67e-18,
		0.11e-18,
		0.2e-19}

	// Start calculation
	X = XVALUE
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	CHEB1 = machine.D1MACH[3]
	T = CHEB1 / ONEHUN
	if X <= FOUR {
		NTERM1 = 13
		NTERM2 = 11
		XLOW = math.Sqrt(EIGHT * CHEB1)
	} else {
		NTERM3 = 24
		XHIGH2 = math.Log(machine.D1MACH[1])
		XHIGH1 = -EIGHT * XHIGH2 / (EIGHT - ONE)
	}

	// Code for 0 <= x <= 4
	if X <= FOUR {
		XPOWTH = math.Pow(X, ONE/THREE)
		if X < XLOW {
			RET = CONLOW * XPOWTH
		} else {
			T = (X*X/EIGHT - HALF) - HALF
			CHEB1 = utils.Cheval(NTERM1, ASYNC1, T)
			CHEB2 = utils.Cheval(NTERM2, ASYNC2, T)
			T = XPOWTH*CHEB1 - math.Pow(XPOWTH, 11)*CHEB2
			RET = T - PIBRT3*X
		}
	} else {
		if X > XHIGH1 {
			RET = ZERO
		} else {
			T = (TWELVE - X) / (X + FOUR)
			CHEB1 = utils.Cheval(NTERM3, ASYNCA, T)
			T = LNRTP2 - X + math.Log(math.Sqrt(X)*CHEB1)
			if T < XHIGH2 {
				RET = ZERO
			} else {
				RET = math.Exp(T)
			}
		}
	}
	return RET
}

// SYNCH2 calculates the synchrotron radiation function
//   x * K(2/3)(x)
// where K(2/3) is a modified Bessel function of order 2/3.
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func SYNCH2(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		TEN    = 10.0e0
		ONEHUN = 100.0e0
		CONLOW = 1.07476412076723931836e0
		LNRTP2 = 0.22579135264472743236e0
	)
	var NTERM1, NTERM2, NTERM3 int
	var CHEB1, CHEB2, RET, T, X, XHIGH1, XHIGH2, XLOW, XPOWTH float64

	var ASYN21 = []float64{
		38.61783992384308548014e0,
		23.03771559496373459697e0,
		5.38024998683357059676e0,
		0.61567938069957107760e0,
		0.4066880046688955843e-1,
		0.172962745526484141e-2,
		0.5106125883657699e-4,
		0.110459595022012e-5,
		0.1823553020649e-7,
		0.23707698034e-9,
		0.248872963e-11,
		0.2152868e-13,
		0.15607e-15,
		0.96e-18,
		0.1e-19}

	var ASYN22 = []float64{
		7.90631482706608042875e0,
		3.13534636128534256841e0,
		0.48548794774537145380e0,
		0.3948166758272372337e-1,
		0.196616223348088022e-2,
		0.6590789322930420e-4,
		0.158575613498559e-5,
		0.2868653011233e-7,
		0.40412023595e-9,
		0.455684443e-11,
		0.4204590e-13,
		0.32326e-15,
		0.210e-17,
		0.1e-19}

	var ASYN2A = []float64{
		2.02033709417071360032e0,
		0.1095623712180740443e-1,
		0.85423847301146755e-3,
		0.7234302421328222e-4,
		0.631244279626992e-5,
		0.56481931411744e-6,
		0.5128324801375e-7,
		0.471965329145e-8,
		0.43807442143e-9,
		0.4102681493e-10,
		0.386230721e-11,
		0.36613228e-12,
		0.3480232e-13,
		0.333010e-14,
		0.31856e-15,
		0.3074e-16,
		0.295e-17,
		0.29e-18,
		0.3e-19}

	X = XVALUE
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	CHEB1 = machine.D1MACH[3]
	T = CHEB1 / ONEHUN
	if X <= FOUR {
		NTERM1 = 14
		NTERM2 = 13
		XLOW = math.Sqrt(EIGHT * CHEB1)
	} else {
		NTERM3 = 18
		XHIGH2 = math.Log(machine.D1MACH[1])
		XHIGH1 = -EIGHT * XHIGH2 / (EIGHT - ONE)
	}

	// Code for 0 <= x <= 4

	if X <= FOUR {
		XPOWTH = math.Pow(X, ONE/THREE)
		if X < XLOW {
			RET = CONLOW * XPOWTH
		} else {
			T = (X*X/EIGHT - HALF) - HALF
			CHEB1 = utils.Cheval(NTERM1, ASYN21, T)
			CHEB2 = utils.Cheval(NTERM2, ASYN22, T)
			RET = XPOWTH*CHEB1 - math.Pow(XPOWTH, 5)*CHEB2
		}
	} else {
		if X > XHIGH1 {
			RET = ZERO
		} else {
			T = (TEN - X) / (X + TWO)
			CHEB1 = utils.Cheval(NTERM3, ASYN2A, T)
			T = LNRTP2 - X + math.Log(math.Sqrt(X)*CHEB1)
			if T < XHIGH2 {
				RET = ZERO
			} else {
				RET = math.Exp(T)
			}
		}
	}
	return RET
}
