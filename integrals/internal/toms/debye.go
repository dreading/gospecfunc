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

// DEBYE1 calculates the Debye function of order 1, defined as
//   ∫ 0 to x of t/(exp(t)-1) dt] / x
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func DEBYE1(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		QUART  = 0.25e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		NINE   = 9.0e0
		THIRT6 = 36.0e0
		ONEHUN = 100.0e0
		DEBINF = 0.60792710185402662866e0
	)

	var I, NEXP, NTERMS int
	var EXPMX, RET, RK, SUM, T, X, XK, XLIM, XLOW, XUPPER float64

	var ADEB1 = []float64{
		2.40065971903814101941e0,
		0.19372130421893600885e0,
		-0.623291245548957703e-2,
		0.35111747702064800e-3,
		-0.2282224667012310e-4,
		0.158054678750300e-5,
		-0.11353781970719e-6,
		0.835833611875e-8,
		-0.62644247872e-9,
		0.4760334890e-10,
		-0.365741540e-11,
		0.28354310e-12,
		-0.2214729e-13,
		0.174092e-14,
		-0.13759e-15,
		0.1093e-16,
		-0.87e-18,
		0.7e-19,
		-0.1e-19}

	X = XVALUE

	// Check XVALUE >= 0.0
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	T = machine.D1MACH[3]
	XLOW = math.Sqrt(T * EIGHT)
	XUPPER = -math.Log(T + T)
	XLIM = -math.Log(machine.D1MACH[1])
	T = T / ONEHUN
	NTERMS = 18
	// Code for x <= 4.0
	if X <= FOUR {
		if X < XLOW {
			RET = ((X-NINE)*X + THIRT6) / THIRT6
		} else {
			T = ((X * X / EIGHT) - HALF) - HALF
			RET = utils.Cheval(NTERMS, ADEB1, T) - QUART*X
		}
	} else {
		// Code for x > 4.0
		RET = ONE / (X * DEBINF)
		if X < XLIM {
			EXPMX = math.Exp(-X)
			if X > XUPPER {
				RET = RET - EXPMX*(ONE+ONE/X)
			} else {
				SUM = ZERO
				RK = float64(int(XLIM / X))
				NEXP = int(RK)
				XK = RK * X
				for I = NEXP; I >= 1; I-- {
					T = (ONE + ONE/XK) / RK
					SUM = SUM*EXPMX + T
					RK = RK - ONE
					XK = XK - X
				}
				RET = RET - SUM*EXPMX
			}
		}
	}
	return RET
}

// DEBYE2 calculates the Debye function of order 2, defined as
//   2* ∫ {0 to x} t^2/(exp(t)-1) dt / x^2
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func DEBYE2(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		TWENT4 = 24.0e0
		ONEHUN = 100.0e0
		DEBINF = 4.80822761263837714160e0
	)

	var I, NEXP, NTERMS int
	var EXPMX, RET, RK, SUM, T, X, XK, XLIM1, XLIM2, XLOW, XUPPER float64

	var ADEB2 = []float64{
		2.59438102325707702826e0,
		0.28633572045307198337e0,
		-0.1020626561580467129e-1,
		0.60491097753468435e-3,
		-0.4052576589502104e-4,
		0.286338263288107e-5,
		-0.20863943030651e-6,
		0.1552378758264e-7,
		-0.117312800866e-8,
		0.8973585888e-10,
		-0.693176137e-11,
		0.53980568e-12,
		-0.4232405e-13,
		0.333778e-14,
		-0.26455e-15,
		0.2106e-16,
		-0.168e-17,
		0.13e-18,
		-0.1e-19}

	X = XVALUE
	if X < ZERO {

		return ZERO
	}

	// Compute the machine-dependent constants.
	T = machine.D1MACH[1]
	XLIM1 = -math.Log(T)
	XLIM2 = math.Sqrt(DEBINF) / math.Sqrt(T)
	T = machine.D1MACH[3]
	XLOW = math.Sqrt(T * EIGHT)
	XUPPER = -math.Log(T + T)
	T = T / ONEHUN
	NTERMS = 18

	// Code for x <= 4.0

	if X <= FOUR {
		if X < XLOW {
			RET = ((X-EIGHT)*X + TWENT4) / TWENT4
		} else {
			T = ((X * X / EIGHT) - HALF) - HALF
			RET = utils.Cheval(NTERMS, ADEB2, T) - X/THREE
		}
	} else {
		// Code for x > 4.0
		if X > XLIM2 {
			RET = ZERO
		} else {
			RET = DEBINF / (X * X)
			if X < XLIM1 {
				EXPMX = math.Exp(-X)
				if X > XUPPER {
					SUM = ((X+TWO)*X + TWO) / (X * X)
				} else {
					SUM = ZERO
					RK = float64(int(XLIM1 / X))
					NEXP = int(RK)
					XK = RK * X
					for I = NEXP; I >= 1; I-- {
						T = (ONE + TWO/XK + TWO/(XK*XK)) / RK
						SUM = SUM*EXPMX + T
						RK = RK - ONE
						XK = XK - X
					}
				}
				RET = RET - TWO*SUM*EXPMX
			}
		}
	}
	return RET
}

// DEBYE3 calculates the Debye function of order 3, defined as
//    3 * ∫ {0 to x} t^3/(exp(t)-1) dt / x^3
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func DEBYE3(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		PT375  = 0.375e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		SIX    = 6.0e0
		SEVP5  = 7.5e0
		EIGHT  = 8.0e0
		TWENTY = 20.0e0
		ONEHUN = 100.0e0
		DEBINF = 0.51329911273421675946e-1
	)

	var I, NEXP, NTERMS int
	var EXPMX, RET, RK, SUM, T, X, XK, XKI, XLIM1, XLIM2, XLOW, XUPPER float64
	var ADEB3 = []float64{
		2.70773706832744094526e0,
		0.34006813521109175100e0,
		-0.1294515018444086863e-1,
		0.79637553801738164e-3,
		-0.5463600095908238e-4,
		0.392430195988049e-5,
		-0.28940328235386e-6,
		0.2173176139625e-7,
		-0.165420999498e-8,
		0.12727961892e-9,
		-0.987963459e-11,
		0.77250740e-12,
		-0.6077972e-13,
		0.480759e-14,
		-0.38204e-15,
		0.3048e-16,
		-0.244e-17,
		0.20e-18,
		-0.2e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	T = machine.D1MACH[1]
	XLIM1 = -math.Log(T)
	XK = ONE / THREE
	XKI = math.Pow(ONE/DEBINF, XK)
	RK = math.Pow(T, XK)
	XLIM2 = XKI / RK
	T = machine.D1MACH[3]
	XLOW = math.Sqrt(T * EIGHT)
	XUPPER = -math.Log(T + T)
	T = T / ONEHUN
	NTERMS = 18

	// Code for x <= 4.0

	if X <= FOUR {
		if X < XLOW {
			RET = ((X-SEVP5)*X + TWENTY) / TWENTY
		} else {
			T = ((X * X / EIGHT) - HALF) - HALF
			RET = utils.Cheval(NTERMS, ADEB3, T) - PT375*X
		}
	} else {

		// Code for x > 4.0

		if X > XLIM2 {
			RET = ZERO
		} else {
			RET = ONE / (DEBINF * X * X * X)
			if X < XLIM1 {
				EXPMX = math.Exp(-X)
				if X > XUPPER {
					SUM = (((X+THREE)*X+SIX)*X + SIX) / (X * X * X)
				} else {
					SUM = ZERO
					RK = float64(int(XLIM1 / X))
					NEXP = int(RK)
					XK = RK * X
					for I = NEXP; I >= 1; I-- {
						XKI = ONE / XK
						T = (((SIX*XKI+SIX)*XKI+THREE)*XKI + ONE) / RK
						SUM = SUM*EXPMX + T
						RK = RK - ONE
						XK = XK - X
					}
				}
				RET = RET - THREE*SUM*EXPMX
			}
		}
	}
	return RET
}

// DEBYE4 calculates the Debye function of order 3, defined as
//    4 * ∫ {0 to x} t^4/(exp(t)-1) dt / x^4
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func DEBYE4(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWOPT5 = 2.5e0

		FOUR = 4.0e0
		FIVE = 5.0e0

		EIGHT  = 8.0e0
		TWELVE = 12.0e0
		EIGHTN = 18.0e0
		TWENT4 = 24.0e0
		FORTY5 = 45.0e0
		ONEHUN = 100.0e0
		DEBINF = 99.54506449376351292781e0
	)

	var I, NEXP, NTERMS int
	var EXPMX, RET, RK, SUM, T, X, XK, XKI, XLIM1, XLIM2, XLOW, XUPPER float64

	var ADEB4 = []float64{
		2.78186941502052346008e0,
		0.37497678352689286364e0,
		-0.1494090739903158326e-1,
		0.94567981143704274e-3,
		-0.6613291613893255e-4,
		0.481563298214449e-5,
		-0.35880839587593e-6,
		0.2716011874160e-7,
		-0.208070991223e-8,
		0.16093838692e-9,
		-0.1254709791e-10,
		0.98472647e-12,
		-0.7772369e-13,
		0.616483e-14,
		-0.49107e-15,
		0.3927e-16,
		-0.315e-17,
		0.25e-18,
		-0.2e-19}

	X = XVALUE

	// Check XVALUE >= 0.0
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.

	T = machine.D1MACH[1]
	XLIM1 = -math.Log(T)
	RK = ONE / FOUR
	XK = math.Pow(DEBINF, RK)
	XKI = math.Pow(T, RK)
	XLIM2 = XK / XKI
	T = machine.D1MACH[3]
	XLOW = math.Sqrt(T * EIGHT)
	XUPPER = -math.Log(T + T)
	T = T / ONEHUN
	NTERMS = 18
	// Code for x <= 4.0

	if X <= FOUR {
		if X < XLOW {
			RET = ((TWOPT5*X-EIGHTN)*X + FORTY5) / FORTY5
		} else {
			T = ((X * X / EIGHT) - HALF) - HALF
			RET = utils.Cheval(NTERMS, ADEB4, T) - (X+X)/FIVE
		}
	} else {
		// Code for x > 4.0
		if X > XLIM2 {
			RET = ZERO
		} else {
			T = X * X
			RET = (DEBINF / T) / T
			if X < XLIM1 {
				EXPMX = math.Exp(-X)
				if X > XUPPER {
					SUM = ((((X+FOUR)*X+TWELVE)*X+TWENT4)*X + TWENT4) / (X * X * X * X)
				} else {
					SUM = ZERO
					RK = float64(int(XLIM1 / X))
					NEXP = int(RK)
					XK = RK * X
					for I = NEXP; I >= 1; I-- {
						XKI = ONE / XK
						T = ((((TWENT4*XKI+TWENT4)*XKI+TWELVE)*XKI+FOUR)*XKI + ONE) / RK
						SUM = SUM*EXPMX + T
						RK = RK - ONE
						XK = XK - X
					}
				}
				RET = RET - FOUR*SUM*EXPMX
			}
		}
	}
	return RET
}
