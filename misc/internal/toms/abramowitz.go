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
	. "github.com/dreading/gospecfunc/machine"
)
 
// ABRAM0 calculates the Abramowitz function of order 0,
//    ∫  0 to infinity exp( -t*t - x/t ) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If XVALUE < 0.0, the function returns NaN
func ABRAM0(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		SIX    = 6.0e0
		ONEHUN = 100.0e0
		RT3BPI = 0.97720502380583984317
		RTPIB2 = 0.88622692545275801365e0
		GVAL0  = 0.13417650264770070909e0
		ONERPI = 0.56418958354775628695e0
	)

	var NTERMA, NTERMF, NTERMG, NTERMH int
	var ASLN, ASVAL, FVAL, GVAL, HVAL, LNXMIN, T, V, X, XLOW1, RET float64

	var AB0F = []float64{-0.68121927093549469816e0,
		-0.78867919816149252495e0,
		0.5121581776818819543e-1,
		-0.71092352894541296e-3,
		0.368681808504287e-5,
		-0.917832337237e-8,
		0.1270202563e-10,
		-0.1076888e-13,
		0.599e-17}

	var AB0G = []float64{-0.60506039430868273190e0,
		-0.41950398163201779803e0,
		0.1703265125190370333e-1,
		-0.16938917842491397e-3,
		0.67638089519710e-6,
		-0.135723636255e-8,
		0.156297065e-11,
		-0.112887e-14,
		0.55e-18}

	var AB0H = []float64{1.38202655230574989705e0,
		-0.30097929073974904355e0,
		0.794288809364887241e-2,
		-0.6431910276847563e-4,
		0.22549830684374e-6,
		-0.41220966195e-9,
		0.44185282e-12,
		-0.30123e-15,
		0.14e-18}

	var AB0AS = []float64{1.97755499723693067407e0,
		-0.1046024792004819485e-1,
		0.69680790253625366e-3,
		-0.5898298299996599e-4,
		0.577164455305320e-5,
		-0.61523013365756e-6,
		0.6785396884767e-7,
		-0.723062537907e-8,
		0.63306627365e-9,
		-0.989453793e-11,
		-0.1681980530e-10,
		0.673799551e-11,
		-0.200997939e-11,
		0.54055903e-12,
		-0.13816679e-12,
		0.3422205e-13,
		-0.826686e-14,
		0.194566e-14,
		-0.44268e-15,
		0.9562e-16,
		-0.1883e-16,
		0.301e-17,
		-0.19e-18,
		-0.14e-18,
		0.11e-18,
		-0.4e-19,
		0.2e-19,
		-0.1e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return math.NaN()
	}

	// Compute the machine-dependent constants.
	T = D1MACH[4] / ONEHUN
	if X <= TWO {
		NTERMF = 8
		NTERMG = 8
		NTERMH = 8
		XLOW1 = math.Sqrt(TWO * D1MACH[3])
	} else {
		NTERMA = 27
		LNXMIN = math.Log(D1MACH[1])
	}

	// Code for 0 <= XVALUE <= 2
	if X <= TWO {
		if X == ZERO {
			return RTPIB2
		}
		if X < XLOW1 {
			return RTPIB2 + X*(math.Log(X)-GVAL0)
		} else {
			T = (X*X/TWO - HALF) - HALF
			FVAL = CHEVAL(NTERMF, AB0F, T)
			GVAL = CHEVAL(NTERMG, AB0G, T)
			HVAL = CHEVAL(NTERMH, AB0H, T)
			return FVAL/ONERPI + X*(math.Log(X)*HVAL-GVAL)
		}
	} else {
		// Code for XVALUE > 2
		V = THREE * math.Pow(X/TWO, TWO/THREE)
		T = (SIX/V - HALF) - HALF
		ASVAL = CHEVAL(NTERMA, AB0AS, T)
		ASLN = math.Log(ASVAL/RT3BPI) - V
		if ASLN < LNXMIN {
			RET = ZERO
		} else {
			RET = math.Exp(ASLN)
		}
		return RET
	}
}

// ABRAM1 calculates the Abramowitz function of order 1,
//    ∫  0 to infinity t * exp( -t*t - x/t ) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If XVALUE < 0.0, the function returns NaN
func ABRAM1(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		SIX    = 6.0e0
		ONEHUN = 100.0e0
		RT3BPI = 0.97720502380583984317
		ONERPI = 0.56418958354775628695
	)

	var NTERMA, NTERMF, NTERMG, NTERMH int
	var ASLN, ASVAL, FVAL, GVAL, HVAL, LNXMIN, RET, T, V, X, XLOW, XLOW1 float64

	var AB1F = []float64{1.47285192577978807369e0,
		0.10903497570168956257e0,
		-0.12430675360056569753e0,
		0.306197946853493315e-2,
		-0.2218410323076511e-4,
		0.6989978834451e-7,
		-0.11597076444e-9,
		0.11389776e-12,
		-0.7173e-16,
		0.3e-19}

	var AB1G = []float64{0.39791277949054503528e0,
		-0.29045285226454720849e0,
		0.1048784695465363504e-1,
		-0.10249869522691336e-3,
		0.41150279399110e-6,
		-0.83652638940e-9,
		0.97862595e-12,
		-0.71868e-15,
		0.35e-18}

	var AB1H = []float64{0.84150292152274947030e0,
		-0.7790050698774143395e-1,
		0.133992455878390993e-2,
		-0.808503907152788e-5,
		0.2261858281728e-7,
		-0.3441395838e-10,
		0.3159858e-13,
		-0.1884e-16,
		0.1e-19}

	var AB1AS = []float64{2.13013643429065549448e0,
		0.6371526795218539933e-1,
		-0.129334917477510647e-2,
		0.5678328753228265e-4,
		-0.279434939177646e-5,
		0.5600214736787e-7,
		0.2392009242798e-7,
		-0.750984865009e-8,
		0.173015330776e-8,
		-0.36648877955e-9,
		0.7520758307e-10,
		-0.1517990208e-10,
		0.301713710e-11,
		-0.58596718e-12,
		0.10914455e-12,
		-0.1870536e-13,
		0.262542e-14,
		-0.14627e-15,
		-0.9500e-16,
		0.5873e-16,
		-0.2420e-16,
		0.868e-17,
		-0.290e-17,
		0.93e-18,
		-0.29e-18,
		0.9e-19,
		-0.3e-19,
		0.1e-19}

	X = XVALUE
	// Error test
	if X < ZERO {
		return math.NaN()
	}

	// Compute the machine-dependent constants.
	T = D1MACH[4] / ONEHUN
	if X <= TWO {
		NTERMF = 9
		NTERMG = 8
		NTERMH = 8
		T = D1MACH[3]
		XLOW1 = math.Sqrt(TWO * T)
		XLOW = T / TWO
	} else {
		NTERMA = 27
		LNXMIN = math.Log(D1MACH[1])
	}

	// Code for 0 <= XVALUE <= 2
	if X <= TWO {
		if X == ZERO {
			return HALF
		}
		if X < XLOW1 {
			if X < XLOW {
				RET = HALF
			} else {
				RET = (ONE - X/ONERPI - X*X*math.Log(X)) * HALF
			}
			return RET
		} else {
			T = (X*X/TWO - HALF) - HALF
			FVAL = CHEVAL(NTERMF, AB1F, T)
			GVAL = CHEVAL(NTERMG, AB1G, T)
			HVAL = CHEVAL(NTERMH, AB1H, T)
			return FVAL - X*(GVAL/ONERPI+X*math.Log(X)*HVAL)
		}
	} else {
		// Code for XVALUE > 2
		V = THREE * math.Pow(X/TWO, TWO/THREE)
		T = (SIX/V - HALF) - HALF
		ASVAL = CHEVAL(NTERMA, AB1AS, T)
		ASLN = math.Log(ASVAL*math.Sqrt(V/THREE)/RT3BPI) - V
		if ASLN < LNXMIN {
			RET = ZERO
		} else {
			RET = math.Exp(ASLN)
		}
		return RET
	}
}

// ABRAM2 calculates the Abramowitz function of order 2,
//    ∫  0 to infinity (t^2) * exp( -t*t - x/t ) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If XVALUE < 0.0, the function returns NaN
func ABRAM2(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		SIX    = 6.0e0
		ONEHUN = 100.0e0
		RT3BPI = 0.97720502380583984317
		RTPIB4 = 0.44311346272637900682
		ONERPI = 0.56418958354775628695
	)

	var NTERMA, NTERMF, NTERMG, NTERMH int
	var ASLN, ASVAL, FVAL, GVAL, HVAL, LNXMIN, RET, T, V, X, XLOW, XLOW1 float64

	var AB2F = []float64{1.03612162804243713846e0,
		0.19371246626794570012e0,
		-0.7258758839233007378e-1,
		0.174790590864327399e-2,
		-0.1281223233756549e-4,
		0.4115018153651e-7,
		-0.6971047256e-10,
		0.6990183e-13,
		-0.4492e-16,
		0.2e-19}
	var AB2G = []float64{1.46290157198630741150e0,
		0.20189466883154014317e0,
		-0.2908292087997129022e-1,
		0.47061049035270050e-3,
		-0.257922080359333e-5,
		0.656133712946e-8,
		-0.914110203e-11,
		0.774276e-14,
		-0.429e-17}
	var AB2H = []float64{0.30117225010910488881e0,
		-0.1588667818317623783e-1,
		0.19295936935584526e-3,
		-0.90199587849300e-6,
		0.206105041837e-8,
		-0.265111806e-11,
		0.210864e-14,
		-0.111e-17}

	var AB2AS = []float64{
		2.46492325304334856893e0,
		0.23142797422248905432e0,
		-0.94068173010085773e-3,
		0.8290270038089733e-4,
		-0.883894704245866e-5,
		0.106638543567985e-5,
		-0.13991128538529e-6,
		0.1939793208445e-7,
		-0.277049938375e-8,
		0.39590687186e-9,
		-0.5408354342e-10,
		0.635546076e-11,
		-0.38461613e-12,
		-0.11696067e-12,
		0.6896671e-13,
		-0.2503113e-13,
		0.785586e-14,
		-0.230334e-14,
		0.64914e-15,
		-0.17797e-15,
		0.4766e-16,
		-0.1246e-16,
		0.316e-17,
		-0.77e-18,
		0.18e-18,
		-0.4e-19,
		0.1e-19}

	X = XVALUE
	// Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	T = D1MACH[4] / ONEHUN
	if X <= TWO {
		NTERMF = 9
		NTERMG = 8
		NTERMH = 7
		XLOW = D1MACH[3]
		XLOW1 = math.Sqrt(TWO * XLOW)
	} else {
		NTERMA = 26
		LNXMIN = math.Log(D1MACH[1])
	}
	// Code for 0 <= XVALUE <= 2
	if X <= TWO {
		if X == ZERO {
			return RTPIB4
		}
		if X < XLOW1 {
			if X < XLOW {
				RET = RTPIB4
			} else {
				RET = RTPIB4 - HALF*X + X*X*X*math.Log(X)/SIX
			}
			return RET
		} else {
			T = (X*X/TWO - HALF) - HALF
			FVAL = CHEVAL(NTERMF, AB2F, T)
			GVAL = CHEVAL(NTERMG, AB2G, T)
			HVAL = CHEVAL(NTERMH, AB2H, T)
			RET = FVAL/ONERPI + X*(X*X*math.Log(X)*HVAL-GVAL)
			return RET
		}
	} else {
		// Code for XVALUE > 2
		V = THREE * math.Pow(X/TWO, TWO/THREE)
		T = (SIX/V - HALF) - HALF
		ASVAL = CHEVAL(NTERMA, AB2AS, T)
		ASLN = math.Log(ASVAL/RT3BPI) + math.Log(V/THREE) - V
		if ASLN < LNXMIN {
			RET = ZERO
		} else {
			RET = math.Exp(ASLN)
		}
		return RET
	}
}
