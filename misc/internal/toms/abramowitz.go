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

// D1MACH PROVIDES DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS
var D1MACH = []float64{math.NaN(), 2.23e-308, 1.79e-308, 1.11e-16, 2.22e-16, 0.30103000998497009}

// ABRAM0 calculates the Abramowitz function of order 0,
//    ∫  0 to infinity exp( -t*t - x/t ) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If XVALUE < 0.0, the function prints a message and returns the value 0.0.
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
			RET = FVAL/ONERPI + X*(math.Log(X)*HVAL-GVAL)
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
	}
	return RET
}
