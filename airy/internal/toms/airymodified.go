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
	. "github.com/dreading/gospecfunc/machine"
	. "github.com/dreading/gospecfunc/misc"
	"math"
)

// AIRINT calculates the integral of the Airy function Ai,
//    ∫  0 to x Ai(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If XVALUE < 0.0, the function returns NaN
func AIRINT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		NINE   = 9.0e0
		FORTY1 = 41.0e0
		ONEHUN = 100.0e0
		NINHUN = 900.0e0
		FR996  = 4996.0e0
		PIBY4  = 0.78539816339744830962e0
		PITIM6 = 18.84955592153875943078e0
		RT2B3P = 0.46065886596178063902e0
		AIRZER = 0.35502805388781723926e0
	)

	var NTERM1, NTERM2, NTERM3, NTERM4, NTERM5 int

	var ARG, GVAL, HVAL, RET, T, TEMP, X, XHIGH1, XLOW1, XNEG1, Z float64

	var AAINT1 = []float64{
		0.37713517694683695526e0,
		-0.13318868432407947431e0,
		0.3152497374782884809e-1,
		-0.318543076436574077e-2,
		-0.87398764698621915e-3,
		0.46699497655396971e-3,
		-0.9544936738983692e-4,
		0.542705687156716e-5,
		0.239496406252188e-5,
		-0.75690270205649e-6,
		0.9050138584518e-7,
		0.320529456043e-8,
		-0.303825536444e-8,
		0.48900118596e-9,
		-0.1839820572e-10,
		-0.711247519e-11,
		0.151774419e-11,
		-0.10801922e-12,
		-0.963542e-14,
		0.313425e-14,
		-0.29446e-15,
		-0.477e-17,
		0.461e-17,
		-0.53e-18,
		0.1e-19,
		0.1e-19}

	var AAINT2 = []float64{
		1.92002524081984009769e0,
		-0.4220049417256287021e-1,
		-0.239457722965939223e-2,
		-0.19564070483352971e-3,
		-0.1547252891056112e-4,
		-0.140490186137889e-5,
		-0.12128014271367e-6,
		-0.1179186050192e-7,
		-0.104315578788e-8,
		-0.10908209293e-9,
		-0.929633045e-11,
		-0.110946520e-11,
		-0.7816483e-13,
		-0.1319661e-13,
		-0.36823e-15,
		-0.21505e-15,
		0.1238e-16,
		-0.557e-17,
		0.84e-18,
		-0.21e-18,
		0.4e-19,
		-0.1e-19}

	var AAINT3 = []float64{
		0.47985893264791052053e0,
		-0.19272375126169608863e0,
		0.2051154129525428189e-1,
		0.6332000070732488786e-1,
		-0.5093322261845754082e-1,
		0.1284424078661663016e-1,
		0.2760137088989479413e-1,
		-0.1547066673866649507e-1,
		-0.1496864655389316026e-1,
		0.336617614173574541e-2,
		0.530851163518892985e-2,
		0.41371226458555081e-3,
		-0.102490579926726266e-2,
		-0.32508221672025853e-3,
		0.8608660957169213e-4,
		0.6671367298120775e-4,
		0.449205999318095e-5,
		-0.670427230958249e-5,
		-0.196636570085009e-5,
		0.22229677407226e-6,
		0.22332222949137e-6,
		0.2803313766457e-7,
		-0.1155651663619e-7,
		-0.433069821736e-8,
		-0.6227777938e-10,
		0.26432664903e-9,
		0.5333881114e-10,
		-0.522957269e-11,
		-0.382229283e-11,
		-0.40958233e-12,
		0.11515622e-12,
		0.3875766e-13,
		0.140283e-14,
		-0.141526e-14,
		-0.28746e-15,
		0.923e-17,
		0.1224e-16,
		0.157e-17,
		-0.19e-18,
		-0.8e-19,
		-0.1e-19}

	var AAINT4 = []float64{
		1.99653305828522730048e0,
		-0.187541177605417759e-2,
		-0.15377536280305750e-3,
		-0.1283112967682349e-4,
		-0.108128481964162e-5,
		-0.9182131174057e-7,
		-0.784160590960e-8,
		-0.67292453878e-9,
		-0.5796325198e-10,
		-0.501040991e-11,
		-0.43420222e-12,
		-0.3774305e-13,
		-0.328473e-14,
		-0.28700e-15,
		-0.2502e-16,
		-0.220e-17,
		-0.19e-18,
		-0.2e-19}
	var AAINT5 = []float64{
		1.13024602034465716133e0,
		-0.464718064639872334e-2,
		-0.35137413382693203e-3,
		-0.2768117872545185e-4,
		-0.222057452558107e-5,
		-0.18089142365974e-6,
		-0.1487613383373e-7,
		-0.123515388168e-8,
		-0.10310104257e-9,
		-0.867493013e-11,
		-0.73080054e-12,
		-0.6223561e-13,
		-0.525128e-14,
		-0.45677e-15,
		-0.3748e-16,
		-0.356e-17,
		-0.23e-18,
		-0.4e-19}

	X = XVALUE
	// Compute the machine-dependent constants.
	Z = D1MACH[3]
	XLOW1 = TWO * Z
	ARG = D1MACH[4]
	XNEG1 = -1 / math.Pow(ARG, TWO/THREE)
	// Error test
	if X < XNEG1 {
		return -TWO / THREE
	}

	T = ARG / ONEHUN
	if X >= ZERO {
		NTERM1 = 25
		NTERM2 = 21
		XHIGH1 = math.Pow(-THREE*math.Log(Z)/TWO, TWO/THREE)
	} else {
		NTERM3 = 40
		NTERM4 = 17
		NTERM5 = 17
	}

	// Code for x >= 0

	if X >= ZERO {
		if X <= FOUR {
			if X < XLOW1 {
				RET = AIRZER * X
			} else {
				T = X/TWO - ONE
				RET = Cheval(NTERM1, AAINT1, T) * X
			}
		} else {
			if X > XHIGH1 {
				TEMP = ZERO
			} else {
				Z = (X + X) * math.Sqrt(X) / THREE
				TEMP = THREE * Z
				T = (FORTY1 - TEMP) / (NINE + TEMP)
				TEMP = math.Exp(-Z) * Cheval(NTERM2, AAINT2, T) / math.Sqrt(PITIM6*Z)
			}
			RET = ONE/THREE - TEMP
		}
	} else {
		// Code for x < 0

		if X >= -EIGHT {
			if X > -XLOW1 {
				RET = AIRZER * X
			} else {
				T = -X/FOUR - ONE
				RET = X * Cheval(NTERM3, AAINT3, T)
			}
		} else {
			Z = -(X + X) * math.Sqrt(-X) / THREE
			ARG = Z + PIBY4
			TEMP = NINE * Z * Z
			T = (FR996 - TEMP) / (NINHUN + TEMP)
			GVAL = Cheval(NTERM4, AAINT4, T)
			HVAL = Cheval(NTERM5, AAINT5, T)
			TEMP = GVAL*math.Cos(ARG) + HVAL*math.Sin(ARG)/Z
			RET = RT2B3P*TEMP/math.Sqrt(Z) - TWO/THREE
		}
	}
	return RET
}
