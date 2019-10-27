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

// AIRINT calculates the integral of the Airy function Ai,
//    ∫  0 to x Ai(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
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
	Z = machine.D1MACH[3]
	XLOW1 = TWO * Z
	ARG = machine.D1MACH[4]
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
				RET = utils.Cheval(NTERM1, AAINT1, T) * X
			}
		} else {
			if X > XHIGH1 {
				TEMP = ZERO
			} else {
				Z = (X + X) * math.Sqrt(X) / THREE
				TEMP = THREE * Z
				T = (FORTY1 - TEMP) / (NINE + TEMP)
				TEMP = math.Exp(-Z) * utils.Cheval(NTERM2, AAINT2, T) / math.Sqrt(PITIM6*Z)
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
				RET = X * utils.Cheval(NTERM3, AAINT3, T)
			}
		} else {
			Z = -(X + X) * math.Sqrt(-X) / THREE
			ARG = Z + PIBY4
			TEMP = NINE * Z * Z
			T = (FR996 - TEMP) / (NINHUN + TEMP)
			GVAL = utils.Cheval(NTERM4, AAINT4, T)
			HVAL = utils.Cheval(NTERM5, AAINT5, T)
			TEMP = GVAL*math.Cos(ARG) + HVAL*math.Sin(ARG)/Z
			RET = RT2B3P*TEMP/math.Sqrt(Z) - TWO/THREE
		}
	}
	return RET
}

// BIRINT calculates the integral of the Airy function Bi,
//    ∫  0 to x Bi(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
// If the function is too large and positive the correct value
// overflows and returns the Inf
func BIRINT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1.0e0
		ONEPT5 = 1.5e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		SEVEN  = 7.0e0
		EIGHT  = 8.0e0
		NINE   = 9.0e0
		SIXTEN = 16.0e0
		ONEHUN = 100.0e0
		NINHUN = 900.0e0
		THR644 = 3644.0e0
		PIBY4  = 0.78539816339744830962e0
		RT2B3P = 0.46065886596178063902e0
		BIRZER = 0.61492662744600073515e0
	)
	var NTERM1, NTERM2, NTERM3, NTERM4, NTERM5 int
	var ARG, F1, F2, RET, T, TEMP, X, XLOW1, XHIGH1, XMAX, XNEG1, Z float64

	var ABINT1 = []float64{
		0.38683352445038543350e0,
		-0.8823213550888908821e-1,
		0.21463937440355429239e0,
		-0.4205347375891315126e-1,
		0.5932422547496086771e-1,
		-0.840787081124270210e-2,
		0.871824772778487955e-2,
		-0.12191600199613455e-3,
		0.44024821786023234e-3,
		0.27894686666386678e-3,
		-0.7052804689785537e-4,
		0.5901080066770100e-4,
		-0.1370862587982142e-4,
		0.505962573749073e-5,
		-0.51598837766735e-6,
		0.397511312349e-8,
		0.9524985978055e-7,
		-0.3681435887321e-7,
		0.1248391688136e-7,
		-0.249097619137e-8,
		0.31775245551e-9,
		0.5434365270e-10,
		-0.4024566915e-10,
		0.1393855527e-10,
		-0.303817509e-11,
		0.40809511e-12,
		0.1634116e-13,
		-0.2683809e-13,
		0.896641e-14,
		-0.183089e-14,
		0.21333e-15,
		0.1108e-16,
		-0.1276e-16,
		0.363e-17,
		-0.62e-18,
		0.5e-19,
		0.1e-19}
	var ABINT2 = []float64{
		2.04122078602516135181e0,
		0.2124133918621221230e-1,
		0.66617599766706276e-3,
		0.3842047982808254e-4,
		0.362310366020439e-5,
		0.50351990115074e-6,
		0.7961648702253e-7,
		0.717808442336e-8,
		-0.267770159104e-8,
		-0.168489514699e-8,
		-0.36811757255e-9,
		0.4757128727e-10,
		0.5263621945e-10,
		0.778973500e-11,
		-0.460546143e-11,
		-0.183433736e-11,
		0.32191249e-12,
		0.29352060e-12,
		-0.1657935e-13,
		-0.4483808e-13,
		0.27907e-15,
		0.711921e-14,
		-0.1042e-16,
		-0.119591e-14,
		0.4606e-16,
		0.20884e-15,
		-0.2416e-16,
		-0.3638e-16,
		0.863e-17,
		0.591e-17,
		-0.256e-17,
		-0.77e-18,
		0.66e-18,
		0.3e-19,
		-0.15e-18,
		0.2e-19,
		0.3e-19,
		-0.1e-19}
	var ABINT3 = []float64{
		0.31076961598640349251e0,
		-0.27528845887452542718e0,
		0.17355965706136543928e0,
		-0.5544017909492843130e-1,
		-0.2251265478295950941e-1,
		0.4107347447812521894e-1,
		0.984761275464262480e-2,
		-0.1555618141666041932e-1,
		-0.560871870730279234e-2,
		0.246017783322230475e-2,
		0.165740392292336978e-2,
		-0.3277587501435402e-4,
		-0.24434680860514925e-3,
		-0.5035305196152321e-4,
		0.1630264722247854e-4,
		0.851914057780934e-5,
		0.29790363004664e-6,
		-0.64389707896401e-6,
		-0.15046988145803e-6,
		0.1587013535823e-7,
		0.1276766299622e-7,
		0.140578534199e-8,
		-0.46564739741e-9,
		-0.15682748791e-9,
		-0.403893560e-11,
		0.666708192e-11,
		0.128869380e-11,
		-0.6968663e-13,
		-0.6254319e-13,
		-0.718392e-14,
		0.115296e-14,
		0.42276e-15,
		0.2493e-16,
		-0.971e-17,
		-0.216e-17,
		-0.2e-19,
		0.6e-19,
		0.1e-19}
	var ABINT4 = []float64{
		1.99507959313352047614e0,
		-0.273736375970692738e-2,
		-0.30897113081285850e-3,
		-0.3550101982798577e-4,
		-0.412179271520133e-5,
		-0.48235892316833e-6,
		-0.5678730727927e-7,
		-0.671874810365e-8,
		-0.79811649857e-9,
		-0.9514271478e-10,
		-0.1137468966e-10,
		-0.136359969e-11,
		-0.16381418e-12,
		-0.1972575e-13,
		-0.237844e-14,
		-0.28752e-15,
		-0.3475e-16,
		-0.422e-17,
		-0.51e-18,
		-0.6e-19,
		-0.1e-19}
	var ABINT5 = []float64{
		1.12672081961782566017,
		-0.671405567525561198e-2,
		-0.69812918017832969e-3,
		-0.7561689886425276e-4,
		-0.834985574510207e-5,
		-0.93630298232480e-6,
		-0.10608556296250e-6,
		-0.1213128916741e-7,
		-0.139631129765e-8,
		-0.16178918054e-9,
		-0.1882307907e-10,
		-0.220272985e-11,
		-0.25816189e-12,
		-0.3047964e-13,
		-0.358370e-14,
		-0.42831e-15,
		-0.4993e-16,
		-0.617e-17,
		-0.68e-18,
		-0.10e-18,
		-0.1e-19}

	X = XVALUE

	// Compute the machine-dependent constants.
	T = machine.D1MACH[3]
	F2 = ONE + ONE
	XNEG1 = -ONE / math.Pow(T, F2/THREE)
	XMAX = machine.D1MACH[2]
	F1 = math.Log(XMAX)
	TEMP = F1 + math.Log(F1)/F2
	XHIGH1 = math.Pow(THREE*TEMP/F2, F2/THREE)

	// Error test
	if X > XHIGH1 {
		return XMAX
	}
	if X < XNEG1 {
		return 0
	}

	XLOW1 = F2 * T
	T = T / ONEHUN
	if X >= ZERO {
		NTERM1 = 36
		NTERM2 = 37
	} else {
		NTERM3 = 37
		NTERM4 = 20
		NTERM5 = 20
	}

	// Code for x >= 0.0
	if X >= ZERO {
		if X < XLOW1 {
			RET = BIRZER * X
		} else {
			if X <= EIGHT {
				T = X/FOUR - ONE
				RET = X * math.Exp(ONEPT5*X) * utils.Cheval(NTERM1, ABINT1, T)
			} else {
				T = SIXTEN*math.Sqrt(EIGHT/X)/X - ONE
				Z = (X + X) * math.Sqrt(X) / THREE
				TEMP = RT2B3P * utils.Cheval(NTERM2, ABINT2, T) / math.Sqrt(Z)
				TEMP = Z + math.Log(TEMP)
				RET = math.Exp(TEMP)
			}
		}
	} else {
		// Code for x < 0.0
		if X >= -SEVEN {
			if X > -XLOW1 {
				RET = BIRZER * X
			} else {
				T = -(X+X)/SEVEN - ONE
				RET = X * utils.Cheval(NTERM3, ABINT3, T)
			}
		} else {
			Z = -(X + X) * math.Sqrt(-X) / THREE
			ARG = Z + PIBY4
			TEMP = NINE * Z * Z
			T = (THR644 - TEMP) / (NINHUN + TEMP)
			F1 = utils.Cheval(NTERM4, ABINT4, T) * math.Sin(ARG)
			F2 = utils.Cheval(NTERM5, ABINT5, T) * math.Cos(ARG) / Z
			RET = (F2 - F1) * RT2B3P / math.Sqrt(Z)
		}
	}
	return RET
}
