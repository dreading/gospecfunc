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

// I0INT computes the integral of the Bessel function I0(x) using the definition
//    ∫  0 to x I0(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
func I0INT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		THREE  = 3.0e0
		ATEEN  = 18.0e0
		THIRT6 = 36.0e0
		ONEHUN = 100.0e0
		LNR2PI = 0.91893853320467274178e0
	)

	var IND, NTERM1, NTERM2 int
	var RET, T, TEMP, X, XHIGH, XLOW float64

	var ARI01 = []float64{
		0.41227906926781516801e0,
		-0.34336345150081519562e0,
		0.22667588715751242585e0,
		-0.12608164718742260032e0,
		0.6012484628777990271e-1,
		-0.2480120462913358248e-1,
		0.892773389565563897e-2,
		-0.283253729936696605e-2,
		0.79891339041712994e-3,
		-0.20053933660964890e-3,
		0.4416816783014313e-4,
		-0.822377042246068e-5,
		0.120059794219015e-5,
		-0.11350865004889e-6,
		0.69606014466e-9,
		0.180622772836e-8,
		-0.26039481370e-9,
		-0.166188103e-11,
		0.510500232e-11,
		-0.41515879e-12,
		-0.7368138e-13,
		0.1279323e-13,
		0.103247e-14,
		-0.30379e-15,
		-0.1789e-16,
		0.673e-17,
		0.44e-18,
		-0.14e-18,
		-0.1e-19}

	var ARI0A = []float64{
		2.03739654571143287070e0,
		0.1917631647503310248e-1,
		0.49923334519288147e-3,
		0.2263187103659815e-4,
		0.158682108285561e-5,
		0.16507855636318e-6,
		0.2385058373640e-7,
		0.392985182304e-8,
		0.46042714199e-9,
		-0.7072558172e-10,
		-0.6747183961e-10,
		-0.2026962001e-10,
		-0.87320338e-12,
		0.175520014e-11,
		0.60383944e-12,
		-0.3977983e-13,
		-0.8049048e-13,
		-0.1158955e-13,
		0.827318e-14,
		0.282290e-14,
		-0.77667e-15,
		-0.48731e-15,
		0.7279e-16,
		0.7873e-16,
		-0.785e-17,
		-0.1281e-16,
		0.121e-17,
		0.214e-17,
		-0.27e-18,
		-0.36e-18,
		0.7e-19,
		0.6e-19,
		-0.2e-19,
		-0.1e-19}

	IND = 1
	X = XVALUE
	if XVALUE < ZERO {
		IND = -1
		X = -X
	}

	// Compute the machine-dependent constants.
	T = math.Log(machine.D1MACH[2])
	XHIGH = T + math.Log(T)*HALF - math.Log(HALF)

	// Error test

	if X > XHIGH {
		RET = math.Exp(XHIGH - LNR2PI - HALF*math.Log(XHIGH))
		if IND == -1 {
			RET = -RET
		}
		return RET
	}

	// Continue with machine-constants

	TEMP = machine.D1MACH[3]
	T = TEMP / ONEHUN
	if X <= ATEEN {
		NTERM1 = 28
		XLOW = math.Sqrt(THIRT6 * TEMP / THREE)
	} else {
		NTERM2 = 33
	}

	// Code for 0 <= |x| <= 18
	if X <= ATEEN {
		if X < XLOW {
			RET = X
		} else {
			T = (THREE*X - ATEEN) / (X + ATEEN)
			RET = X * math.Exp(X) * utils.Cheval(NTERM1, ARI01, T)
		}
	} else {
		// Code for |x| > 18
		T = (THIRT6/X - HALF) - HALF
		TEMP = X - HALF*math.Log(X) - LNR2PI + math.Log(utils.Cheval(NTERM2, ARI0A, T))
		RET = math.Exp(TEMP)
	}
	if IND == -1 {
		RET = -RET
	}
	return RET
}

// J0INT computes the integral of the Bessel function J0(x) using the definition
//    ∫  0 to x J0(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
func J0INT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1e0
		TWELVE = 12.0e0
		SIXTEN = 16.0e0
		ONEHUN = 100.0e0
		ONE28  = 128.0e0
		FIVE12 = 512.0e0
		RT2BPI = 0.79788456080286535588e0
		PIB411 = 201.0e0
		PIB412 = 256.0e0
		PIB42  = 0.24191339744830961566e-3
	)

	var IND, NTERM1, NTERM2, NTERM3 int
	var PIB41, RET, T, TEMP, X, XHIGH, XLOW, XMPI4 float64

	var ARJ01 = []float64{
		0.38179279321690173518e0,
		-0.21275636350505321870e0,
		0.16754213407215794187e0,
		-0.12853209772196398954e0,
		0.10114405455778847013e0,
		-0.9100795343201568859e-1,
		0.6401345264656873103e-1,
		-0.3066963029926754312e-1,
		0.1030836525325064201e-1,
		-0.255670650399956918e-2,
		0.48832755805798304e-3,
		-0.7424935126036077e-4,
		0.922260563730861e-5,
		-0.95522828307083e-6,
		0.8388355845986e-7,
		-0.633184488858e-8,
		0.41560504221e-9,
		-0.2395529307e-10,
		0.122286885e-11,
		-0.5569711e-13,
		0.227820e-14,
		-0.8417e-16,
		0.282e-17,
		-0.9e-19}

	var ARJ0A1 = []float64{
		1.24030133037518970827e0,
		-0.478125353632280693e-2,
		0.6613148891706678e-4,
		-0.186042740486349e-5,
		0.8362735565080e-7,
		-0.525857036731e-8,
		0.42606363251e-9,
		-0.4211761024e-10,
		0.488946426e-11,
		-0.64834929e-12,
		0.9617234e-13,
		-0.1570367e-13,
		0.278712e-14,
		-0.53222e-15,
		0.10844e-15,
		-0.2342e-16,
		0.533e-17,
		-0.127e-17,
		0.32e-18,
		-0.8e-19,
		0.2e-19,
		-0.1e-19}

	var ARJ0A2 = []float64{
		1.99616096301341675339e0,
		-0.190379819246668161e-2,
		0.1539710927044226e-4,
		-0.31145088328103e-6,
		0.1110850971321e-7,
		-0.58666787123e-9,
		0.4139926949e-10,
		-0.365398763e-11,
		0.38557568e-12,
		-0.4709800e-13,
		0.650220e-14,
		-0.99624e-15,
		0.16700e-15,
		-0.3028e-16,
		0.589e-17,
		-0.122e-17,
		0.27e-18,
		-0.6e-19,
		0.1e-19}

	X = XVALUE
	IND = 1
	if X < ZERO {
		X = -X
		IND = -1
	}

	// Compute the machine-dependent constants.
	TEMP = machine.D1MACH[3]
	XHIGH = ONE / TEMP

	// Error test
	if X > XHIGH {
		RET = ONE
		if IND == -1 {
			RET = -RET
		}
		return RET
	}

	// continue with constants
	T = TEMP / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 23
		XLOW = math.Sqrt(TWELVE * TEMP)
	} else {
		NTERM2 = 21
		NTERM3 = 18
	}

	// Code for 0 <= |x| <= 16
	if X <= SIXTEN {
		if X < XLOW {
			RET = X
		} else {
			T = X*X/ONE28 - ONE
			RET = X * utils.Cheval(NTERM1, ARJ01, T)
		}
	} else {
		// Code for |x| > 16
		T = FIVE12/(X*X) - ONE
		PIB41 = PIB411 / PIB412
		XMPI4 = (X - PIB41) - PIB42
		TEMP = math.Cos(XMPI4) * utils.Cheval(NTERM2, ARJ0A1, T) / X
		TEMP = TEMP - math.Sin(XMPI4)*utils.Cheval(NTERM3, ARJ0A2, T)
		RET = ONE - RT2BPI*TEMP/math.Sqrt(X)
	}
	if IND == -1 {
		RET = -RET
	}
	return RET
}

// K0INT computes the integral of the modified Bessel
// function K0(x) using the definition
//    ∫  0 to x K0(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
func K0INT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		SIX    = 6.0e0
		TWELVE = 12.0e0
		EIGHTN = 18.0e0
		ONEHUN = 100.0e0
		CONST1 = 1.11593151565841244881e0
		CONST2 = -0.11593151565841244881e0
		PIBY2  = 1.57079632679489661923e0
		RT2BPI = 0.79788456080286535588e0
	)

	var NTERM1, NTERM2, NTERM3 int
	var FVAL, RET, T, TEMP, X, XHIGH, XLOW float64

	var AK0IN1 = []float64{
		16.79702714464710959477e0,
		9.79134687676889407070e0,
		2.80501316044337939300e0,
		0.45615620531888502068e0,
		0.4716224457074760784e-1,
		0.335265148269698289e-2,
		0.17335181193874727e-3,
		0.679951889364702e-5,
		0.20900268359924e-6,
		0.516603846976e-8,
		0.10485708331e-9,
		0.177829320e-11,
		0.2556844e-13,
		0.31557e-15,
		0.338e-17,
		0.3e-19}

	var AK0IN2 = []float64{
		10.76266558227809174077e0,
		5.62333479849997511550e0,
		1.43543664879290867158e0,
		0.21250410143743896043e0,
		0.2036537393100009554e-1,
		0.136023584095623632e-2,
		0.6675388699209093e-4,
		0.250430035707337e-5,
		0.7406423741728e-7,
		0.176974704314e-8,
		0.3485775254e-10,
		0.57544785e-12,
		0.807481e-14,
		0.9747e-16,
		0.102e-17,
		0.1e-19}

	var AK0INA = []float64{
		1.91172065445060453895e0,
		-0.4183064565769581085e-1,
		0.213352508068147486e-2,
		-0.15859497284504181e-3,
		0.1497624699858351e-4,
		-0.167955955322241e-5,
		0.21495472478804e-6,
		-0.3058356654790e-7,
		0.474946413343e-8,
		-0.79424660432e-9,
		0.14156555325e-9,
		-0.2667825359e-10,
		0.528149717e-11,
		-0.109263199e-11,
		0.23518838e-12,
		-0.5247991e-13,
		0.1210191e-13,
		-0.287632e-14,
		0.70297e-15,
		-0.17631e-15,
		0.4530e-16,
		-0.1190e-16,
		0.319e-17,
		-0.87e-18,
		0.24e-18,
		-0.7e-19,
		0.2e-19,
		-0.1e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return ZERO
	}
	// Compute the machine-dependent constants.
	TEMP = machine.D1MACH[3]
	T = TEMP / ONEHUN
	if X <= SIX {
		NTERM1 = 15
		NTERM2 = 15
		XLOW = math.Sqrt(EIGHTN * TEMP)
	} else {
		NTERM3 = 27
		XHIGH = -math.Log(TEMP + TEMP)
	}
	// Code for 0 <= XVALUE <= 6
	if X <= SIX {
		if X < XLOW {
			FVAL = X
			if X > ZERO {
				FVAL = FVAL * (CONST1 - math.Log(X))
			}
			RET = FVAL
		} else {
			T = ((X*X)/EIGHTN - HALF) - HALF
			FVAL = (CONST2 + math.Log(X)) * utils.Cheval(NTERM2, AK0IN2, T)
			RET = X * (utils.Cheval(NTERM1, AK0IN1, T) - FVAL)
		}
		// Code for x > 6
	} else {
		FVAL = PIBY2
		if X < XHIGH {
			T = (TWELVE/X - HALF) - HALF
			TEMP = math.Exp(-X) * utils.Cheval(NTERM3, AK0INA, T)
			FVAL = FVAL - TEMP/(math.Sqrt(X)*RT2BPI)
		}
		RET = FVAL
	}
	return RET
}

// Y0INT computes the integral of the Bessel
// function Y0(x) using the definition
//    ∫  0 to x Y0(t) dt
// The code uses Chebyshev expansions with the coefficients
// given to an accuracy of 20 decimal places.
func Y0INT(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1.0e0
		NINE   = 9.0e0
		SIXTEN = 16.0e0

		ONEHUN = 100.0e0
		ONE28  = 128.0e0
		FIVE12 = 512.0e0
		RT2BPI = 0.79788456080286535588e0
		PIB411 = 201.0e0
		PIB412 = 256.0e0
		PIB42  = 0.24191339744830961566e-3
		TWOBPI = 0.63661977236758134308e0
		GAL2M1 = -1.11593151565841244881e0
		GAMLN2 = -0.11593151565841244881e0
	)

	var NTERM1, NTERM2, NTERM3, NTERM4 int
	var PIB41, RET, T, TEMP, X, XHIGH, XLOW, XMPI4 float64

	var ARJ01 = []float64{
		0.38179279321690173518e0,
		-0.21275636350505321870e0,
		0.16754213407215794187e0,
		-0.12853209772196398954e0,
		0.10114405455778847013e0,
		-0.9100795343201568859e-1,
		0.6401345264656873103e-1,
		-0.3066963029926754312e-1,
		0.1030836525325064201e-1,
		-0.255670650399956918e-2,
		0.48832755805798304e-3,
		-0.7424935126036077e-4,
		0.922260563730861e-5,
		-0.95522828307083e-6,
		0.8388355845986e-7,
		-0.633184488858e-8,
		0.41560504221e-9,
		-0.2395529307e-10,
		0.122286885e-11,
		-0.5569711e-13,
		0.227820e-14,
		-0.8417e-16,
		0.282e-17,
		-0.9e-19}

	var ARY01 = []float64{
		0.54492696302724365490e0,
		-0.14957323588684782157e0,
		0.11085634486254842337e0,
		-0.9495330018683777109e-1,
		0.6820817786991456963e-1,
		-0.10324653383368200408e0,
		0.10625703287534425491e0,
		-0.6258367679961681990e-1,
		0.2385645760338293285e-1,
		-0.644864913015404481e-2,
		0.131287082891002331e-2,
		-0.20988088174989640e-3,
		0.2716042484138347e-4,
		-0.291199114014694e-5,
		0.26344333093795e-6,
		-0.2041172069780e-7,
		0.137124781317e-8,
		-0.8070680792e-10,
		0.419883057e-11,
		-0.19459104e-12,
		0.808782e-14,
		-0.30329e-15,
		0.1032e-16,
		-0.32e-18,
		0.1e-19}

	var ARY0A1 = []float64{
		1.24030133037518970827e0,
		-0.478125353632280693e-2,
		0.6613148891706678e-4,
		-0.186042740486349e-5,
		0.8362735565080e-7,
		-0.525857036731e-8,
		0.42606363251e-9,
		-0.4211761024e-10,
		0.488946426e-11,
		-0.64834929e-12,
		0.9617234e-13,
		-0.1570367e-13,
		0.278712e-14,
		-0.53222e-15,
		0.10844e-15,
		-0.2342e-16,
		0.533e-17,
		-0.127e-17,
		0.32e-18,
		-0.8e-19,
		0.2e-19,
		-0.1e-19}

	var ARY0A2 = []float64{
		1.99616096301341675339e0,
		-0.190379819246668161e-2,
		0.1539710927044226e-4,
		-0.31145088328103e-6,
		0.1110850971321e-7,
		-0.58666787123e-9,
		0.4139926949e-10,
		-0.365398763e-11,
		0.38557568e-12,
		-0.4709800e-13,
		0.650220e-14,
		-0.99624e-15,
		0.16700e-15,
		-0.3028e-16,
		0.589e-17,
		-0.122e-17,
		0.27e-18,
		-0.6e-19,
		0.1e-19}

	X = XVALUE

	// First error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	TEMP = machine.D1MACH[3]
	XHIGH = ONE / TEMP

	// Second error test
	if X > XHIGH {
		return ZERO
	}

	// continue with machine constants
	T = TEMP / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 23
		NTERM2 = 24
		XLOW = math.Sqrt(NINE * TEMP)
	} else {
		NTERM3 = 21
		NTERM4 = 18
	}

	// Code for 0 <= x <= 16
	if X <= SIXTEN {
		if X < XLOW {
			if X == ZERO {
				RET = ZERO
			} else {
				RET = (math.Log(X) + GAL2M1) * TWOBPI * X
			}
		} else {
			T = X*X/ONE28 - ONE
			TEMP = (math.Log(X) + GAMLN2) * utils.Cheval(NTERM1, ARJ01, T)
			TEMP = TEMP - utils.Cheval(NTERM2, ARY01, T)
			RET = TWOBPI * X * TEMP
		}
	} else {
		// Code for x > 16
		T = FIVE12/(X*X) - ONE
		PIB41 = PIB411 / PIB412
		XMPI4 = (X - PIB41) - PIB42
		TEMP = math.Sin(XMPI4) * utils.Cheval(NTERM3, ARY0A1, T) / X
		TEMP = TEMP + math.Cos(XMPI4)*utils.Cheval(NTERM4, ARY0A2, T)
		RET = -RT2BPI * TEMP / math.Sqrt(X)
	}
	return RET
}
