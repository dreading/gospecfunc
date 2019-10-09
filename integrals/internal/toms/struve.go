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

// STRVH0 calculates the Struve function of order 0, denoted H0(x),
//   (2/pi) ∫ {0 to pi/2} sin(x cos(t)) dt
// H0 also satisfies the second-order equation
//    x*D(Df) + Df + x*f = 2x/pi
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func STRVH0(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		EIGHT  = 8.0e0
		ELEVEN = 11.0e0
		TWENTY = 20.0e0
		SIXTP5 = 60.5e0
		TWO62  = 262.0e0
		THR2P5 = 302.5e0
		ONEHUN = 100.0e0
		PIBY4  = 0.78539816339744830962e0
		RT2BPI = 0.79788456080286535588e0
		TWOBPI = 0.63661977236758134308e0
	)

	var INDSGN, NTERM1, NTERM2, NTERM3, NTERM4 int
	var H0AS, RET, T, X, XHIGH, XLOW, XMP4, XSQ, Y0P, Y0Q, Y0VAL float64
	var ARRH0 = []float64{
		0.28696487399013225740e0,
		-0.25405332681618352305e0,
		0.20774026739323894439e0,
		-0.20364029560386585140e0,
		0.12888469086866186016e0,
		-0.4825632815622261202e-1,
		0.1168629347569001242e-1,
		-0.198118135642418416e-2,
		0.24899138512421286e-3,
		-0.2418827913785950e-4,
		0.187437547993431e-5,
		-0.11873346074362e-6,
		0.626984943346e-8,
		-0.28045546793e-9,
		0.1076941205e-10,
		-0.35904793e-12,
		0.1049447e-13,
		-0.27119e-15,
		0.624e-17,
		-0.13e-18}

	var ARRH0A = []float64{
		1.99291885751992305515e0,
		-0.384232668701456887e-2,
		-0.32871993712353050e-3,
		-0.2941181203703409e-4,
		-0.267315351987066e-5,
		-0.24681031075013e-6,
		-0.2295014861143e-7,
		-0.215682231833e-8,
		-0.20303506483e-9,
		-0.1934575509e-10,
		-0.182773144e-11,
		-0.17768424e-12,
		-0.1643296e-13,
		-0.171569e-14,
		-0.13368e-15,
		-0.2077e-16,
		0.2e-19,
		-0.55e-18,
		0.10e-18,
		-0.4e-19,
		0.1e-19}

	var AY0ASP = []float64{
		1.99944639402398271568e0,
		-0.28650778647031958e-3,
		-0.1005072797437620e-4,
		-0.35835941002463e-6,
		-0.1287965120531e-7,
		-0.46609486636e-9,
		-0.1693769454e-10,
		-0.61852269e-12,
		-0.2261841e-13,
		-0.83268e-15,
		-0.3042e-16,
		-0.115e-17,
		-0.4e-19}

	var AY0ASQ = []float64{
		1.99542681386828604092e0,
		-0.236013192867514472e-2,
		-0.7601538908502966e-4,
		-0.256108871456343e-5,
		-0.8750292185106e-7,
		-0.304304212159e-8,
		-0.10621428314e-9,
		-0.377371479e-11,
		-0.13213687e-12,
		-0.488621e-14,
		-0.15809e-15,
		-0.762e-17,
		-0.3e-19,
		-0.3e-19}

	// Start computation
	X = XVALUE
	INDSGN = 1
	if X < ZERO {
		X = -X
		INDSGN = -1
	}

	// Compute the machine-dependent constants.
	H0AS = machine.D1MACH[3]
	XHIGH = ONE / machine.D1MACH[4]

	// Error test
	if math.Abs(XVALUE) > XHIGH {
		return ZERO
	}

	// continue with machine constants
	T = H0AS / ONEHUN
	if X <= ELEVEN {
		NTERM1 = 19
		Y0P = math.Sqrt(H0AS)
		XLOW = Y0P + Y0P + Y0P
	} else {
		NTERM2 = 20
		NTERM3 = 12
		NTERM4 = 13
	}

	// Code for abs(x) <= 11
	if X <= ELEVEN {
		if X < XLOW {
			RET = TWOBPI * X
		} else {
			T = ((X*X)/SIXTP5 - HALF) - HALF
			RET = TWOBPI * X * utils.Cheval(NTERM1, ARRH0, T)
		}
	} else {
		// Code for abs(x) > 11
		XSQ = X * X
		T = (TWO62 - XSQ) / (TWENTY + XSQ)
		Y0P = utils.Cheval(NTERM3, AY0ASP, T)
		Y0Q = utils.Cheval(NTERM4, AY0ASQ, T) / (EIGHT * X)
		XMP4 = X - PIBY4
		Y0VAL = Y0P*math.Sin(XMP4) - Y0Q*math.Cos(XMP4)
		Y0VAL = Y0VAL * RT2BPI / math.Sqrt(X)
		T = (THR2P5 - XSQ) / (SIXTP5 + XSQ)
		H0AS = TWOBPI * utils.Cheval(NTERM2, ARRH0A, T) / X
		RET = Y0VAL + H0AS
	}
	if INDSGN == -1 {
		RET = -RET
	}
	return RET
}

// STRVH1 calculates the Struve function of order 1, denoted H1(x),
//   (2/pi) ∫ {0 to pi/2} sin( x cos(t))*sin^2 t dt
// H1 also satisfies the second-order differential equation
//                 x^2 * D^2 f  +  x * Df  +  (x^2 - 1)f  =  2x^2 / pi
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func STRVH1(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		EIGHT  = 8.0e0
		NINE   = 9.0e0
		ELEVEN = 11.0e0
		FIFTEN = 15.0e0
		TWENTY = 20.0e0
		ONEHUN = 100.0e0
		FORTP5 = 40.5e0
		ONE82  = 182.0e0
		TW02P5 = 202.5e0
		SIXTP5 = 60.5e0
		TWO6   = 262.0e0
		THR2P5 = 302.5e0
		RT2BPI = 0.79788456080286535588e0
		THPBY4 = 2.35619449019234492885e0
		TWOBPI = 0.63661977236758134308e0
	)

	var NTERM1, NTERM2, NTERM3, NTERM4 int
	var H1AS, RET, T, X, XHIGH, XLOW1, XLOW2, XM3P4, XSQ, Y1P, Y1Q, Y1VAL float64

	var ARRH1 = []float64{
		0.17319061083675439319e0,
		-0.12606917591352672005e0,
		0.7908576160495357500e-1,
		-0.3196493222321870820e-1,
		0.808040581404918834e-2,
		-0.136000820693074148e-2,
		0.16227148619889471e-3,
		-0.1442352451485929e-4,
		0.99219525734072e-6,
		-0.5441628049180e-7,
		0.243631662563e-8,
		-0.9077071338e-10,
		0.285926585e-11,
		-0.7716975e-13,
		0.180489e-14,
		-0.3694e-16,
		0.67e-18,
		-0.1e-19}

	var ARRH1A = []float64{2.01083504951473379407e0,
		0.592218610036099903e-2,
		0.55274322698414130e-3,
		0.5269873856311036e-4,
		0.506374522140969e-5,
		0.49028736420678e-6,
		0.4763540023525e-7,
		0.465258652283e-8,
		0.45465166081e-9,
		0.4472462193e-10,
		0.437308292e-11,
		0.43568368e-12,
		0.4182190e-13,
		0.441044e-14,
		0.36391e-15,
		0.5558e-16,
		-0.4e-19,
		0.163e-17,
		-0.34e-18,
		0.13e-18,
		-0.4e-19,
		0.1e-19}

	var AY1ASP = []float64{2.00135240045889396402e0,
		0.71104241596461938e-3,
		0.3665977028232449e-4,
		0.191301568657728e-5,
		0.10046911389777e-6,
		0.530401742538e-8,
		0.28100886176e-9,
		0.1493886051e-10,
		0.79578420e-12,
		0.4252363e-13,
		0.227195e-14,
		0.12216e-15,
		0.650e-17,
		0.36e-18,
		0.2e-19}
	var AY1ASQ = []float64{
		5.99065109477888189116e0,
		-0.489593262336579635e-2,
		-0.23238321307070626e-3,
		-0.1144734723857679e-4,
		-0.57169926189106e-6,
		-0.2895516716917e-7,
		-0.147513345636e-8,
		-0.7596537378e-10,
		-0.390658184e-11,
		-0.20464654e-12,
		-0.1042636e-13,
		-0.57702e-15,
		-0.2550e-16,
		-0.210e-17,
		0.2e-19,
		-0.2e-19}

	X = math.Abs(XVALUE)

	// Compute the machine-dependent constants.
	XHIGH = (HALF + HALF) / machine.D1MACH[4]

	// Error test

	if X > XHIGH {
		return ZERO
	}

	// continue with machine constants
	H1AS = machine.D1MACH[3]
	T = H1AS / ONEHUN
	if X <= NINE {
		NTERM1 = 17
		XLOW1 = HALF * math.Sqrt(machine.D1MACH[1])
		XLOW1 = XLOW1 + XLOW1 + XLOW1
		XLOW2 = math.Sqrt(FIFTEN * H1AS)
	} else {
		NTERM2 = 21
		NTERM3 = 14
		NTERM4 = 15
	}
	// Code for abs(x) <= 9
	if X <= NINE {
		if X < XLOW1 {
			RET = ZERO
		} else {
			XSQ = X * X
			if X < XLOW2 {
				RET = TWOBPI * XSQ
			} else {
				T = (XSQ/FORTP5 - HALF) - HALF
				RET = TWOBPI * XSQ * utils.Cheval(NTERM1, ARRH1, T)
			}
		}
	} else {
		// Code for abs(x) > 9
		XSQ = X * X
		T = (ONE82 - XSQ) / (TWENTY + XSQ)
		Y1P = utils.Cheval(NTERM3, AY1ASP, T)
		Y1Q = utils.Cheval(NTERM4, AY1ASQ, T) / (EIGHT * X)
		XM3P4 = X - THPBY4
		Y1VAL = Y1P*math.Sin(XM3P4) + Y1Q*math.Cos(XM3P4)
		Y1VAL = Y1VAL * RT2BPI / math.Sqrt(X)
		T = (TW02P5 - XSQ) / (FORTP5 + XSQ)
		H1AS = TWOBPI * utils.Cheval(NTERM2, ARRH1A, T)
		RET = Y1VAL + H1AS
	}
	return RET
}

// STRVL0 calculates the modified Struve function of order 0, denoted L0(x),
// defined as the solution of the second-order equation
//                x*D(Df) + Df - x*f  =  2x/pi
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
//
// If the value of |XVALUE| is too large, the result
// would cause an floating-pt overflow and the function returns NaN
func STRVL0(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		FOUR   = 4.0e0
		SIXTEN = 16.0e0
		TWENT4 = 24.0e0
		TWENT8 = 28.0e0
		TWO88  = 288.0e0
		ATEHUN = 800.0e0
		ONEHUN = 100.0e0
		LNR2PI = 0.91893853320467274178e0
		TWOBPI = 0.63661977236758134308e0
	)

	var INDSGN, NTERM1, NTERM2, NTERM3 int
	var CH1, CH2, RET, T, TEST, X, XHIGH1, XHIGH2, XLOW, XMAX, XSQ float64

	var ARL0 = []float64{
		0.42127458349979924863e0,
		-0.33859536391220612188e0,
		0.21898994812710716064e0,
		-0.12349482820713185712e0,
		0.6214209793866958440e-1,
		-0.2817806028109547545e-1,
		0.1157419676638091209e-1,
		-0.431658574306921179e-2,
		0.146142349907298329e-2,
		-0.44794211805461478e-3,
		0.12364746105943761e-3,
		-0.3049028334797044e-4,
		0.663941401521146e-5,
		-0.125538357703889e-5,
		0.20073446451228e-6,
		-0.2588260170637e-7,
		0.241143742758e-8,
		-0.10159674352e-9,
		-0.1202430736e-10,
		0.262906137e-11,
		-0.15313190e-12,
		-0.1574760e-13,
		0.315635e-14,
		-0.4096e-16,
		-0.3620e-16,
		0.239e-17,
		0.36e-18,
		-0.4e-19}

	var ARL0AS = []float64{
		2.00861308235605888600e0,
		0.403737966500438470e-2,
		-0.25199480286580267e-3,
		0.1605736682811176e-4,
		-0.103692182473444e-5,
		0.6765578876305e-7,
		-0.444999906756e-8,
		0.29468889228e-9,
		-0.1962180522e-10,
		0.131330306e-11,
		-0.8819190e-13,
		0.595376e-14,
		-0.40389e-15,
		0.2651e-16,
		-0.208e-17,
		0.11e-18}

	var AI0ML0 = []float64{
		2.00326510241160643125e0,
		0.195206851576492081e-2,
		0.38239523569908328e-3,
		0.7534280817054436e-4,
		0.1495957655897078e-4,
		0.299940531210557e-5,
		0.60769604822459e-6,
		0.12399495544506e-6,
		0.2523262552649e-7,
		0.504634857332e-8,
		0.97913236230e-9,
		0.18389115241e-9,
		0.3376309278e-10,
		0.611179703e-11,
		0.108472972e-11,
		0.18861271e-12,
		0.3280345e-13,
		0.565647e-14,
		0.93300e-15,
		0.15881e-15,
		0.2791e-16,
		0.389e-17,
		0.70e-18,
		0.16e-18}

	// Start computation
	X = XVALUE
	INDSGN = 1
	if X < ZERO {
		X = -X
		INDSGN = -1
	}

	// Compute the machine-dependent constants.
	TEST = machine.D1MACH[3]
	T = TEST / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 27
		XLOW = (ONE + TWO) * math.Sqrt(TEST)
	} else {
		NTERM2 = 15
		NTERM3 = 23
		XMAX = machine.D1MACH[2]
		XHIGH1 = math.Sqrt((TWENT8 + TWO) / TEST)
		XHIGH2 = TWENT8 / TEST
	}
	// Code for |xvalue| <= 16
	if X <= SIXTEN {
		if X < XLOW {
			RET = TWOBPI * X
		} else {
			T = (FOUR*X - TWENT4) / (X + TWENT4)
			RET = TWOBPI * X * utils.Cheval(NTERM1, ARL0, T) * math.Exp(X)
		}
	} else {
		// Code for |xvalue| > 16
		if X > XHIGH2 {
			CH1 = ONE
		} else {
			T = (X - TWENT8) / (FOUR - X)
			CH1 = utils.Cheval(NTERM2, ARL0AS, T)
		}
		if X > XHIGH1 {
			CH2 = ONE
		} else {
			XSQ = X * X
			T = (ATEHUN - XSQ) / (TWO88 + XSQ)
			CH2 = utils.Cheval(NTERM3, AI0ML0, T)
		}
		TEST = math.Log(CH1) - LNR2PI - math.Log(X)/TWO + X
		if TEST > math.Log(XMAX) {
			RET = math.NaN()
		} else {
			RET = math.Exp(TEST) - TWOBPI*CH2/X
		}
	}
	if INDSGN == -1 {
		RET = -RET
	}
	return RET
}

// STRVL1 calculates the modified Struve function of order 1, denoted L1(x),
// defined as the solution of the second-order equation
//               x^2*D(Df) + x*Df - (x^2+1)f = 2*x^2/pi
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
//
// If the value of |XVALUE| is too large, the result
// would cause an floating-pt overflow and the function returns NaN
func STRVL1(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1.0e0
		TWO    = 2.0e0
		FOUR   = 4.0e0
		SIXTEN = 16.0e0
		TWENT4 = 24.0e0
		THIRTY = 30.0e0
		TWO88  = 288.0e0
		ATEHUN = 800.0e0
		ONEHUN = 100.0e0
		LNR2PI = 0.91893853320467274178e0
		PI3BY2 = 4.71238898038468985769e0
		TWOBPI = 0.63661977236758134308e0
	)

	var NTERM1, NTERM2, NTERM3 int
	var CH1, CH2, RET, T, TEST, X, XHIGH1, XHIGH2, XLOW1, XLOW2, XMAX, XSQ float64

	var ARL1 = []float64{
		0.38996027351229538208e0,
		-0.33658096101975749366e0,
		0.23012467912501645616e0,
		-0.13121594007960832327e0,
		0.6425922289912846518e-1,
		-0.2750032950616635833e-1,
		0.1040234148637208871e-1,
		-0.350532294936388080e-2,
		0.105748498421439717e-2,
		-0.28609426403666558e-3,
		0.6925708785942208e-4,
		-0.1489693951122717e-4,
		0.281035582597128e-5,
		-0.45503879297776e-6,
		0.6090171561770e-7,
		-0.623543724808e-8,
		0.38430012067e-9,
		0.790543916e-11,
		-0.489824083e-11,
		0.46356884e-12,
		0.684205e-14,
		-0.569748e-14,
		0.35324e-15,
		0.4244e-16,
		-0.644e-17,
		-0.21e-18,
		0.9e-19}

	var ARL1AS = []float64{
		1.97540378441652356868e0,
		-0.1195130555088294181e-1,
		0.33639485269196046e-3,
		-0.1009115655481549e-4,
		0.30638951321998e-6,
		-0.953704370396e-8,
		0.29524735558e-9,
		-0.951078318e-11,
		0.28203667e-12,
		-0.1134175e-13,
		0.147e-17,
		-0.6232e-16,
		-0.751e-17,
		-0.17e-18,
		0.51e-18,
		0.23e-18,
		0.5e-19}

	var AI1ML1 = []float64{
		1.99679361896789136501e0,
		-0.190663261409686132e-2,
		-0.36094622410174481e-3,
		-0.6841847304599820e-4,
		-0.1299008228509426e-4,
		-0.247152188705765e-5,
		-0.47147839691972e-6,
		-0.9020819982592e-7,
		-0.1730458637504e-7,
		-0.332323670159e-8,
		-0.63736421735e-9,
		-0.12180239756e-9,
		-0.2317346832e-10,
		-0.439068833e-11,
		-0.82847110e-12,
		-0.15562249e-12,
		-0.2913112e-13,
		-0.543965e-14,
		-0.101177e-14,
		-0.18767e-15,
		-0.3484e-16,
		-0.643e-17,
		-0.118e-17,
		-0.22e-18,
		-0.4e-19,
		-0.1e-19}

	// START CALCULATION
	X = math.Abs(XVALUE)
	// Compute the machine-dependent constants.

	TEST = machine.D1MACH[3]
	T = TEST / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 26
		XLOW1 = math.Sqrt(THIRTY * TEST / TWO)
		XLOW2 = math.Sqrt((FOUR + ONE) * machine.D1MACH[1])
	} else {
		NTERM2 = 16
		NTERM3 = 25
		XMAX = machine.D1MACH[2]
		XHIGH2 = THIRTY / TEST
		XHIGH1 = math.Sqrt(XHIGH2)
	}

	// CODE FOR |XVALUE| <= 16
	if X <= SIXTEN {
		if X <= XLOW2 {
			RET = ZERO
		} else {
			XSQ = X * X
			if X < XLOW1 {
				RET = XSQ / PI3BY2
			} else {
				T = (FOUR*X - TWENT4) / (X + TWENT4)
				RET = XSQ * utils.Cheval(NTERM1, ARL1, T) * math.Exp(X) / PI3BY2
			}
		}
	} else {
		// CODE FOR |XVALUE| > 16
		if X > XHIGH2 {
			CH1 = ONE
		} else {
			T = (X - THIRTY) / (TWO - X)
			CH1 = utils.Cheval(NTERM2, ARL1AS, T)
		}
		if X > XHIGH1 {
			CH2 = ONE
		} else {
			XSQ = X * X
			T = (ATEHUN - XSQ) / (TWO88 + XSQ)
			CH2 = utils.Cheval(NTERM3, AI1ML1, T)
		}
		TEST = math.Log(CH1) - LNR2PI - math.Log(X)/TWO + X
		if TEST > math.Log(XMAX) {
			RET = math.NaN()
		} else {
			RET = math.Exp(TEST) - TWOBPI*CH2
		}
	}
	return RET
}

// I0ML0 calculates I0(x) - L0(x)
// where I0(x) is the modified Bessel function of the first kind of
// order 0, and L0(x) is the modified Struve function of order 0.
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func I0ML0(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		ONE    = 1.0e0
		SIX    = 6.0e0
		SIXTEN = 16.0e0
		FORTY  = 40.0e0
		TWO88  = 288.0e0
		ATEHUN = 800.0e0
		ONEHUN = 100.0e0
		TWOBPI = 0.63661977236758134308e0
	)

	var NTERM1, NTERM2 int
	var RET, T, X, XHIGH, XLOW, XSQ float64

	var AI0L0 = []float64{
		0.52468736791485599138e0,
		-0.35612460699650586196e0,
		0.20487202864009927687e0,
		-0.10418640520402693629e0,
		0.4634211095548429228e-1,
		-0.1790587192403498630e-1,
		0.597968695481143177e-2,
		-0.171777547693565429e-2,
		0.42204654469171422e-3,
		-0.8796178522094125e-4,
		0.1535434234869223e-4,
		-0.219780769584743e-5,
		0.24820683936666e-6,
		-0.2032706035607e-7,
		0.90984198421e-9,
		0.2561793929e-10,
		-0.710609790e-11,
		0.32716960e-12,
		0.2300215e-13,
		-0.292109e-14,
		-0.3566e-16,
		0.1832e-16,
		-0.10e-18,
		-0.11e-18}

	var AI0L0A = []float64{
		2.00326510241160643125e0,
		0.195206851576492081e-2,
		0.38239523569908328e-3,
		0.7534280817054436e-4,
		0.1495957655897078e-4,
		0.299940531210557e-5,
		0.60769604822459e-6,
		0.12399495544506e-6,
		0.2523262552649e-7,
		0.504634857332e-8,
		0.97913236230e-9,
		0.18389115241e-9,
		0.3376309278e-10,
		0.611179703e-11,
		0.108472972e-11,
		0.18861271e-12,
		0.3280345e-13,
		0.565647e-14,
		0.93300e-15,
		0.15881e-15,
		0.2791e-16,
		0.389e-17,
		0.70e-18,
		0.16e-18}

	X = XVALUE
	// Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.

	XSQ = machine.D1MACH[3]
	T = XSQ / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 23
		XLOW = XSQ
	} else {
		NTERM2 = 23
		XHIGH = math.Sqrt(ATEHUN / XSQ)
	}
	// Code for x <= 16
	if X <= SIXTEN {
		if X < XLOW {
			return ONE
		}
		T = (SIX*X - FORTY) / (X + FORTY)
		return utils.Cheval(NTERM1, AI0L0, T)
	} else {
		// Code for x > 16
		if X > XHIGH {
			RET = TWOBPI / X
		} else {
			XSQ = X * X
			T = (ATEHUN - XSQ) / (TWO88 + XSQ)
			RET = utils.Cheval(NTERM2, AI0L0A, T) * TWOBPI / X
		}
	}
	return RET
}

// I1ML1 calculates I1(x) - L1(x)
// where I1(x) is the modified Bessel function of the first kind of
// order 1, and L1(x) is the modified Struve function of order 1.
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func I1ML1(XVALUE float64) float64 {

	const (
		ZERO = 0.0e0
		ONE  = 1.0e0
		TWO  = 2.0e0

		SIX    = 6.0e0
		SIXTEN = 16.0e0
		FORTY  = 40.0e0
		TWO88  = 288.0e0
		ATEHUN = 800.0e0
		ONEHUN = 100.0e0
		TWOBPI = 0.63661977236758134308e0
	)

	var NTERM1, NTERM2 int
	var RET, T, X, XHIGH, XLOW, XSQ float64

	var AI1L1 = []float64{
		0.67536369062350576137e0,
		-0.38134971097266559040e0,
		0.17452170775133943559e0,
		-0.7062105887235025061e-1,
		0.2517341413558803702e-1,
		-0.787098561606423321e-2,
		0.214814368651922006e-2,
		-0.50862199717906236e-3,
		0.10362608280442330e-3,
		-0.1795447212057247e-4,
		0.259788274515414e-5,
		-0.30442406324667e-6,
		0.2720239894766e-7,
		-0.158126144190e-8,
		0.1816209172e-10,
		0.647967659e-11,
		-0.54113290e-12,
		-0.308311e-14,
		0.305638e-14,
		-0.9717e-16,
		-0.1422e-16,
		0.84e-18,
		0.7e-19,
		-0.1e-19}

	var AI1L1A = []float64{
		1.99679361896789136501e0,
		-0.190663261409686132e-2,
		-0.36094622410174481e-3,
		-0.6841847304599820e-4,
		-0.1299008228509426e-4,
		-0.247152188705765e-5,
		-0.47147839691972e-6,
		-0.9020819982592e-7,
		-0.1730458637504e-7,
		-0.332323670159e-8,
		-0.63736421735e-9,
		-0.12180239756e-9,
		-0.2317346832e-10,
		-0.439068833e-11,
		-0.82847110e-12,
		-0.15562249e-12,
		-0.2913112e-13,
		-0.543965e-14,
		-0.101177e-14,
		-0.18767e-15,
		-0.3484e-16,
		-0.643e-17,
		-0.118e-17,
		-0.22e-18,
		-0.4e-19,
		-0.1e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return ZERO
	}
	// Compute the machine-dependent constants.
	XSQ = machine.D1MACH[3]
	T = XSQ / ONEHUN
	if X <= SIXTEN {
		NTERM1 = 23
		XLOW = XSQ + XSQ
	} else {
		NTERM2 = 25
		XHIGH = math.Sqrt(ATEHUN / XSQ)
	}

	// Code for x <= 16

	if X <= SIXTEN {
		if X < XLOW {
			return X / TWO
		}
		T = (SIX*X - FORTY) / (X + FORTY)
		return utils.Cheval(NTERM1, AI1L1, T) * X / TWO

	} else {
		// Code for x > 16
		if X > XHIGH {
			RET = TWOBPI
		} else {
			XSQ = X * X
			T = (ATEHUN - XSQ) / (TWO88 + XSQ)
			RET = utils.Cheval(NTERM2, AI1L1A, T) * TWOBPI
		}
	}
	return RET
}
