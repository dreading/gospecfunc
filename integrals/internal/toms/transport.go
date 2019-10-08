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

// TRAN02 calculates the transport integral of order 2
//  ∫ 0 to x {t^2 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN02(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		THREE  = 3.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 2
		RNUMJN = 2.0e0
		VALINF = 0.32898681336964528729e1
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XLOW1, XHIGH1, XHIGH2, XHIGH3, XK, XK1 float64

	var ATRAN = []float64{
		1.67176044643453850301e0,
		-0.14773535994679448986e0,
		0.1482138199469363384e-1,
		-0.141953303263056126e-2,
		0.13065413244157083e-3,
		-0.1171557958675790e-4,
		0.103334984457557e-5,
		-0.9019113042227e-7,
		0.781771698331e-8,
		-0.67445656840e-9,
		0.5799463945e-10,
		-0.497476185e-11,
		0.42596097e-12,
		-0.3642189e-13,
		0.311086e-14,
		-0.26547e-15,
		0.2264e-16,
		-0.193e-17,
		0.16e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = ONE / (HALF * XK)
		XHIGH3 = math.Log(XK)
	}

	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW1 {
			RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
		} else {
			T = (((X * X) / EIGHT) - HALF) - HALF
			RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN03 calculates the transport integral of order 3
//  ∫ 0 to x {t^3 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN03(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 3
		RNUMJN = 3.0e0
		VALINF = 0.72123414189575657124e1
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.76201254324387200657e0,
		-0.10567438770505853250e0,
		0.1197780848196578097e-1,
		-0.121440152036983073e-2,
		0.11550997693928547e-3,
		-0.1058159921244229e-4,
		0.94746633853018e-6,
		-0.8362212128581e-7,
		0.731090992775e-8,
		-0.63505947788e-9,
		0.5491182819e-10,
		-0.473213954e-11,
		0.40676948e-12,
		-0.3489706e-13,
		0.298923e-14,
		-0.25574e-15,
		0.2186e-16,
		-0.187e-17,
		0.16e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XLOW2 = math.Sqrt(machine.D1MACH[1] / HALF)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}
	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN04 calculates the transport integral of order 4
//  ∫ 0 to x {t^4 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN04(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 4
		RNUMJN = 4.0e0
		VALINF = 0.25975757609067316596e2
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64
	var ATRAN = []float64{
		0.48075709946151105786e0,
		-0.8175378810321083956e-1,
		0.1002700665975162973e-1,
		-0.105993393598201507e-2,
		0.10345062450304053e-3,
		-0.964427054858991e-5,
		0.87455444085147e-6,
		-0.7793212079811e-7,
		0.686498861410e-8,
		-0.59995710764e-9,
		0.5213662413e-10,
		-0.451183819e-11,
		0.38921592e-12,
		-0.3349360e-13,
		0.287667e-14,
		-0.24668e-15,
		0.2113e-16,
		-0.181e-17,
		0.15e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}

	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN05 calculates the transport integral of order 5
//  ∫ 0 to x {t^5 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN05(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 5
		RNUMJN = 5.0e0
		VALINF = 0.12443133061720439116e3
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.34777777713391078928e0,
		-0.6645698897605042801e-1,
		0.861107265688330882e-2,
		-0.93966822237555384e-3,
		0.9363248060815134e-4,
		-0.885713193408328e-5,
		0.81191498914503e-6,
		-0.7295765423277e-7,
		0.646971455045e-8,
		-0.56849028255e-9,
		0.4962559787e-10,
		-0.431093996e-11,
		0.37310094e-12,
		-0.3219769e-13,
		0.277220e-14,
		-0.23824e-15,
		0.2044e-16,
		-0.175e-17,
		0.15e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}
	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN06 calculates the transport integral of order 6
//  ∫ 0 to x {t^6 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN06(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 6
		RNUMJN = 6.0e0
		VALINF = 0.73248700462880338059e3
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.27127335397840008227e0,
		-0.5588610553191453393e-1,
		0.753919513290083056e-2,
		-0.84351138579211219e-3,
		0.8549098079676702e-4,
		-0.818715493293098e-5,
		0.75754240427986e-6,
		-0.6857306541831e-7,
		0.611700376031e-8,
		-0.54012707024e-9,
		0.4734306435e-10,
		-0.412701055e-11,
		0.35825603e-12,
		-0.3099752e-13,
		0.267501e-14,
		-0.23036e-15,
		0.1980e-16,
		-0.170e-17,
		0.15e-18,
		-0.1e-19}

	X = XVALUE

	// Error test
	if X < ZERO {
		return ZERO
	}

	// Compute the machine-dependent constants.
	XK = machine.D1MACH[3]
	T = XK / ONEHUN
	if X <= FOUR {
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}

	// Code for x < =  4 .0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4 .0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN07 calculates the transport integral of order 7
//  ∫ 0 to x {t^7 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN07(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 7
		RNUMJN = 7.0e0
		VALINF = 0.50820803580048910473e4
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.22189250734010404423e0,
		-0.4816751061177993694e-1,
		0.670092448103153629e-2,
		-0.76495183443082557e-3,
		0.7863485592348690e-4,
		-0.761025180887504e-5,
		0.70991696299917e-6,
		-0.6468025624903e-7,
		0.580039233960e-8,
		-0.51443370149e-9,
		0.4525944183e-10,
		-0.395800363e-11,
		0.34453785e-12,
		-0.2988292e-13,
		0.258434e-14,
		-0.22297e-15,
		0.1920e-16,
		-0.165e-17,
		0.14e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}

	// Code for x <= 4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN08 calculates the transport integral of order 8
//  ∫ 0 to x {t^8 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN08(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 8
		RNUMJN = 8.0e0
		VALINF = 0.40484399001901115764e5
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.18750695774043719233e0,
		-0.4229527646093673337e-1,
		0.602814856929065592e-2,
		-0.69961054811814776e-3,
		0.7278482421298789e-4,
		-0.710846250050067e-5,
		0.66786706890115e-6,
		-0.6120157501844e-7,
		0.551465264474e-8,
		-0.49105307052e-9,
		0.4335000869e-10,
		-0.380218700e-11,
		0.33182369e-12,
		-0.2884512e-13,
		0.249958e-14,
		-0.21605e-15,
		0.1863e-16,
		-0.160e-17,
		0.14e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}
	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		//  Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}

// TRAN09 calculates the transport integral of order 9
//  ∫ 0 to x {t^9 exp(t)/[exp(t)-1]^2 } dt
// The code uses Chebyshev expansions with the coefficients
// given to 20 decimal places
func TRAN09(XVALUE float64) float64 {

	const (
		ZERO   = 0.0e0
		HALF   = 0.5e0
		ONE    = 1.0e0
		FOUR   = 4.0e0
		EIGHT  = 8.0e0
		ONEHUN = 100.0e0
		NUMJN  = 9
		RNUMJN = 9.0e0
		VALINF = 0.36360880558872871397e6
	)

	var K1, K2, NTERMS, NUMEXP int
	var RET, RK, SUMEXP, SUM2, T, X, XHIGH1, XHIGH2, XHIGH3, XK, XK1, XLOW1, XLOW2 float64

	var ATRAN = []float64{
		0.16224049991949846835e0,
		-0.3768351452195937773e-1,
		0.547669715917719770e-2,
		-0.64443945009449521e-3,
		0.6773645285280983e-4,
		-0.666813497582042e-5,
		0.63047560019047e-6,
		-0.5807478663611e-7,
		0.525551305123e-8,
		-0.46968861761e-9,
		0.4159395065e-10,
		-0.365808491e-11,
		0.32000794e-12,
		-0.2787651e-13,
		0.242017e-14,
		-0.20953e-15,
		0.1810e-16,
		-0.156e-17,
		0.13e-18,
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
		NTERMS = 19
		XLOW1 = math.Sqrt(EIGHT * XK)
		XK1 = RNUMJN - ONE
		XLOW2 = math.Pow(XK1*machine.D1MACH[1], ONE/XK1)
	} else {
		XHIGH1 = -math.Log(machine.D1MACH[4])
		XHIGH2 = RNUMJN / XK
		XHIGH3 = math.Log(XK)
	}
	// Code for x < =  4.0
	if X <= FOUR {
		if X < XLOW2 {
			RET = ZERO
		} else {
			if X < XLOW1 {
				RET = math.Pow(X, NUMJN-1) / (RNUMJN - ONE)
			} else {
				T = (((X * X) / EIGHT) - HALF) - HALF
				RET = math.Pow(X, NUMJN-1) * utils.Cheval(NTERMS, ATRAN, T)
			}
		}
	} else {
		// Code for x > 4.0
		if X > XHIGH2 {
			SUMEXP = ONE
		} else {
			if X <= XHIGH1 {
				NUMEXP = int(XHIGH1/X) + 1
				T = math.Exp(-X)
			} else {
				NUMEXP = 1
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
				for K2 = 1; K2 <= NUMJN; K2++ {
					SUM2 = SUM2*XK1*XK + ONE
					XK1 = XK1 + ONE
				}
				SUMEXP = SUMEXP*T + SUM2
				RK = RK - ONE
			}
		}
		T = RNUMJN*math.Log(X) - X + math.Log(SUMEXP)
		if T < XHIGH3 {
			RET = VALINF
		} else {
			RET = VALINF - math.Exp(T)
		}
	}
	return RET
}
