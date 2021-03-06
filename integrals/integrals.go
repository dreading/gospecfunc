// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc

import (
	"github.com/dreading/gospecfunc/integrals/internal/toms"
)

// Abramowitz computes the  Abramowitz function
//   ∫ 0 to ∞ of t^order exp( -t*t - x/t ) dt
func Abramowitz(order int, x float64) float64 {
	switch order {
	case 0:
		return toms.ABRAM0(x)
	case 1:
		return toms.ABRAM1(x)
	case 2:
		return toms.ABRAM2(x)
	default:
		panic("Order must be 0, 1 or 2")
	}
}

// Clausen calculates the Clausen's integral,
//   ∫ 0 to x of (-ln(2*sin(t/2))) dt
func Clausen(x float64) float64 {
	return toms.CLAUSN(x)
}

// Debye calculates the Debye function, defined as
//   order*∫ 0 to x of t^order/(exp(t)-1) dt] / x^order
func Debye(order int, x float64) float64 {
	switch order {
	case 1:
		return toms.DEBYE1(x)
	case 2:
		return toms.DEBYE2(x)
	case 3:
		return toms.DEBYE3(x)
	case 4:
		return toms.DEBYE4(x)
	default:
		panic("Order must be between 1 and 4")
	}
}

// Goodst calculates the function defined as
//   ∫ 0 to ∞ of ( exp(-u*u)/(u+x) ) du
func Goodst(x float64) float64 {
	return toms.GOODST(x)
}

// Lobach calculates the the Lobachewsky function L(x)
//   ∫ 0 to x of -ln | cos t | dt
func Lobach(x float64) float64 {
	return toms.LOBACH(x)
}

// Strom calculates Stromgren's integral
//   ∫ 0 to x { t^7 exp(2t)/[exp(t)-1]^3 } dt
func Strom(x float64) float64 {
	return toms.STROM(x)
}

// Synch calculates the synchrotron radiation function
// For order 1:
//   x ∫ 0 to infinity { K(5/3)(t) } dt
//   where K(5/3) is a modified Bessel function of order 5/3.
// For order 2:
//   x * K(2/3)(x)
//   where K(2/3) is a modified Bessel function of order 2/3.
func Synch(order int, x float64) float64 {
	switch order {
	case 1:
		return toms.SYNCH1(x)
	case 2:
		return toms.SYNCH2(x)
	default:
		panic("Order must be 1 or 2")
	}
}

// Transport calculates the transport integral of order n
//  ∫ 0 to x {t^n exp(t)/[exp(t)-1]^2 } dt
func Transport(order int, x float64) float64 {
	switch order {
	case 2:
		return toms.TRAN02(x)
	case 3:
		return toms.TRAN03(x)
	case 4:
		return toms.TRAN04(x)
	case 5:
		return toms.TRAN05(x)
	case 6:
		return toms.TRAN06(x)
	case 7:
		return toms.TRAN07(x)
	case 8:
		return toms.TRAN08(x)
	case 9:
		return toms.TRAN09(x)
	default:
		panic("order must be between 2 and 9")
	}
}

// Struve calculates the Struve function of order 0 and 1, H0(x) and H1(x) respectively.
// 		H0(x) = (2/pi) ∫ {0 to pi/2} sin(x cos(t)) dt
// and
//     H1(x)  = (2/pi) ∫ {0 to pi/2} sin( x cos(t))*sin^2 t dt
func Struve(order int, x float64) float64 {
	switch order {
	case 0:
		return toms.STRVH0(x)
	case 1:
		return toms.STRVH1(x)
	default:
		panic("order must be 1 or 2")
	}
}

// StruveModified calculates the modified Struve function of order 0 and 1, L0(x) and L1(x) respectively.
// defined as the solution of the second-order equation
//     x*D(Df) + Df - x*f  =  2x/pi
// and
//     x^2*D(Df) + x*Df - (x^2+1)f = 2*x^2/pi
func StruveModified(order int, x float64) float64 {
	switch order {
	case 0:
		return toms.STRVL0(x)
	case 1:
		return toms.STRVL1(x)
	default:
		panic("order must be 1 or 2")
	}
}

// BesselMinusStruveModified calculates
// 		Ii(x) - Li(x) for i =1,2
// where Ii(x) is the modified Bessel function of the first kind of
// order i, and Li(x) is the modified Struve function of order
func BesselMinusStruveModified(order int, x float64) float64 {
	switch order {
	case 0:
		return toms.I0ML0(x)
	case 1:
		return toms.I1ML1(x)
	default:
		panic("order must be 1 or 2")
	}
}

// AtnInt calculates the value of the inverse-tangent integral defined by
//   ∫ 0 to x ( (arctan t)/t ) dt
func AtnInt(x float64) float64 {
	return toms.ATNINT(x)
}

// Exp3 calculates the value of ∫ 0 to x (exp(-t*t*t)) dt
func Exp3(x float64) float64 {
	return toms.EXP3(x)
}

// I0Int calculates the value of ∫ 0 to x I0(t) dt
func I0Int(x float64) float64 {
	return toms.I0INT(x)
}

// J0Int calculates the value of ∫ 0 to x J0(t) dt
func J0Int(x float64) float64 {
	return toms.J0INT(x)
}

// Y0Int calculates the value of ∫ 0 to x Y0(t) dt
func Y0Int(x float64) float64 {
	return toms.Y0INT(x)
}

// K0Int calculates the value of ∫ 0 to x K0(t) dt
func K0Int(x float64) float64 {
	return toms.K0INT(x)
}

// AiInt calculates the integral of the Airy function Ai,
//    ∫  0 to x Ai(t) dt
func AiInt(x float64) float64 {
	return toms.AIRINT(x)
}

// BiInt calculates the integral of the Airy function Bi,
//    ∫  0 to x Bi(t) dt
func BiInt(x float64) float64 {
	return toms.BIRINT(x)
}
