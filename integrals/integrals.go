// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc

import (
	"github.com/dreading/gospecfunc/integrals/internal/toms"
)

// Abramowitz0 computes the  Abramowitz function of order 0
//   ∫ 0 to ∞ of exp( -t*t - x/t ) dt
func Abramowitz0(x float64) float64 {
	return toms.ABRAM0(x)
}

// Abramowitz1 computes the  Abramowitz function of order 0
//   ∫ 0 to ∞ of t * exp( -t*t - x/t ) dt
func Abramowitz1(x float64) float64 {
	return toms.ABRAM1(x)
}

// Abramowitz2 calculates the Abramowitz function of order 2,
//    ∫  0 to ∞ of (t^2) * exp( -t*t - x/t ) dt
func Abramowitz2(x float64) float64 {
	return toms.ABRAM2(x)
}

// Clausen calculates the Clausen's integral,
//   ∫ 0 to x of (-ln(2*sin(t/2))) dt
func Clausen(x float64) float64 {
	return toms.CLAUSN(x)
}

// Debye1 calculates the Debye function of order 1, defined as
//   ∫ 0 to x of t/(exp(t)-1) dt] / x
func Debye1(x float64) float64 {
	return toms.DEBYE1(x)
}

// Debye2 calculates the Debye function of order 2, defined as
//   2 * ∫ {0 to x} t^2/(exp(t)-1) dt / x^2
func Debye2(x float64) float64 {
	return toms.DEBYE2(x)
}

// Debye3 calculates the Debye function of order 3, defined as
//    3 * ∫ {0 to x} t^3/(exp(t)-1) dt / x^3
func Debye3(x float64) float64 {
	return toms.DEBYE3(x)
}

// Debye4 calculates the Debye function of order 4, defined as
//    4 * ∫ {0 to x} t^4/(exp(t)-1) dt / x^4
func Debye4(x float64) float64 {
	return toms.DEBYE4(x)
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

// Synch1 calculates the synchrotron radiation function
//   x ∫ 0 to infinity { K(5/3)(t) } dt
// where K(5/3) is a modified Bessel function of order 5/3.
func Synch1(x float64) float64 {
	return toms.SYNCH1(x)
}

// Synch2 calculates the synchrotron radiation function
//   x * K(2/3)(x)
// where K(2/3) is a modified Bessel function of order 2/3.
func Synch2(x float64) float64 {
	return toms.SYNCH2(x)
}

// Tran02 calculates the transport integral of order 2
//  ∫ 0 to x {t^2 exp(t)/[exp(t)-1]^2 } dt
func Tran02(x float64) float64 {
	return toms.TRAN02(x)
}

// Tran03 calculates the transport integral of order 3
//  ∫ 0 to x {t^3 exp(t)/[exp(t)-1]^2 } dt
func Tran03(x float64) float64 {
	return toms.TRAN03(x)
}

// Tran04 calculates the transport integral of order 4
//  ∫ 0 to x {t^4 exp(t)/[exp(t)-1]^2 } dt
func Tran04(x float64) float64 {
	return toms.TRAN04(x)
}

// Tran05 calculates the transport integral of order 5
//  ∫ 0 to x {t^5 exp(t)/[exp(t)-1]^2 } dt
func Tran05(x float64) float64 {
	return toms.TRAN05(x)
}

// Tran06 calculates the transport integral of order 6
//  ∫ 0 to x {t^6 exp(t)/[exp(t)-1]^2 } dt
func Tran06(x float64) float64 {
	return toms.TRAN06(x)
}

// Tran07 calculates the transport integral of order 7
//  ∫ 0 to x {t^7 exp(t)/[exp(t)-1]^2 } dt
func Tran07(x float64) float64 {
	return toms.TRAN07(x)
}

// Tran08 calculates the transport integral of order 8
//  ∫ 0 to x {t^8 exp(t)/[exp(t)-1]^2 } dt
func Tran08(x float64) float64 {
	return toms.TRAN08(x)
}

// Tran09 calculates the transport integral of order 9
//  ∫ 0 to x {t^9 exp(t)/[exp(t)-1]^2 } dt
func Tran09(x float64) float64 {
	return toms.TRAN09(x)
}
