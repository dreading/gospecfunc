// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc

import (
	"github.com/dreading/gospecfunc/integrals/internal/toms"
)

// Abramowitz0 computes the  Abramowitz function of order 0
//   ∫ 0 to infinity exp( -t*t - x/t ) dt
func Abramowitz0(x float64) float64 {
	return toms.ABRAM0(x)
}

// Abramowitz1 computes the  Abramowitz function of order 0
//   ∫ 0 to infinity t * exp( -t*t - x/t ) dt
func Abramowitz1(x float64) float64 {
	return toms.ABRAM1(x)
}

// Abramowitz2 calculates the Abramowitz function of order 2,
//    ∫  0 to infinity (t^2) * exp( -t*t - x/t ) dt
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
