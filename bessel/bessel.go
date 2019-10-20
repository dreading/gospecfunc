// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package bessel

import (
	"github.com/dreading/gospecfunc/bessel/internal/amos"
	"math"
)

// I computes the bessel function I for complex arguments
func I(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, CYR, CYI, NZ, IERR = amos.ZBESI(real(z), imag(z), α, 1, 1, CYR, CYI, NZ, IERR)
	return complex(CYR[1], CYI[1])
}

// J computes the bessel function J for complex arguments
func J(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, CYR, CYI, NZ, IERR = amos.ZBESJ(real(z), imag(z), α, 1, 1, CYR, CYI, NZ, IERR)
	return complex(CYR[1], CYI[1])
}

// K computes the bessel function K for complex arguments
func K(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, CYR, CYI, NZ, IERR = amos.ZBESK(real(z), imag(z), α, 1, 1, CYR, CYI, NZ, IERR)
	return complex(CYR[1], CYI[1])
}

// Y computes the bessel function Y for complex arguments
func Y(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	CWRKR := []float64{math.NaN(), 0}
	CWRKI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, CYR, CYI, NZ, _, _, IERR = amos.ZBESY(real(z), imag(z), α, 1, 1, CYR, CYI, NZ, CWRKR, CWRKI, IERR)
	return complex(CYR[1], CYI[1])
}

// H1 computes the hankel function of order 1 for complex arguments
func H1(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, _, CYR, CYI, NZ, IERR = amos.ZBESH(real(z), imag(z), α, 1, 1, 1, CYR, CYI, NZ, IERR)
	return complex(CYR[1], CYI[1])
}

// H2 computes the hankel function of order 1 for complex arguments
func H2(α float64, z complex128) complex128 {
	CYR := []float64{math.NaN(), 0}
	CYI := []float64{math.NaN(), 0}
	var NZ, IERR int
	_, _, _, _, _, _, CYR, CYI, NZ, IERR = amos.ZBESH(real(z), imag(z), α, 1, 2, 1, CYR, CYI, NZ, IERR)
	return complex(CYR[1], CYI[1])
}
