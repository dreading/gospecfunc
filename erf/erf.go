// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package erf

import (
	"github.com/dreading/gospecfunc/erf/internal/toms"
	"math"
	"math/cmplx"
)

// Erf computes approximate values for the error function
func Erf(z complex128) complex128 {
	return 1 - toms.Faddeyeva(1i*z)*cmplx.Exp(-z*z)
}

// Erfc computes approximate values for the complementary error function erfc(z) = 1 - erf(z)
func Erfc(z complex128) complex128 {
	return toms.Faddeyeva(1i*z) * cmplx.Exp(-z*z)
}

// Erfcx computes approximate values for the scaled complementary error function erfcx(z) = exp(z^2) * erfc(z)
func Erfcx(z complex128) complex128 {
	return toms.Faddeyeva(1i * z)
}

// Erfi computes approximate values for the imaginary error function erfi(z) = -i*erf(iz).
func Erfi(z complex128) complex128 {
	return -1i * Erf(1i*z)
}

// Dawson computes approximate values for the Dawson function (integral). Dawson function is the one-sided Fourier–Laplace sine transform of the Gaussian function.
func Dawson(z complex128) complex128 {
	return complex(0, 0.5*math.Sqrt(math.Pi)) * (cmplx.Exp(-z*z) - toms.Faddeyeva(z))
}

// Fresnel computes approximate values for the cos and sin Fresnel integral int_0^x cos(t^2) dt and integral int_0^x sin(t^2) dt
func Fresnel(z complex128) (complex128, complex128) {
	var halfSqrtPi = 0.5 * math.Sqrt(math.Pi)
	var erfPlus = complex(0.5, 0.5) * Erf(z*complex(halfSqrtPi, -halfSqrtPi))
	var erfNeg = complex(0.5, -0.5) * Erf(z*complex(halfSqrtPi, halfSqrtPi))
	return complex(0.5, 0) * (erfNeg + erfPlus), complex(0, 0.5) * (erfNeg - erfPlus)
}

// Voigt computes approximate values for the real and imaginary Voigt functions: https://dlmf.nist.gov/7.19
// Here we use Faddeyeva to provide analytical continuation for all t in R
func Voigt(x float64, t float64) (float64, float64) {
	//Limiting values as t -> 0
	if t <= math.Nextafter(0, 1) {
		var onexx = 1 + x*x
		return 1 / onexx, x / onexx
	}
	var one2SqrtT = 1 / (2 * math.Sqrt(t))
	var z = complex(one2SqrtT, -x*one2SqrtT)
	var c = complex(math.Sqrt(math.Pi/(4*t)), 0) * toms.Faddeyeva(1i*z)
	return real(c), imag(c)
}

// Faddeyeva computes the plasma dispersion Faddeyeva function, w(z) = exp(-z^2) * erfc(-i*z)
// where z=x+iy and erfc(z) is the complex complementary error function of z
func Faddeyeva(z complex128) complex128 {
	return toms.Faddeyeva(z)
}
