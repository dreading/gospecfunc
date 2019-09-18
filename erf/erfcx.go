// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package erf

import (
	"math"
	"math/cmplx"
)

// Erfc computes approximate values for the complementary error function erfc(z) = 1 - erf(z)
func Erf(z complex128) complex128 {
	return 1 - Faddeyeva(1i*z)*cmplx.Exp(-z*z)
}

// Erfc computes approximate values for the complementary error function erfc(z) = 1 - erf(z)
func Erfc(z complex128) complex128 {
	return Faddeyeva(1i*z)*cmplx.Exp(-z*z)
}

// Erfcx computes approximate values for the scaled complementary error function erfcx(z) = exp(z^2) * erfc(z)
func Erfcx(z complex128) complex128 {
	return Faddeyeva(1i*z)
}

// Erfi computes approximate values for the imaginary error function erfi(z) = -i*erf(iz).
func Erfi(z complex128) complex128 {
	return -1i*Erf(1i*z)
}

// Dawson computes approximate values for the Dawson function (integral). Dawson function is the one-sided Fourier–Laplace sine transform of the Gaussian function.
func Dawson(z complex128) complex128 {
	return  complex(0,0.5*math.Sqrt(math.Pi)) * (cmplx.Exp(-z*z) - Faddeyeva(z))
}

// FresnelC computes approximate values for the Fresnel integral int_0^x cos(t^2) dt
func FresnelC(z complex128) complex128 {
	var halfSqrtPi = 0.5 * math.Sqrt(math.Pi)
	var erfPlus = complex(0.5,0.5) * Erf(z*complex(halfSqrtPi,-halfSqrtPi))
	var erfNeg =  complex(0.5,-0.5) * Erf(z*complex(halfSqrtPi, halfSqrtPi))
	return complex(0.5,0)*(erfNeg + erfPlus)
}

// FresnelS computes approximate values for the Fresnel integral int_0^x sin(t^2) dt
func FresnelS(z complex128) complex128 {
	var halfSqrtPi = 0.5 * math.Sqrt(math.Pi)
	var erfPlus = complex(0.5,0.5) * Erf(z*complex(halfSqrtPi,-halfSqrtPi))
	var erfNeg =  complex(0.5,-0.5) * Erf(z*complex(halfSqrtPi, halfSqrtPi))
	return complex(0,0.5)*(erfNeg - erfPlus)
}

