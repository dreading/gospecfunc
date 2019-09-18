package erf

import (
	"math"
	"math/cmplx"
)


func Erf(x complex128) complex128 {
	// Erfc computes the  complementary error function erfc(z) = 1 - erf(z).
	return 1 - Faddeyeva(1i*x)*cmplx.Exp(-x*x)
}

func Erfc(x complex128) complex128 {
	// Erfc computes the  complementary error function erfc(z) = 1 - erf(z).
	return Faddeyeva(1i*x)*cmplx.Exp(-x*x)
}

func Erfcx(x complex128) complex128 {
	// Erfcx computes the scaled complementary error function erfcx(z) = exp(z^2) * erfc(z).
	return Faddeyeva(1i*x)
}

func Erfi(x complex128) complex128 {
	// Erfc computes the  complementary error function erfc(z) = 1 - erf(z).
	return -1i*Erf(1i*x)
}

func Dawson(x complex128) complex128 {
	// Erfc computes the  complementary error function erfc(z) = 1 - erf(z).
	return  complex(0,0.5*math.Sqrt(math.Pi)) * (cmplx.Exp(-x*x) - Faddeyeva(x))
}
