// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// The code below are based on
// Algorithm 916: Computing the Faddeyeva and Voigt Functions
// by Zaghloul, Mofreh R. and Ali, Ahmed N.
// ACM Trans. Math. Softw. December 2011 Vol 38 No 2 Jan 2012
// http://doi.acm.org/10.1145/2049673.2049679.
//
// The go code is a port of the original matlab code.

package toms

import (
	"github.com/dreading/gospecfunc/erf/internal/libcerf"
	"math"
	"math/cmplx"
)

// Faddeyeva computes the plasma dispersion Faddeyeva function, w(z) = exp(-z^2) * erfc(-i*z)
// where z=x+iy and erfc(z) is the complex complementary error function of z
func Faddeyeva(z complex128) complex128 {

	//For purely imaginary input values
	if IsReal(1i * z) {
		return complex(libcerf.Erfcx(imag(z)), 0)
	}

	var Rmin float64 = 1 / math.MaxFloat64 //Nextafter(0, 1)
	var sqrtLogRmin float64 = math.Sqrt(-math.Log(Rmin))

	//Faddeyeva cannot calculate w for y negative & exp(y^2-x^2)>=the largest
	//floating point number due to unavoidable over flow problems
	if imag(z) < 0 && (imag(z)*imag(z)-real(z)*real(z)) >= sqrtLogRmin*sqrtLogRmin {
		return cmplx.Inf()
	}

	var eps float64 = 1e-16
	var tiny float64 = 0.6447 * eps
	var a float64 = math.Sqrt(-1 * math.Pi * math.Pi / math.Log(tiny/2))

	var halfA float64 = 0.5 * a
	var aSqr float64 = a * a
	var twoA float64 = 2 * a
	var twoASqr float64 = 2 * aSqr
	var fourASqr float64 = 2 * twoASqr
	var aPi float64 = a / math.Pi
	var twoAPi float64 = 2 * aPi
	var oneSqrtPi float64 = 1 / math.Sqrt(math.Pi)
	var erfcxY float64 = libcerf.Erfcx(math.Abs(imag(z)))

	var x float64 = real(z)
	var y float64 = imag(z)
	var xsign float64 = Sign(x)
	var ysign float64 = Sign(y)
	x = xsign * x
	y = math.Max(Rmin, ysign*y)

	var xSqr float64 = x * x
	var ySqr float64 = y * y
	var twoYx float64 = 2 * y * x
	var twoAX float64 = twoA * x
	var expXSqr float64 = math.Exp(-xSqr)
	var cos2yx float64 = math.Cos(twoYx)
	var sin2yx float64 = math.Sin(twoYx)

	var w complex128
	if x == 0 {
		//% For x=0, use the asymptotic expressions for x--->0 from Eq. (6) in the manuscript
		w = complex(erfcxY, 0) * (complex(1, xsign*(x/y)))
	} else {

		var vOld float64 = expXSqr * (erfcxY*cos2yx + (twoAPi/y)*(math.Sin(twoYx/2)*math.Sin(twoYx/2)))
		var lOld float64 = -erfcxY + aPi/y
		var Sigma3 float64 = Rmin
		var Sigma5 float64 = Rmin
		var Sigma45 float64
		var delta3 float64 = 1
		var delta5 float64 = 1
		var n float64
		var n3 float64 = math.Ceil(x / a)
		var n33 = n3 - 1

		if (sqrtLogRmin - x) > 0 {
			var Sigma1 float64 = Rmin
			var Sigma2 float64 = Rmin
			var Sigma4 float64 = Rmin
			var exp1 float64 = math.Exp(-twoAX)
			var exp2 float64 = math.Exp((fourASqr*n3 - 2*twoAX) - twoASqr)
			var exp3 float64 = math.Exp(-((twoASqr*n3 - twoAX) - twoASqr))
			var del2tmp float64 = 1.0
			var del3tmp float64 = math.Exp(-((aSqr*n3*n3 - twoAX*n3 - twoASqr*n3) + xSqr + twoAX + aSqr))
			var del33tmp float64 = math.Exp(aSqr - (twoASqr*n3 - twoAX))
			var del3 float64
			var del5 float64
			var den1 float64
			var expDel1 float64
			var exp3Den float64
			var exp33Den float64

			for delta3 >= tiny || (delta5 >= tiny && n <= 50) {
				n = n + 1
				den1 = aSqr*n*n + ySqr
				expDel1 = math.Exp(-(aSqr * n * n)) / den1
				del2tmp = del2tmp * exp1
				if n33 >= 1 {
					del3tmp = del3tmp * exp3
					exp3Den = del3tmp * expDel1 * (den1 / (aSqr*n3*n3 + ySqr))
					del33tmp = del33tmp * exp2
					exp33Den = exp3Den * del33tmp * ((aSqr*n3*n3 + ySqr) / (aSqr*n33*n33 + ySqr))
					del5 = (n33*exp33Den + n3*exp3Den)
					del3 = (exp33Den + exp3Den)
				} else {
					del3tmp = del3tmp * exp3
					del3 = del3tmp * expDel1 * (den1 / (aSqr*n3*n3 + ySqr))
					del5 = (n3 * del3)
				}
				delta3 = del3 / Sigma3
				delta5 = del5 / Sigma5

				Sigma1 = Sigma1 + expDel1
				Sigma2 = Sigma2 + del2tmp*expXSqr*expDel1
				Sigma3 = Sigma3 + del3
				Sigma4 = Sigma4 + n*del2tmp*expXSqr*expDel1
				Sigma5 = Sigma5 + del5

				if x >= 5e-4 {
					Sigma45 = -Sigma4 + Sigma5
				} else {
					var twoAXn2 float64 = twoAX * n * twoAX * n
					Sigma45 = Sigma45 + 2*n*n*twoAX*expXSqr*expDel1*(1+1.666666666666667e-001*(twoAXn2)+8.333333333333333e-003*(twoAXn2)*(twoAXn2))
				}
				n3 = n3 + 1
				n33 = n33 - 1
			}

			if y <= 5e0 && twoYx > Rmin {
				var wReal = vOld + y*twoAPi*(-cos2yx*expXSqr*Sigma1+0.5*(Sigma2+Sigma3))
				var wImag = xsign * (sin2yx*expXSqr*(lOld+twoAPi*y*Sigma1) + twoAPi*halfA*Sigma45)
				w = complex(wReal, wImag)
			} else if y <= 5e0 && twoYx <= Rmin {
				var wReal = vOld + y*twoAPi*(-cos2yx*expXSqr*Sigma1+.5*(Sigma2+Sigma3))
				var wImag = xsign * (2*y*expXSqr*(x*lOld+x*twoAPi*y*Sigma1) + twoAPi*halfA*Sigma45)
				w = complex(wReal, wImag)
			} else {
				var wReal = vOld + y*twoAPi*(-cos2yx*expXSqr*Sigma1+0.5*(Sigma2+Sigma3))
				var wImag = xsign * (sin2yx*expXSqr*math.Min(0, math.Abs(lOld+(twoAPi*y*Sigma1))) + twoAPi*halfA*Sigma45)
				w = complex(wReal, wImag)
			}

		} else if x >= sqrtLogRmin && x < 1e15 {
			var exp2 = math.Exp((fourASqr*n3 - 2*twoAX) - twoASqr)
			var del33tmp = math.Exp(aSqr + (twoAX - twoASqr*n3))
			var del3 float64
			var del5 float64
			for delta3 >= tiny || delta5 >= tiny && n <= 50 {
				n = n + 1
				if n33 >= 1 {
					var exp3Den = math.Exp(-(a*n3-x)*(a*n3-x)) / (aSqr*n3*n3 + ySqr)
					del33tmp = del33tmp * exp2
					var exp33Den = exp3Den * del33tmp * ((aSqr*n3*n3 + ySqr) / (aSqr*n33*n33 + ySqr))
					del5 = (n33*exp33Den + n3*exp3Den)
					del3 = (exp33Den + exp3Den)
				} else {
					del3 = math.Exp(-(a*n3-x)*(a*n3-x)) / (aSqr*n3*n3 + ySqr)
					del5 = n3 * del3
				}
				delta3 = del3 / Sigma3
				delta5 = del5 / Sigma5
				Sigma3 = Sigma3 + del3
				Sigma5 = Sigma5 + del5
				n3 = n3 + 1
				n33 = n33 - 1
			}

			var wReal = vOld + y*aPi*Sigma3
			var wImag = xsign * (sin2yx*expXSqr*lOld + twoAPi*halfA*Sigma5)
			w = complex(wReal, wImag)

		} else {
			var wReal = oneSqrtPi * (y / (xSqr + ySqr))
			var wImag = oneSqrtPi * (xsign * x / (xSqr + ySqr))
			w = complex(wReal, wImag)
		}
	}

	if ysign < 0 {
		var twoExpXSqrYsqr = 2 * math.Exp(-xSqr+ySqr)
		var wReal = twoExpXSqrYsqr*cos2yx - real(w)
		var wImag = xsign*twoExpXSqrYsqr*sin2yx + imag(w)
		w = complex(wReal, wImag)
	}

	return w
}

// IsReal returns true if the complex number is real
func IsReal(cmplx complex128) bool {
	return imag(cmplx) == 0
}

// Sign returns the sign of float
func Sign(x float64) float64 {
	switch {
	case math.IsNaN(x):
		return math.NaN()
	case math.IsInf(x, 1):
		return 1
	case math.IsInf(x, -1):
		return -1
	case x == 0:
		return 0
	case x < 0:
		return -1
	}
	return 1
}
