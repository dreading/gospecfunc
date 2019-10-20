// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package bessel

import (
	"github.com/dreading/gospecfunc/bessel/internal/amos"
	"github.com/dreading/gospecfunc/bessel/internal/toms"
)

// Ai computes the Airy Ai function
func Ai(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 0, 1)
	return complex(air, aii)
}

// Aix computes the exponentially scaled Airy Ai function
// Aix = e^((2/3)*z^(3/2)) * AI(z)
func Aix(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 0, 2)
	return complex(air, aii)
}

// Aid computes the first derivative of the Airy Ai function dAi(z)/dz
func Aid(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 1, 1)
	return complex(air, aii)
}

// Aidx computes the exponentially scaled first derivative of the Airy Ai function
// Aidx = e^((2/3)*z^(3/2)) * dAi(z)/dz
func Aidx(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 1, 2)
	return complex(air, aii)
}

// Bi computes the Airy Bi function
func Bi(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 0, 1)
	return complex(bir, bii)
}

// Bix computes the expotentially scaled Airy Bi function
// Bix = e^(-|Real((2/3)*z^(3/2))|) * BI(z)
func Bix(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 0, 2)
	return complex(bir, bii)
}

// Bid computes the first derivative of the Airy Bi function dBi(z)/dz
func Bid(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 1, 1)
	return complex(bir, bii)
}

// Bidx computes the expotentially scaled first derivative of the Airy Bi function
// Bidx = e^(-|Real((2/3)*z^(3/2))|) * dBi(z)/dz
func Bidx(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 1, 2)
	return complex(bir, bii)
}

// Gi calculates the modified Airy function Gi
//    ∫  0 to infinity sin(x*t+t^3/3) dt  / pi
func Gi(x float64) float64 {
	return toms.AIRYGI(x)
}

// Hi calculates the  modified Airy function Hi
//    ∫  0 to infinity exp(x*t-t^3/3) dt  / pi
func Hi(x float64) float64 {
	return toms.AIRYHI(x)
}
