// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package airy

import (
	"github.com/dreading/gospecfunc/airy/internal/amos"
)

// Ai computes the Airy Ai function
func Ai(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 0, 1)
	return complex(air, aii)
}

// AiE computes the exponentially scaled Airy Ai function 
// AiE =e^((2/3)*z^(3/2)) * AI(z)
func AiE(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 0, 2)
	return complex(air, aii)
}

// AiD computes the first derivate of the Airy Ai function dAi(z)/dz
func AiD(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 1, 1)
	return complex(air, aii)
}

// AiDE computes the exponentially scaled first derivate of the Airy Ai function 
// AiDE = e^((2/3)*z^(3/2)) * dAi(z)/dz
func AiDE(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 1, 2)
	return complex(air, aii)
}

// Bi computes the Airy Bi function
func Bi(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 0, 1)
	return complex(bir, bii)
}

// Bi computes the expotentially scaled Airy Bi function
// BiE =e^(-|Real((2/3)*z^(3/2))|) * BI(z)
func BiE(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 0, 2)
	return complex(bir, bii)
}

// BiD computes the first derivate of the Airy Bi function dBi(z)/dz
func BiD(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 1, 1)
	return complex(bir, bii)
}

// BiD computes the expotentially scaled first derivate of the Airy Bi function 
// BiE =e^(-|Real((2/3)*z^(3/2))|) * dBi(z)/dz
func BiDE(z complex128) complex128 {
	bir, bii, _ := amos.ZBIRY(real(z), imag(z), 1, 2)
	return complex(bir, bii)
}
