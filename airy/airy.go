// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package airy

import (
	"github.com/dreading/gospecfunc/airy/internal"
)

func AiryAi(z complex128) complex128 {
	air, aii, _, _ := amos.ZAIRY(real(z), imag(z), 0, 1)
	return complex(air, aii)
}
