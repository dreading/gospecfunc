// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
package misc

import (
	"github.com/dreading/gospecfunc/misc/internal/toms"
)

// Cheval evaluates a Chebyshev series, using the Clenshaw method
// with Reinsch modification
func Cheval(n int, A []float64, t float64) float64 {
	return toms.CHEVAL(n, A, t)
}
