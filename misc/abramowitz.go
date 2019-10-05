// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc

import (
	"github.com/dreading/gospecfunc/misc/internal/toms"
)

// Abramowitz0 computes the  Abramowitz function of order 0
//   ∫ 0 to infinity exp( -t*t - x/t ) dt
func Abramowitz0(x float64) float64 {
	return toms.ABRAM0(x)
}
