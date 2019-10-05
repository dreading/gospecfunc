// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package machine

import (
	"math"
)

// I1MACH PROVIDES INTEGER MACHINE DEPENDENT CONSTANTS
var I1MACH = []int{-0, 5, 6, 0, 0, 32, 4, 2, 31, 2147483647, 2, 24, -125, 127, 53, -1021, 1023}

// D1MACH PROVIDES DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS
var D1MACH = []float64{math.NaN(), 2.23e-308, 1.79e-308, 1.11e-16, 2.22e-16, 0.30103000998497009}
