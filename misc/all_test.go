// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc_test

import (
	. "github.com/dreading/gospecfunc/misc"
	"testing"
)

func TestAbramowitz0(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{1.0e0, 512.0e0, 0.87377726306985360531e0},
		{1.0e0, 128.0e0, 0.84721859650456925922e0},
		{1.0e0, 32.0e0, 0.77288934483988301615e0},
		{1.0e0, 8.0e0, 0.59684345853450151603e0},
		{1.0e0, 2.0e0, 0.29871735283675888392e0},
		{1.0e0, 1.0e0, 0.15004596450516388138e0},
		{5.0e0, 4.0e0, 0.11114662419157955096e0},
		{3.0e0, 2.0e0, 0.83909567153151897766e-1},
		{15.0e0, 8.0e0, 0.56552321717943417515e-1},
		{2.0e0, 1.0e0, 0.49876496603033790206e-1},
		{17.0e0, 8.0e0, 0.44100889219762791328e-1},
		{3.0e0, 1.0e0, 0.19738535180254062496e-1},
		{4.0e0, 1.0e0, 0.86193088287161479900e-2},
		{5.0e0, 1.0e0, 0.40224788162540127227e-2},
		{6.0e0, 1.0e0, 0.19718658458164884826e-2},
		{7.0e0, 1.0e0, 0.10045868340133538505e-2},
		{10.0e0, 1.0e0, 0.15726917263304498649e-3},
		{15.0e0, 1.0e0, 0.10352666912350263437e-4},
		{20.0e0, 1.0e0, 0.91229759190956745069e-6},
		{40.0e0, 1.0e0, 0.25628287737952698742e-9},
	}

	for _, tc := range testCases {
		ζ := Abramowitz0(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Abramowitz0(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
		}

	}
}

// The floating point comparison tests are copied from from math/all_test.go.
func tolerance(a, b, e float64) bool {
	// Multiplying by e here can underflow denormal values to zero.
	// Check a==b so that at least if a and b are small and identical
	// we say they match.
	if a == b {
		return true
	}
	d := a - b
	if d < 0 {
		d = -d
	}
	// note: b is correct (expected) value, a is actual value.
	// make error tolerance a fraction of b, not a.
	if b != 0 {
		e = e * b
		if e < 0 {
			e = -e
		}
	}
	return d < e
}

func close(a, b float64) bool      { return tolerance(a, b, 1e-14) }
func veryclose(a, b float64) bool  { return tolerance(a, b, 5e-16) }
func soclose(a, b, e float64) bool { return tolerance(a, b, e) }
