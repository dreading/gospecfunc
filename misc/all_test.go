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

func TestAbramowitz1(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{1.0e0, 512.0e0, 0.49828219848799921792e0},
		{1.0e0, 128.0e0, 0.49324391773047288556e0},
		{1.0e0, 32.0e0, 0.47431612784691234649e0},
		{1.0e0, 8.0e0, 0.41095983258760410149e0},
		{1.0e0, 2.0e0, 0.25317617388227035867e0},
		{1.0e0, 1.0e0, 0.14656338138597777543e0},
		{5.0e0, 4.0e0, 0.11421547056018366587e0},
		{3.0e0, 2.0e0, 0.90026307383483764795e-1},
		{15.0e0, 8.0e0, 0.64088214170742303375e-1},
		{2.0e0, 1.0e0, 0.57446614314166191085e-1},
		{17.0e0, 8.0e0, 0.51581624564800730959e-1},
		{3.0e0, 1.0e0, 0.25263719555776416016e-1},
		{4.0e0, 1.0e0, 0.11930803330196594536e-1},
		{5.0e0, 1.0e0, 0.59270542280915272465e-2},
		{6.0e0, 1.0e0, 0.30609215358017829567e-2},
		{7.0e0, 1.0e0, 0.16307382136979552833e-2},
		{10.0e0, 1.0e0, 0.28371851916959455295e-3},
		{15.0e0, 1.0e0, 0.21122150121323238154e-4},
		{20.0e0, 1.0e0, 0.20344578892601627337e-5},
		{40.0e0, 1.0e0, 0.71116517236209642290e-9},
	}

	for _, tc := range testCases {
		ζ := Abramowitz1(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Abramowitz1(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
		}

	}
}

func TestAbramowitz2(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{1.0e0, 512.0e0, 0.44213858162107913430e0},
		{1.0e0, 128.0e0, 0.43923379545684026308e0},
		{1.0e0, 32.0e0, 0.42789857297092602234e0},
		{1.0e0, 8.0e0, 0.38652825661854504406e0},
		{1.0e0, 2.0e0, 0.26538204413231368110e0},
		{1.0e0, 1.0e0, 0.16848734838334595000e0},
		{5.0e0, 4.0e0, 0.13609200032513227112e0},
		{3.0e0, 2.0e0, 0.11070330027727917352e0},
		{15.0e0, 8.0e0, 0.82126019995530382267e-1},
		{2.0e0, 1.0e0, 0.74538781999594581763e-1},
		{17.0e0, 8.0e0, 0.67732034377612811390e-1},
		{3.0e0, 1.0e0, 0.35641808698811851022e-1},
		{4.0e0, 1.0e0, 0.17956589956618269083e-1},
		{5.0e0, 1.0e0, 0.94058737143575370625e-2},
		{6.0e0, 1.0e0, 0.50809356204299213556e-2},
		{7.0e0, 1.0e0, 0.28149565414209719359e-2},
		{10.0e0, 1.0e0, 0.53808696422559303431e-3},
		{15.0e0, 1.0e0, 0.44821756380146327259e-4},
		{20.0e0, 1.0e0, 0.46890678427324100410e-5},
		{40.0e0, 1.0e0, 0.20161544850996420504e-8},
	}

	for _, tc := range testCases {
		ζ := Abramowitz2(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Abramowitz2(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
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
