// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package misc_test

import (
	. "github.com/dreading/gospecfunc/integrals"
	"testing"
)

// Global exported variables are used to store the
// return values of functions measured in the benchmarks.
// Storing the results in these variables prevents the compiler
// from completely optimizing the benchmarked functions away.
var (
	GlobalF float64
)

func BenchmarkAbramowitz0(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Abramowitz(0, 0.8)
	}
	GlobalF = ζ
}

func BenchmarkAbramowitz1(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Abramowitz(1, 0.8)
	}
	GlobalF = ζ
}

func BenchmarkAbramowitz2(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Abramowitz(2, 0.8)
	}
	GlobalF = ζ
}

func BenchmarkClausen(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Clausen(0.8)
	}
	GlobalF = ζ
}
