// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package airy_test

import (
	. "github.com/dreading/gospecfunc/airy"
	"testing"
)

// Global exported variables are used to store the
// return values of functions measured in the benchmarks.
// Storing the results in these variables prevents the compiler
// from completely optimizing the benchmarked functions away.
var (
	GlobalF float64
	GlobalC complex128
)

func BenchmarkAiryAi(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Ai(0.8 + 5i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiNearZero(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Ai(0.001 + 0.00005i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiLarge(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Ai(20.2 + 33.22i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiLargeNegative(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Ai(-20.2 - 33.22i)
	}
	GlobalC = ζ
}

func BenchmarkAiryBi(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Bi(0.8 + 5i)
	}
	GlobalC = ζ
}

func BenchmarkAiryBiNearZero(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Bi(0.001 + 0.00005i)
	}
	GlobalC = ζ
}

func BenchmarkAiryBiLarge(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Bi(20.2 + 33.22i)
	}
	GlobalC = ζ
}

func BenchmarkAiryBiLargeNegative(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = Bi(-20.2 - 33.22i)
	}
	GlobalC = ζ
}

func BenchmarkAiryIntegral(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = AiInt(0.12)
	}
	GlobalF = ζ
}

func BenchmarkAiryGi(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Gi(3.2)
	}
	GlobalF = ζ
}

func BenchmarkAiryHi(b *testing.B) {
	var ζ float64
	for n := 0; n < b.N; n++ {
		ζ = Hi(3.2)
	}
	GlobalF = ζ
}
