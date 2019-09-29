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
	GlobalC complex128
)

func BenchmarkAiryAi(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = AiryAi(0.8 + 5i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiNearZero(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = AiryAi(0.001 + 0.00005i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiLarge(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = AiryAi(20.2 + 33.22i)
	}
	GlobalC = ζ
}

func BenchmarkAiryAiLargeNegative(b *testing.B) {
	var ζ complex128
	for n := 0; n < b.N; n++ {
		ζ = AiryAi(-20.2 - 33.22i)
	}
	GlobalC = ζ
}
