[![Build Status](https://travis-ci.org/dreading/gospecfunc.svg?branch=master)](https://travis-ci.org/dreading/gospecfunc) 
[![Go Report Card](https://goreportcard.com/badge/github.com/dreading/gospecfunc)](https://goreportcard.com/report/github.com/dreading/gospecfunc)
[![stability-unstable](https://img.shields.io/badge/stability-unstable-yellow.svg)](https://github.com/emersion/stability-badges#unstable)

# Installation
The core packages of the gospecfunc suite are written in pure Go. Installation is done using go get.
```
go get -u github.com/dreading/gospecfunc
```

# Special Functions

## Bessel

Functions  | Domain | Description |
---------- | ------ | ----------- |
Ai     |  ℂ  | Airy Ai  function |
AiE    |  ℂ   | Exponentially scaled Airy Ai function|
AiD    |  ℂ   | First derivative of the Airy Ai function|
AiDE    |  ℂ   | Exponentially scaled first derivative of the Airy Ai function|
Bi     |  ℂ  | Biry Bi function |
BiE    |  ℂ   | Exponentially scaled Biry Bi function|
BiD    |  ℂ   | First derivative of the Biry Bi function|
BiDE    |  ℂ   | Exponentially scaled first derivative of the Biry Bi function|
Gi     |  ℝ  | Modified Airy Gi  function |
Hi     |  ℝ  | Modified Airy Hi  function |
I      | ℂ  | Modified Bessel function of the first kind  |
J      | ℂ  | Bessel function of the first kind |
K      | ℂ  | Modified Bessel function of the second kind  |
Y      | ℂ  | Bessel function of the second kind  |
H1      | ℂ  | Hankel fucntion of of the first kind  |
H2      | ℂ  | Hankel fucntion of of the second kind  |

## Erf

## Integrals


# Testing 
```
 go test ./*/. 
```
# Benchmarking
```
 go test -bench=. ./*/.
```
# License
Original code is licensed under the Gonum License found in the LICENSE file. Portions of the code are subject to the additional licenses found in THIRD_PARTY_LICENSES.  
