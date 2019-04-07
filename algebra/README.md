# algebra
 sounds cooler than math
(and math is too general)

## assumptions
you already know basic concepts: gcd, lcm, prime, modulo, factorization

## covered

### modular arithmetic
- difficulty: easy
- content: primes, polynomial hash, extended gcd, modular inverse,
  binary exponention, chinese remainder theorem
### Fast Fourier Transform
- difficulty: medium
- content: 3 variants of fft
### Möbius function
- difficulty: medium
- content: Möbius function [wiki][mob], prime density [wiki][dense]

## not covered

### diophantine equation
- find all ![eq][dioph]
- hint: compute gcd
### Number Theoretic Transform
- fft without floating point numbers
- todo
### faster primality tests
- naive in O(sqrt n)
- todo Miller-Rabin O(log^2 n)
### linear time prime sieve
- using least prime factor
- todo
### discrete log
- find ![eq][dislog1]
- meet in the middle: ![eq][dislog2]
- todo + discrete root n stuff

[mob]: https://en.wikipedia.org/wiki/M%C3%B6bius_function "Möbius function"
[dense]: https://en.wikipedia.org/wiki/Prime-counting_function "Prime count"
[dioph]: https://latex.codecogs.com/gif.latex?%5Cinline%20%28x%2Cy%29%20%5C%20%7C%20%5C%20ax&plus;by%3Dc
[dislog1]: https://latex.codecogs.com/gif.latex?%5Cinline%20x%20%5C%20%7C%20%5C%20a%5Ex%20%5Cequiv%20b%20%5Cpmod%20m%2C%20%5C%20%5Cgcd%28a%2Cm%29%3D1
[dislog2]: https://latex.codecogs.com/gif.latex?%5Cinline%20a%5E%7Bp%20%5Csqrt%7Bm%7D%7D%20%5Cequiv%20a%5E%7Bq%7Db%20%5Cpmod%20m

