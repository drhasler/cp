/* git@drhasler - 04/19 */
#include <bits/stdc++.h>
using namespace std;

#define ll long long
#define pb push_back
const int N = 1e6;

/* Möbius function mu(n)
 * O(nlog2n) need to check it
 *
 * applications:
 * 1. inverse möbius formula
 * g(n) = sum_{d|n} f(d) implies
 * f(n) = sum_{d|n} mu(d) g(n/d)
 * 2. problem 1139D
 * expectation and inclusion-exclusion
 * over GCD distinct prime factors */

int mu[N];
vector<int> primes;
bool prime[N];
void mobius() {
    fill(prime+2,prime+N,1);
    for (int i=2;i<N;i++) {
        if (prime[i]) primes.pb(i), mu[i]=-1;
        for (auto &p:primes) {
            if (i*p>=N) break;
            prime[i*p]=0; // not prime
            if (i%p==0) { mu[i*p]=0; break; } // mult > 1
            else mu[i*p]=mu[i]*mu[p]; // distinct or 0
        }
    }
}

int main() {
    // pi(N) < N/10 < 2 pi(N) for N in 1e5..1e9
    primes.reserve(N/10);
    mobius();
}
