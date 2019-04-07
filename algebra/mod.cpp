/* git@drhasler - 04/19 */
#include <bits/stdc++.h>
using namespace std;

#define ll long long
#define pii pair<int,int>
#define xx first
#define yy second
const int MOD = 1e9+7;

/* Modular arithmetic
 * Summary:
 * - prime numbers
 * - modular inverse
 * - CRT */

/* Prime numbers
 * they appear a lot in modular arithmetic
 * because they provide low collision rates
 * in hashtables and yield interesting
 * algebraic properties esp in group theory
 * here are some:
 * 31, 31013, 998244353, 1e9+7, 1e9+9, (1<<31)-1
 * more at: http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
 * they are commonly used in hashing */

/* Polynomial hash
 * simplest hashing function i can think of
 * warning: it is easy to generate collisions */
int p_hash(string s) {
    ll res = 0;
    for (int a:s) res = (res*31+a) % MOD;
    return res;
}

/* Modular inverse
 * solutions:
 * 1. using Fermat's little theorem,
 * we can find the inverse of any integer
 * 0 < a < MOD: a^{-1} = a^{p-2} [MOD]
 * 2. with extended GCD we get the coefficients
 * satisfying ax + by = gcd(a,b) (1 when coprimes) */

/* modpow, using binary exponentiation */
int modpow(ll x, int p, int mod) {
    ll ans = 1;
    while (p) {
        if (p&1) ans = ans*x % mod;
        x = x*x % mod;
        p/=2;
    }
    return ans;
}

/* xgcd(a,b,x,y) = a*x+b*y
 * we can often omit the return value */
int xgcd(int a, int b, int &x, int &y) {
    if (a<b) { swap(a,b); swap(x,y); }
    if (!b) { x=1, y=0; return a; }
    int g = xgcd(b,a%b,y,x);
    y -= (a/b)*x;
    return g;
}

int modinv(int a) {
    /* sol 1 */
    return modpow(a, MOD-2, MOD);
    /* sol 2
    int x,y;
    xgcd(MOD, a, y, x);
    return y;
    */
}


/* Chinese Remainder Theorem
 * backstory:
 * a chinese general could count his soldiers crazy fast
 * by looking at the remainder of the integer division
 * of his troups in various group sizes
 * math:
 * find X such that
 * X = a_k [k], for pairwise coprime ki
 * the solution is unique modulo k1*..*kn
 * warning: be aware of overflows */

pii CRT(vector<pii>& ak) {
    ll a=0,mod=1;
    for (auto& k:ak) {
        int x,y; xgcd(mod,k.yy,x,y);
        /* x*mod = 1[kyy], y*kyy = 1[mod] */
        a = k.xx*x*mod + a*y*k.yy;
        mod *= k.yy;
        a %= mod;
    }
    if (a<0) a+=mod;
    return {a,mod};
}


int main() {
    string lb = "longboiiiiiiiiiiiiii";
    printf("hash of %s: %d\n", lb.c_str(), p_hash(lb));
    int x = 3;
    printf("%d * %d = 1[%d]\n",x,modinv(3),MOD);
    vector<pii> ak = {{1,3},{4,5},{2,7}};
    int a,m; tie(a,m) = CRT(ak);
    printf("CRT: %d [%d]\n",a,m);
}

