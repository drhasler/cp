/* git@drhasler - 02/19 */
#include <bits/stdc++.h>
using namespace std;

#define cd complex<double>
const double PI = acos(-1);

/* Fast Fourier Transform
 * O(nlogn)
 * FFT[j] = sum x[k] e^{ik2pi/N}
 * iFFT[j] = 1/N sum X[k] e^{-ik2pi/N}
 *
 * recursion (divide and conquer):
 * - split odd and even coefficients
 * - fft each array
 * - merge
 * FFT_n[  i  ] = E[i] + e^{i2pi/n} O[i]
 * FFT_n[i+n/2] = E[i] - e^{i2pi/n} O[i]
 *
 * applications:
 * 1. polynomial multiplication
 * store coefficients in an array
 * A*B = invFFT(FFT(A).*FFT(B))
 * (works integer multiplication as well)
 */


// Recursive FFT
void recfft(vector<cd>& a,bool invert) {
    int n = a.size();
    int n2 = n/2;
    vector<cd> a0(n2),a1(n2);
    for (int i=0;i<n2;i++) {
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }
    if (n2>1) {
        recfft(a0,invert);
        recfft(a1,invert);
    }
    double ang = PI/n2*(invert?-1:1);
    cd w(1), wn(cos(ang),sin(ang));
    for (int i=0;i<n2;i++) {
        a[i]    = a0[i] + w*a1[i];
        a[i+n2] = a0[i] - w*a1[i];
        w *= wn;
    }
    if (invert) for (int i=0;i<n;i++) a[i]/=2;
}

/* in-place: Fast FFT
 * each index bitwise backward
 * such that bigger clusters for lower significant bit */
int reverse(int num,int logn) {
    int res = 0;
    for (int i=0;i<logn;i++) if (num&(1<<i)) res |= 1<<(logn-1-i);
    return res;
}
void ipfft(cd* a,int n,bool invert) {
    int logn = 0;
    while ((1 << logn) < n) logn++;
    for (int i=0,r;i<n;i++) if (i < (r=reverse(i,logn)))
        swap(a[i],a[r]);
    for (int len=1;len<n;len*=2) {
        double ang = PI/len * (invert ? -1 : 1);
        cd wlen(cos(ang),sin(ang));
        for (int i=0;i<n;i+=2*len) {
            cd w(1);
            for (int j=0;j<len;j++) {
                cd u = a[i+j], v = a[i+j+len]*w;
                a[i+j] = u+v;
                a[i+j+len] = u-v;
                w *= wlen;
            }
        }
    }
    if (invert) for (int i=0;i<n;i++) a[i]/=n;
}

/* Transformation Fourier Plus Rapide Que Rapide
 * computes the reverse more efficiently
 * by simulating all the swaps */
void tfprqr(cd* a,int n, bool invert) {
    for (int i=1,j=0;i<n;i++) {
        int bit = n/2;
        for (;j&bit;bit/=2) j^=bit;
        j^=bit;
        if (i<j) swap(a[i],a[j]);
    }
    for (int len=1;len<n;len*=2) {
        double ang = PI/len * (invert ? -1 : 1);
        cd wlen(cos(ang),sin(ang));
        for (int i=0;i<n;i+=2*len) {
            cd w(1);
            for (int j=0;j<len;j++) {
                cd u = a[i+j], v = a[i+j+len]*w;
                a[i+j] = u+v;
                a[i+j+len] = u-v;
                w*=wlen;
            }
        }
    }
    if (invert) for (int i=0;i<n;i++) a[i]/=n;
}

int main() {
    int n = 8;
    cd a[8];
    cd b[8];
    a[2] = b[3] = 2.;
    tfprqr(a,8,0);
    tfprqr(b,8,0);
    for (int i=0;i<n;i++) a[i] *= b[i];
    tfprqr(a,8,1);
    for (int i=0;i<n;i++) cout << a[i].real() << ' ';
    cout << endl;
}
