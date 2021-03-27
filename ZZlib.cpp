// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

long IsPrime(long n)
// primality proving by trial division
// return 1 or 0 if n is prime or not, respectively
// assume n>=2
{
    long p,m;
    if(n<2) return 0;
    m = SqrRoot(n);
    PrimeSeq ps;
    while((p = ps.next()) && p<=m)
        if(n%p == 0) return 0;
    if(p) return 1;
    else Error("too large n");
}

long Jacobi(long a, long b) 
// Jacobi Symbol
// 0<=a<b, b odd
// if((a%=b) < 0) a+=b;
{
	long c,t(0);
	while(a) {
		for(c=0; (a&1)==0; c++) a>>=1;
		if(c&1) {
			c = b&7;
			if(c==3 || c==5) t ^= 1;
		}
		if(a&b&2) t ^= 1;
		c = b%a;
		b = a;
		a = c;
	}
	if(b!=1) return 0;
	else if(t) return -1;
	else return 1;
}

void factor(Vec<Pair<long, long> >& f, long n)
// input:
//   n = integer, |n| < 2^{60}
// output:
//   f = prime factorization of |n| by trial division
//       vector of (prime, exponent) pair
//       in increasing order of primes
{
    long i(0),j,m,p;
    
    if(n<0) n=-n;
    if(n==0 || n==1) {
        f.SetLength(0);
        return;
    }
    m = SqrRoot(n);
    PrimeSeq ps;
    while((p = ps.next()) <= m) {
        if(p==0) Error("too large n");
        for(j=0; n%p==0; j++) n/=p;
        if(j==0) continue;
        f.SetLength(i+1);
        f[i].a = p;
        f[i].b = j;
        if(n==1) return;
        i++;
        m = SqrRoot(n);
    }
    f.SetLength(i+1);
    f[i].a = n;
    f[i].b = 1;
}

void CRT(long& x, long& m, long a, long p)
// incremental chinese remaindering
// compute y such that y=x (mod m), y=a (mod p), 0<=y<mp
//   and set x:=y, m:=mp
// assume 0<=a<p, p is coprime to m
{
    long b(m);
    a = SubMod(a, x%p, p);
    b *= MulMod(a, InvMod(m%p, p), p);
    x += b;
    m *= p;
}