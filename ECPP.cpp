// uses NTL
//   http://www.shoup.net/ntl

#include "ECPP.h"

long ECPP_TRYDIV_MAX = 1<<16;
double ECPP_RHO_TIMEOUT = 0.01;// for Pollard rho method

long ProbPrimeFactor(ZZ& d, const ZZ& n, const ZZ& t)
// d = probable prime factor of n such that t<d<n
// return 0 or 1 if d is found or not, respectively
// return 1 if n is probable prime
{
    ZZ q;
    long p;
    PrimeSeq ps;
    if(ProbPrime(d=n)) return 1;
    while((p = ps.next()) && p < ECPP_TRYDIV_MAX) {
        if(!divide(d,d,p)) continue;
        while(divide(d,d,p));
        if(d<=t) return 1;
        if(ProbPrime(d)) return 0;
    }
    for(;;) {
        if(brent_rho(d, q=d, ECPP_RHO_TIMEOUT))
            return 1;
        if(d<=t) {
            div(d,q,d);
            if(d<=t) return 1;
        }
        if(ProbPrime(d)) return 0;
    }    
}

void set(CertElem& e, const EC_p& P, const ZZ& m1) {
    e.P = P;
    e.m = m1;
    e.n = ZZ_p::modulus();
    e.a = P.a();
    e.b = P.b();
}

long Certify(const ECPPCert& c)
// return 1 if certificate c is correct
// return 0 if c is not correct
{
    long i;
    ZZ t,q1;
    EC_p P;
    if(c.length()==0) return 0;
    for(i=c.length()-1; i>0; i--) {
        if(!IsOdd(c[i].n)) return 0;
        if(c[i].n % 3 ==0) return 0;// gcd(n,6)==1
        ZZ_p::init(c[i].n);
        EC_p::init(c[i].a, c[i].b);
        if(EC_p::IsSingular()) return 0;
        if(!IsOnEC(c[i].P)) return 0;
        SqrRoot(q1, c[i].n); LeftShift(t,q1,2);
        SqrRoot(t,t);
        q1+=t; q1++;// q1 = (n^{1/4}+1)^2
        if(c[i-1].n <= q1) return 0;/// q>q1
        mul(P, c[i].P, c[i].m);
        if(!IsZero(P)) return 0;// mP==O
        div(t, c[i].m, c[i-1].n);
        mul(P, c[i].P, t);
        if(!z_coprime(P)) return 0;// (m/q)P!=O
    }
    conv(i, c[0].n);
    return IsPrime(i);
}