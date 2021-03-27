// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ_pXFactoring.h>
#include "EC_pE.h"

long EC_p::schoof(ZZ& o)
// set o = number of points on Elliptic Curve
//   by Schoof method
// reference: R.Crandall and C.Pomenrace
//   "Prime Numbers: A Computational Perspective"
//   second edition, Algorithm 7.5.6
{
    const ZZ& p(ZZ_p::modulus());// large prime
    long j,k,l,q;
    EC_pE P,Q,R;
    ZZ M,m,t,p1;
    add(o,p,1);
    sub(p1,p,1); p1>>=1;
    LeftShift(M,p,4);
    SqrRoot(M,M);// M=4*sqrt(p)
    PrimeSeq ps;
    m = ps.next();
    t = DetIrredTest(EC_p::f);//t%2 in CRT
    while(m<=M) {
        l = ps.next();
        if((q = p%l) == 0) break;// p is not prime
        EC_pE::init(l);
        set(P);
//      FrobConj(P,P);
//      FrobConj(Q,P);// Q=Phi(Phi(P))
        power(P.x, P.x, p);
        power(P.y, EC_pE::f(), p1);
        power(Q.x, P.x, p);
        power(Q.y, P.y, o);// Q=Phi(Phi(P))
        set(Q.z);
        set(R,q);// R=[p%l]P
        add(Q,Q,R);
        clear(R);
        for(j=0, k=(l>>1); j<=k; j++) {
            if(compare_x(Q,R)) {
                if(!compare_y(Q,R)) j=l-j;
                break;
            }
            add(R,R,P);//R=[j]P
        }
        CRT(t,m,j,l);
    }
    o -= t;// o=p+1-t
    // kill to avoid "internal error: can't grow this _ntl_gbigint"
    EC_pE::kill();
    return !q;
}