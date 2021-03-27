// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
using namespace NTL;

#define RHO_GCD_INTVL 100

long pollard_rho(ZZ& d, const ZZ& n, double T)
// input:
//   n = composite integer, n>=4
//   T = timeout in seconds
// output:
//   d = divisor of n, 1 < d < n
//       by Pollard rho method
// return:
//   0 if successful, -1 if failure
// reference:
//   J. M. Pollard "A Monte Carlo Method for Factorization"
//     BIT Numerical Mathematics 15 (1975) 331
{
    long a(1),j;
    ZZ u(2),v(2),q(1),s,t;

    if(&d==&n) return pollard_rho(d,s=n,T);
    T += GetTime();
    for(;;) {
        s=u; t=v;
        for(j=0; j<RHO_GCD_INTVL; j++) {
            SqrMod(u,u,n); AddMod(u,u,a,n);
            SqrMod(v,v,n); AddMod(v,v,a,n);
            SqrMod(v,v,n); AddMod(v,v,a,n);
            SubMod(d,u,v,n);
            MulMod(q,q,d,n);
        }
        GCD(d,q,n);
        if(!IsOne(d)) break;
        if(GetTime() > T) return -1;
    }
    if(d==n) do {
        SqrMod(s,s,n); AddMod(s,s,a,n);
        SqrMod(t,t,n); AddMod(t,t,a,n);
        SqrMod(t,t,n); AddMod(t,t,a,n);
        sub(q,s,t);
        GCD(d,q,n);
    } while(IsOne(d));
    if(IsOne(d) || d==n) return 1;
    else return 0;
}

long brent_rho(ZZ& d, const ZZ& n, double T)
// input:
//   n = composite integer, n>=4
//   T = timeout in seconds
// output:
//   d = divisor of n, 1 < d < n
//       by modified Pollard rho method
// return:
//   0 if successful, -1 if failure
// reference:
//   R. P. Brent "An Improved Monte Carlo Factorization Algorithm"
//     BIT Numerical Mathematics 20 (1980) 176
{
    ZZ u(2),q,s,t;
    long a,r,i,j;
    
    if(&d==&n) return brent_rho(d,s=n,T);
    T += GetTime();
    for(a=1;; a++) {
        set(q);
        for(r=1; r>0; r<<=1) {
            s=u;
            for(i=0; i<r; i++) {
                SqrMod(u,u,n);
                AddMod(u,u,a,n);
            }
            for(i=j=0; i<r;) {
                t=u;
                j += RHO_GCD_INTVL;
                if(j>r) j=r;
                for(; i<j; i++) {
                    SqrMod(u,u,n);
                    AddMod(u,u,a,n);
                    SubMod(d,s,u,n);
                    MulMod(q,q,d,n);
                }
                GCD(d,q,n);
                if(!IsOne(d)) goto a;
                if(GetTime() > T) return -1;
            }
        }
    a:  ;
        if(d<n) return 0;
        do {
            SqrMod(t,t,n);
            AddMod(t,t,a,n);
            sub(q,s,t);
            GCD(d,q,n);
        } while(IsOne(d));
        if(d<n) return 0;
    }
}