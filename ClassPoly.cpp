// uses NTL
//   http://www.shoup.net/ntl
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    Algorithm 7.6.1

#include<NTL/ZZX.h>
#include "RRX.h"
#include "CC.h"
using namespace NTL;

void Delta_q(CC& d, const CC& q)
// d = q*( 1 + sum_{n=1}^infty (-1)^n
//      *( q^{n(3n-1)/2} + q^{n(3n+1)/2} ) )^24
{
    long i;
    CC z1,z2,q1,q2,q3,d1;
    sqr(q2,q);
    q1=z1=q;
    q3=z2=q2;
    q3*=q;
    set(d);
    d-=z1;
    d-=z2;
    for(i=0; d!=d1; i^=1) {
        q1*=q3;
        q2*=q3;
        z1*=q1;
        z2*=q2;
        d1=d;
        if(i&1) { d-=z1; d-=z2; }
        else    { d+=z1; d+=z2; }
    }
    pow(d,d1,24);
    d*=q;
}

void j_inv(CC& j, const CC& q)
// j = (256*f + 1)^3 / f
//     where f = Delta(q^2) / Delta(q)
{
    CC z,w;
    sqr(z,q);
    Delta_q(w,z);
    Delta_q(z,q);
    w/=z;
    z=w;
    z*=256.;
    z+=1.;
    sqr(j,z);
    j*=z;
    j/=w;
}

long ClassPoly(ZZX& f, long d)
// d<0, d = fundamental discriminant of quadratic field
// f = Hilbert class polynomial, return class number
{
    static double PL2(atan(1)*4/log(2));
    double sum(0), rd(sqrt(-d));
    long i,a,b,ac,a2,h(0),bmax(long(rd/sqrt(3)));
    CC j;
    RR RD,PI;
    RRX u,v,w;
    u.SetLength(2);
    v.SetLength(3);
    w.SetLength(1);
    for(b=(d&1); b<=bmax; b+=2) {
        ac = (b*b-d)>>2;
        for(a=max(b,1); (a2=a*a)<=ac; a++) {
            if(ac%a) continue;
            if(a==b || a2==ac || b==0) sum += 1./a;
            else sum += 2./a;
        }
    }
    RR::SetPrecision(long(PL2*rd*sum) + 48);
    conv(RD,-d);
    SqrRoot(RD,RD);
    ComputePi(PI);
    set(u[1]);
    set(v[2]);
    set(w[0]);
    for(b=(d&1); b<=bmax; b+=2) {
        ac = (b*b-d)>>2;
        for(a=max(b,1); (a2=a*a)<=ac; a++) {
            if(ac%a) continue;
            negate(real(j), RD);
            conv(imag(j), -b);
            j /= double(a);
            j *= PI;
            exp(j,j);
            j_inv(j,j);
            if(a==b || a2==ac || b==0) {
                h++;
                negate(u[0], real(j));
                w *= u;
            }
            else {
                h+=2;
                negate(v[1], real(j));
                v[1] *= 2.;
                norm(v[0], j);
                w *= v;
            }
        }
    }
    f.SetLength(h+1);
    for(i=h; i>=0; i--) RoundToZZ(f[i], w[i]);
    return h;
}