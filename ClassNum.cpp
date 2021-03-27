#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

void factor(Vec<Pair<long, long> >&, long);

long conductor(long d)
// input:
//   d = discriminant, d=0 or 1 (mod 4)
// return:
//   largest f such that f^2 divides d or d/4
//   d/4 if s odd or t=3 (mod 4), where d = 2^s t (t odd)
//   d otherwise
{
    long i,j(d),f(1);
    Vec<Pair<long, long> > p;
    factor(p,-d);
    for(i=0; (j&1)==0; i++) j>>=1;
    if(i&1 || j&2) p[0].b -= 2;
    for(i=0; i<p.length(); i++) {
        if(p[i].b >= 2)
            f *= power_long(p[i].a, p[i].b>>1);
    }
    return f;
}

long ClassNum(long d)
// input: discriminant d<0 of imaginary quadratic fields
// return: class number h
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    Algorithm 5.3.5
{
    long a,b,f,q,a2,h(0),B(sqrt(-d/3.));
    for(b=(d&1); b<=B; b+=2) {
        q = b*b-d; q>>=2;
        for(a=max(b,1); (a2=a*a)<=q; a++) {
            if(q%a) continue;
            if(GCD(GCD(a,b),q/a)>1) continue;
            if(a==b || a2==q || b==0) h++;
            else h+=2;
        }        
    }
    return h;
}