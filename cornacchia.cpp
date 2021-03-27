// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
using namespace NTL;

long IsSquare(ZZ&, const ZZ&);

long cornacchia(ZZ& x, ZZ& y, const ZZ& d, const ZZ& p)
// solve x^2-dy^2=4p by Cornacchia method
//   where p=prime, -4p<d<0, d=0,1(mod4)
// if solution is found, return 0, else return 1
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    Algorithm 1.5.3
{
    ZZ q,r,t;
    if(p==2) {
        add(t,d,8);
        if(!IsSquare(x,t)) return 1;
        set(y);
        return 0;
    }
    rem(t,d,p);
    if(Jacobi(t,p)<0) return 1;
    SqrRootMod(x,t,p);
    if(IsOdd(x)^IsOdd(d)) sub(x,p,x);
    LeftShift(y,p,1);
    LeftShift(q,y,1);
    SqrRoot(t,q);
    while(x>t) {
        rem(r,y,x);
        y=x;
        x=r;
    }
    sqr(t,x);
    t-=q;
    if(!divide(t,t,d)) return 1;
    if(!IsSquare(y,t)) return 1;
    return 0;
}