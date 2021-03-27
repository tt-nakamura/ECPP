// uses NTL
//   http://www.shoup.net/ntl
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    Algorithm 1.7.3

#include<NTL/ZZ.h>
using namespace NTL;

static bool q64[] = {
    1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,
    1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
};
static bool q63[] = {
    1,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,
    1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,
    0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,
    0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,
};
static bool q65[] = {
    1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,
    1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,
    0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,
    0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,
    1,
};
static bool q11[] = {
    1,1,0,1,1,1,0,0,0,1,0,
};

long IsSquare(ZZ& q, const ZZ& n)
// if n is square, set q=[sqrt(n)] and return 1
// assume n>0
{
    static long m(63*65*11);
    if(!q64[trunc_long(n,6)]) return 0;
    long r(n%m);
    if(!q63[r%63]) return 0;
    if(!q65[r%65]) return 0;
    if(!q11[r%11]) return 0;
    if(&q!=&n) {
        ZZ s;
        SqrRoot(q,n);
        sqr(s,q);
        return s==n;
    }
    else {
        ZZ s,m(n);
        SqrRoot(q,m);
        sqr(s,q);
        return s==m;
    }
}

long IsSquare(long n)
// if n is square, return [sqrt(n)] else return 0
// assume n>0
{
    if(!q64[n&63]) return 0;
    if(!q63[n%63]) return 0;
    if(!q65[n%65]) return 0;
    if(!q11[n%11]) return 0;
    long s = SqrRoot(n);
    if(s*s==n) return s;
    else return 0;
}