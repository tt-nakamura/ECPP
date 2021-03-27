// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   R. Crandall and C. Pomerance
//    "Prime Numbers: A Computational Perspective"
//     2nd edition, section 7.5.2

#include "EC_pE.h"

Vec<ZZ_pE> EC_pE::ff;
Vec<ZZ_pX> EC_pE::g;

void EC_pE::init(long n)
// set degree n of torsion points so that [n]P = O for all P
// division polynomials are unchanged from previous init
//   but are appended if n is increased
// assume EC_p has been initialized
{
    ZZ_pX h;
    DivPoly(n+1);
    div(h, g[n], n);
    ZZ_pE::init(h);
    ff.SetLength(2);
    conv(ff[0], EC_p::f);
    sqr(ff[1], ff[0]);
}

void set(EC_pE& a)
// set a=(X,1)
{
    ZZ_pX X;
    SetX(X);
    conv(a.x, X);
    set(a.y);
    set(a.z);
}

void set(EC_pE& a, long n)
// set a = [n](X,1)
// assume |n| <= g.length()-3
// reference: Crandall & Pomerance, Theorem 7.5.5
{
    bool s(n<0);
    if(s) n=-n;
    if(n==0) clear(a);
    else if(n==1) set(a);
    else {
        ZZ_pE p1,p2,q1,q2;
        ZZ_pX X; SetX(X);
        conv(p1, EC_pE::g[n]);
        sqr(p2,p1);
        mul(a.z, p1, p2);
        if(!(n&1)) a.z *= EC_pE::f2();
        conv(q1,X);
        mul(a.x, q1, a.z);
        conv(q1, EC_pE::g[n-1]);
        conv(q2, EC_pE::g[n+1]);
        mul(p2,q1,q2);
        p2 *= p1;
        p2 *= EC_pE::f();
        a.x -= p2;
        sqr(q1,q1);
        sqr(q2,q2);
        conv(p1, EC_pE::g[n+2]);
        conv(p2, EC_pE::g[n-2]);
        q1 *= p1;
        q2 *= p2;
        sub(a.y, q1, q2);
        a.y /= 4;
    }
    if(s) negate(a,a);
}

void negate(EC_pE& p, const EC_pE& a) {// p=-a
    negate(p.y, a.y);
    if(&p!=&a) { p.x = a.x; p.z = a.z; }
}

void doubleh(EC_pE& p, const EC_pE& a)
// set p=a+a
// reference: Crandall & Pomerance, p326
{
    if(IsZero(a.y) || IsZero(a.z)) clear(p);
    else {
        ZZ_pE l,m,n,s;
        sqr(s, a.y);
        mul(l, s, a.x);
        l *= a.z;
        l *= EC_pE::f();
        l *= 4;
        sqr(m, a.x);
        sqr(n, a.z);
        m *= 3;
        n *= EC_p::a();
        m += n;
        mul(n, a.y, a.z);
        n += n;
        sqr(p.x, m);
        p.x -= l;
        p.x -= l;
        l -= p.x;
        p.x *= n;
        p.x *= EC_pE::f();
        mul(p.y, m, l);
        sqr(m,n);
        s *= m;
        s *= EC_pE::f2();
        p.y -= s;
        p.y -= s;
        mul(p.z, m, n);
        p.z *= EC_pE::f2();
    }
}

void add(EC_pE& p, const EC_pE& a, const EC_pE& b)
// set p=a+b
// reference: Crandall & Pomerance, p326
{
    if(IsZero(a)) p=b;
    else if(IsZero(b)) p=a;
    else {
        ZZ_pE s,t,u,v,w,m;
        mul(u, b.x, a.z);
        mul(v, a.x, b.z);
        mul(w, b.y, a.z);
        mul(m, a.y, b.z);
        sub(s,u,v);
        add(t,u,v);
        sub(u,w,m);
        if(IsZero(s)) {
            if(IsZero(u)) doubleh(p,a);
            else clear(p);
            return;
        }
        add(v,w,m);
        mul(w, a.z, b.z);
        sqr(m,s);
        t *= m;
        sqr(p.x, u);
        p.x *= w;
        p.x *= EC_pE::f();
        p.x -= t;
        t -= p.x;
        t -= p.x;
        p.x *= s;
        mul(p.y, u, t);
        m *= s;
        v *= m;
        p.y -= v;
        p.y /= 2;
        mul(p.z, m, w);
    }
}

void sub(EC_pE& p, const EC_pE& a, const EC_pE& b)
// set p=a-b
{
    if(&a==&b) clear(p);
    else if(&p==&a) {
        EC_pE c;
        negate(c,b);
        add(p,a,c);
    }
    else {
        negate(p,b);
        add(p,a,p);
    }
}

void mul(EC_pE& p, const EC_pE& a, long n)
// set p = [n]a
{
    bool s(n<0);
    if(s) n=-n;
    if(n==0) clear(p);
    else if(n==1) p=a;
    else if(n==2) doubleh(p,a);
    else if(&p==&a) {
        EC_pE b(a);
        mul(p,b,n);
    }
    else {
        unsigned long k,m(n);
        m <<= 1;
        m += n;
        k = (1L<<(NumBits(m)-1));
        p = a;
        for(k>>=1; k>1; k>>=1) {
            doubleh(p,p);
            if((m&k) && !(n&k)) add(p,p,a);
            else if(!(m&k) && (n&k)) sub(p,p,a);
        }
    }
    if(s) negate(p,p);
}

long compare_x(const EC_pE& a, const EC_pE& b)
// return (a.x/a.z == b.x/b.z)
{
    ZZ_pE s1,s2;
    mul(s1, a.z, b.x);
    mul(s2, b.z, a.x);
    return s1 == s2;
}

long compare_y(const EC_pE& a, const EC_pE& b)
// return (a.y/a.z == b.y/b.z)
{
    ZZ_pE s1,s2;
    mul(s1, a.z, b.y);
    mul(s2, b.z, a.y);
    return s1 == s2;
}

void FrobConj(EC_pE& p, const EC_pE& a)
// set p=a^p (Frobenius Conjugate)
{
    ZZ p1;
    ZZ_pE t;
    power(p.x, a.x, ZZ_p::modulus());
    power(p.y, a.y, ZZ_p::modulus());
    power(p.z, a.z, ZZ_p::modulus());
    sub(p1, ZZ_p::modulus(), 1); p1>>=1;
    power(t, EC_pE::f(), p1);
    p.y *= t;
}

void EC_pE::DivPoly(long n)
// construct division polynomials g
// append g[i] (i=m,...,n) where m=g.length
// reference: Crandall & Pomerance, Definition 7.5.4
{
    long m(g.length());
    if(n>=m) g.SetLength(n+1);
    if(m<=1 && n>=1) set(g[1]);
    if(m<=2 && n>=2) SetCoeff(g[2], 0, 2);
    if(m<=3 && n>=3) {
        SetCoeff(g[3], 4, 3);
        mul(g[3][2], EC_p::a(), 6);
        mul(g[3][1], EC_p::b(), 12);
        sqr(g[3][0], EC_p::a());
        negate(g[3][0], g[3][0]);
    }
    if(m<=4 && n>=4) {
        ZZ_p t;
        SetCoeff(g[4], 6, 1);
        mul(g[4][4], EC_p::a(), 5);
        mul(g[4][3], EC_p::b(), 20);
        sqr(g[4][2], EC_p::a());
        g[4][2] *= -5;
        mul(g[4][1], EC_p::a(), EC_p::b());
        g[4][1] *= -4;
        sqr(g[4][0], EC_p::b());
        g[4][0] *= -8;
        power(t, EC_p::a(), 3);
        g[4][0] -= t;
        g[4] *= 4;
    }
    long i,j;
    ZZ_pX t,f;
    sqr(f, EC_p::f);
    for(i=m; i<=n; i++) {
        if(i<5) continue;
        j = (i>>1);
        if(i&1) {
            power(g[i], g[j], 3);
            g[i] *= g[j+2];
            power(t, g[j+1], 3);
            t *= g[j-1];
            if(j&1) t *= f;
            else g[i] *= f;
            g[i] -= t;
        }
        else {
            sqr(g[i], g[j-1]);
            g[i] *= g[j+2];
            sqr(t, g[j+1]);
            t *= g[j-2];
            g[i] -= t;
            g[i] *= g[j];
            g[i] /= 2;
        }
    }
}