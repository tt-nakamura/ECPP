// uses NTL
//   http://www.shoup.net/ntl

#ifndef __Poly_H__
#define __Poly_H__

#include<NTL/vector.h>
using namespace NTL;

template<class T>
struct Poly {// polynomial
    Vec<T> c;// coefficients
    inline long length() const { return c.length(); }
    inline void SetLength(long i) { c.SetLength(i); }
    inline T& operator[](long i) { return c[i]; }
    inline const T& operator[](long i) const { return c[i]; }
    void operator=(const T&);
    void normalize();
};

template<class T> long operator==(const Poly<T>&, const Poly<T>&);
template<class T> long operator==(const Poly<T>&, const T&);
template<class T> void operator*=(Poly<T>&, const Poly<T>&);
template<class T> void operator%=(Poly<T>&, const Poly<T>&);
template<class T> void operator/=(Poly<T>&, const Poly<T>&);

template<class T> void add(Poly<T>&, const Poly<T>&, const T&);
template<class T> void sub(Poly<T>&, const Poly<T>&, const T&);
template<class T> void add(Poly<T>&, const Poly<T>&, const Poly<T>&);
template<class T> void sub(Poly<T>&, const Poly<T>&, const Poly<T>&);
template<class T> void mul(Poly<T>&, const Poly<T>&, const Poly<T>&);
template<class T> void mul(Poly<T>&, const Poly<T>&, const T&);
template<class T> void div(Poly<T>&, const Poly<T>&, const T&);
template<class T> void rem(Poly<T>&, const Poly<T>&, const T&);
template<class T> void sqr(Poly<T>&, const Poly<T>&);

template<class T> void negate(Poly<T>&, const Poly<T>&);
template<class T> void DivRem(Poly<T>&, Poly<T>&, const Poly<T>&, const Poly<T>&);
template<class T> void PseudoDivRem(Poly<T>&, Poly<T>&, const Poly<T>&, const Poly<T>&);
template<class T> void translate(Poly<T>&, const Poly<T>&, const T&);
template<class T> void diff(Poly<T>&, const Poly<T>&);
template<class T> void MakeMonic(Poly<T>&);
template<class T> void LeftShift(Poly<T>&, const Poly<T>&, long);
template<class T> void RightShift(Poly<T>&, const Poly<T>&, long);
template<class T> void scale(Poly<T>&, const Poly<T>&, const T&);
template<class T> void eval(T&, const Poly<T>&, const T&);

template<class T> inline void operator+=(Poly<T>& x, const T& a) { add(x,x,a); }
template<class T> inline void operator-=(Poly<T>& x, const T& a) { sub(x,x,a); }
template<class T> inline void operator+=(Poly<T>& x, const Poly<T>& a) { add(x,x,a); }
template<class T> inline void operator-=(Poly<T>& x, const Poly<T>& a) { sub(x,x,a); }
template<class T> inline void operator*=(Poly<T>& x, const Poly<T>& a) { mul(x,x,a); }
template<class T> inline void operator*=(Poly<T>& x, const T& a) { mul(x,x,a); }
template<class T> inline void operator/=(Poly<T>& x, const T& a) { div(x,x,a); }
template<class T> inline void operator%=(Poly<T>& x, const T& a) { rem(x,x,a); }
template<class T> inline void operator<<=(Poly<T>& x, long k) { LeftShift(x,x,k); }
template<class T> inline void operator>>=(Poly<T>& x, long k) { RightShift(x,x,k); }
template<class T> inline long operator!=(const Poly<T>& a, const Poly<T>& b) { return !(a==b); }
template<class T> inline long operator!=(const Poly<T>& a, const T& b) { return !(a==b); }

template<class T> inline void mul(Poly<T>& x, const T& a, const Poly<T>& b) { mul(x,b,a); }
template<class T> inline void clear(Poly<T>& a) { a.c.SetLength(0); }
template<class T> inline void set(Poly<T>& a) { a.c.SetLength(1), set(a.c[0]); }
template<class T> inline long deg(const Poly<T>& a) { return a.c.length()-1; }
template<class T> inline long IsZero(const Poly<T>& a) { return a.c.length()==0; }
template<class T> inline long IsOne(const Poly<T>& a) { return a.c.length()==1 && IsOne(a.c[0]); }
template<class T> inline long IsMonic(const Poly<T>& a) { return IsOne(a.c[a.c.length()-1]); }
template<class T> inline const T& LeadCoeff(const Poly<T>& a) { return a.c[a.c.length()-1]; }
template<class T> inline T& LeadCoeff(Poly<T>& a) { return a.c[a.c.length()-1]; }
template<class T> inline std::ostream& operator<<(std::ostream& s, const Poly<T>& a) { return s << a.c; }

template<class T>
void Poly<T>::operator=(const T& a) {
    if(IsZero(a)) c.SetLength(0);
    else { c.SetLength(1); c[0] = a; }
}

template<class T>
void add(Poly<T>& x, const Poly<T>& a, const T& b) {
    if(!IsZero(a)) {
        if(&x!=&a) {
            x.SetLength(a.length());
            for(long i=1; i<a.length(); i++) x[i] = a[i];
        }
        add(x[0], a[0], b);
        if(x.length()==1 && IsZero(x[0])) clear(x);
    }
    else if(!IsZero(b)) {
        x.SetLength(1);
        x[0] = b;
    }
    else if(&x!=&a) clear(x);
}

template<class T>
void sub(Poly<T>& x, const Poly<T>& a, const T& b) {
    if(!IsZero(a)) {
        if(&x!=&a) {
            x.SetLength(a.length());
            for(long i=1; i<a.length(); i++) x[i] = a[i];
        }
        sub(x[0], a[0], b);
        if(x.length()==1 && IsZero(x[0])) clear(x);
    }
    else if(!IsZero(b)) {
        x.SetLength(1);
        negate(x[0], b);
    }
    else if(&x!=&a) clear(x);
}

template<class T>
void add(Poly<T>& x, const Poly<T>& a, const Poly<T>& b) {
    long i, m(a.length()), n(b.length());
    if(m>=n) {
        if(&x!=&a) {
            x.SetLength(m);
            for(i=n; i<m; i++) x[i] = a[i];
        }
        for(i=0; i<n; i++) add(x[i], a[i], b[i]);
        if(m==n) x.normalize();
    }
    else {
        if(&x!=&b) {
            x.SetLength(n);
            for(i=m; i<n; i++) x[i] = b[i];
        }
        for(i=0; i<m; i++) add(x[i], a[i], b[i]);
    }
}

template<class T>
void sub(Poly<T>& x, const Poly<T>& a, const Poly<T>& b) {
    long i, m(a.length()), n(b.length());
    if(m>=n) {
        if(&x!=&a) {
            x.SetLength(m);
            for(i=n; i<m; i++) x[i] = a[i];
        }
        for(i=0; i<n; i++) sub(x[i], a[i], b[i]);
        if(m==n) x.normalize();
    }
    else {
        if(&x!=&b) x.SetLength(n);
        for(i=m; i<n; i++) negate(x[i], b[i]);
        for(i=0; i<m; i++) sub(x[i], a[i], b[i]);
    }
}

template<class T>
void mul(Poly<T>& x, const Poly<T>& a, const T& b) {
    if(&x!=&a) x.SetLength(a.length());
    for(long i=0; i<a.length(); i++) mul(x[i], a[i], b);
    x.normalize();
}

template<class T>
void div(Poly<T>& x, const Poly<T>& a, const T& b) {
    if(&x!=&a) x.SetLength(a.length());
    for(long i=0; i<a.length(); i++) div(x[i], a[i], b);
}

template<class T>
void rem(Poly<T>& x, const Poly<T>& a, const T& b) {
    if(&x!=&a) x.SetLength(a.length());
    for(long i=0; i<a.length(); i++) rem(x[i], a[i], b);
    x.normalize();
}

template<class T>
void operator%=(Poly<T>& x, const Poly<T>& a) {
    if(&x==&a) { clear(x); return; }
    static T t,r;
    long i,j,k,m(!IsMonic(a));
//  if(m) inv(r, LeadCoeff(a));
    i = x.length() - a.length();
    k = x.length() - 1;
    while(i>=0) {
//      if(m) x[k] *= r;
        if(m) x[k] /= LeadCoeff(a);
        j = a.length() - 2;
        while(j>=0) {
            mul(t, x[k], a[j]);
            x[i+j] -= t;
            j--;
        }
        i--; k--;
    }
    if(x.length() >= a.length()) {
        x.SetLength(a.length() - 1);
        x.normalize();
    }
}

template<class T>
void operator/=(Poly<T>& x, const Poly<T>& a) {
    if(&x==&a) set(x);
    else if(x.length() < a.length()) clear(x);
    else {
        long n(x.length()), m(a.length()-1), i;
        x%=a;
        x.SetLength(n-m);
        for(i=0; i<x.length(); i++) x[i] = x[m+i];
    }
}

template<class T>
long operator==(const Poly<T>& a, const Poly<T>& b) {
    if(a.length() != b.length()) return false;
    for(long i=0; i<a.length(); i++)
        if(a[i] != b[i]) return false;
    return true;
}

template<class T>
long operator==(const Poly<T>& a, const T& b) {
    if(a.length() == 0) return IsZero(b);
    else return a.length()==1 && a[0] == b;
}

template<class T>
void Poly<T>::normalize() {
    long i(c.length()-1);
    if(i<0 || !IsZero(c[i])) return;
    for(i--; i>=0; i--)
        if(!IsZero(c[i])) break;
    c.SetLength(i+1);
}

template<class T>
void mul(Poly<T>& x, const Poly<T>& a, const Poly<T>& b) {
    if(IsZero(a) || IsZero(b)) clear(x);
    else if(&a==&b) sqr(x,a);
    else if(&x==&b) mul(x,b,a);
    else {
        long i,j;
        static T t;
        i = a.length()-1;
        x.SetLength(i + b.length());
        for(j=b.length()-1; j>=0; j--) mul(x[i+j], a[i], b[j]);
        for(i--; i>=0; i--) {
            for(j=b.length()-1; j>0; j--) {
                mul(t, a[i], b[j]);
                x[i+j] += t;
            }
            mul(x[i], a[i], b[0]);
        }
    }
}

template<class T>
void sqr(Poly<T>& x, const Poly<T>& a) {
    if(IsZero(a)) { clear(x); return; }
    long i,j;
    static T t;
    i = a.length()-1; j=(i<<1);
    x.SetLength(j+1);
    sqr(x[j], a[i]);
    if(i==0) return;
    for(j=i-1; j>=0; j--) mul(x[i+j], a[i], a[j]);
    for(i--; i>0; i--) {
        sqr(t, a[i]); j=(i<<1);
        x[j] += x[j];
        x[j] += t; j++;
        x[j] += x[j];
        for(j=i-1; j>0; j--) {
            mul(t, a[i], a[j]);
            x[i+j] += t;
        }
        mul(x[i], a[i], a[0]);
    }
    x[1] += x[1];
    sqr(x[0], a[0]);
}

template<class T>
void negate(Poly<T>& x, const Poly<T>& a) {
    x.SetLength(a.length());
    for(long i=0; i<a.length(); i++) negate(x[i], a[i]);
}

template<class T>
void DivRem(Poly<T>& q, Poly<T>& r, const Poly<T>& a, const Poly<T>& b) {
    if(&r==&b) { Poly<T> t(b); DivRem(q,r,a,t); return; }
    if(&r!=&a) r=a;
    if(a.length() < b.length()) clear(q);
    else {
        long i, n(b.length()-1), m(a.length()-n);
        r %= b;
        q.SetLength(m);
        for(i=0; i<m; i++) q[i] = r[n+i];
    }
}

template<class T>
void PseudoDivRem(Poly<T>& q, Poly<T>& r, const Poly<T>& a, const Poly<T>& b) {
    if(&r==&b) { Poly<T> t(b); PseudoDivRem(q,r,a,t); return; }
    if(&r!=&a) r=a;
    if(a.length() < b.length()) { clear(q); return; }
    q.SetLength(a.length() - b.length() + 1);
    long m(!IsMonic(b)),i,j,k;
    static T t;
    for(i=q.length()-1, k=r.length()-1; i>=0; i--, k--) {
        q[i] = r[k];
        if(m) {// pseudo division
            for(j=q.length()-1; j>i; j--) q[j] *= LeadCoeff(b);
            for(j=k-1; j>=0; j--) r[j] *= LeadCoeff(b);
        }
        for(j=b.length()-2; j>=0; j--) {
            mul(t, q[i], b[j]);
            r[i+j] -= t;
        }
    }
    if(a.length() >= b.length()) {
        r.SetLength(b.length()-1);
        r.normalize();
    }
}

template<class T>
void MakeMonic(Poly<T>& a) {
    static T t;
    inv(t, LeadCoeff(a));
    for(long i=a.length()-2; i>=0; i--) a[i] *= t;
    set(LeadCoeff(a));
}

template<class T>
void diff(Poly<T>& x, const Poly<T>& a) {
    long i,m(a.length());
    x.SetLength(m-1);
    for(long i=1; i<m; i++) mul(x[i-1], a[i], i);
}

template<class T>
void translate(Poly<T>& x, const Poly<T>& a, const T& c) {
    long i,j;
    static T t;
    if(&x!=&a) x=a;
    for(i=0; i<x.length()-1; i++) {
        for(j=x.length()-2; j>=i; j--) {
            mul(t, c, x[j+1]);
            x[j] += t;
        }
    }
}

template<class T>
void RightShift(Poly<T>& x, const Poly<T>& a, long k) {// k>0
    long i, m(a.length());
    if(m<=k) clear(x);
    else {
        x.SetLength(m-k);
        for(i=k; i<m; i++) x[i-k] = a[i];
    }
}

template<class T>
void LeftShift(Poly<T>& x, const Poly<T>& a, long k) {// k>0
    long i, m(a.length());
    x.SetLength(m+k);
    for(i=m-1; i>=0; i--) x[i+k] = a[i];
    for(i=0; i<k; i++) clear(x[i]);
}

template<class T>
void scale(Poly<T>& f, const Poly<T>& g, const T& a) {// f(x)=g(ax)
    T t(a);
    if(&f!=&g) {
        f.SetLength(g.length());
        if(!IsZero(g)) f[0] = g[0];
    }
    for(long i=1; i<g.length(); i++) {
        mul(f[i], g[i], t);
        t *= a;
    }
}

template<class T>
void eval(T& x, const Poly<T>& f, const T& a) {
    if(IsZero(f)) clear(x);
    else if(&x==&a) { T t(a); eval(x,f,t); }
    else {
        long i(f.length()-1);
        x = f[i];
        for(i--; i>=0; i--) {
            x *= a;
            x += f[i];
        }
    }
}

#endif // __Poly_h__