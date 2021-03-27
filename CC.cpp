#include "CC.h"

static RR s,t;

void exp(CC& w, const CC& z) {// w=exp(z)
    exp(s, real(z));
    cos(real(w), imag(z));
    sin(imag(w), imag(z));
    w*=s;
}

void norm(RR& x, const CC& z) {// x=|z|^2
    sqr(x, real(z));
    sqr(s, imag(z));
    x += s;
}

void abs(RR& x, const CC& z) {// x=|z|
    sqr(x, real(z));
    sqr(s, imag(z));
    x += s;
    SqrRoot(x,x);
}

void sqr(CC& w, const CC& z) {// w=z^2
    sqr(s, imag(z));
    mul(t, real(z), imag(z));
    sqr(real(w), real(z));
    real(w) -= s;
    mul(imag(w), t, 2.);
}

void pow(CC& w, const CC& z, long n) {// w=z^n
    if(n==0) set(w);
    else if(n<0) {
        pow(w,z,-n);
        norm(t,w);
        negate(imag(w), imag(w));
        w/=t;
    }
    else if(&w==&z) {
        CC u(z);
        long m(1L<<(NumBits(n)-1));
        for(m>>=1; m; m>>=1) {
            sqr(w,w);
            if(n&m) w*=u;
        }
    }
    else {
        w=z;
        long m(1L<<(NumBits(n)-1));
        for(m>>=1; m; m>>=1) {
            sqr(w,w);
            if(n&m) w*=z;
        }
    }
}