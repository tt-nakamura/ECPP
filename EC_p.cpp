// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   R. Crandall and C. Pomerance
//    "Prime Numbers: A Computational Perspective"
//     2nd edition, section 7.2

#include "EC_p.h"

ZZ EC_p::d;
ZZ_pX EC_p::f;

void EC_p::init(const ZZ_p& a, const ZZ_p& b)
// y^2 = x^3 + ax + b
// assume ZZ_p has been initialized
// kill f to avoid "internal error: can't grow this _ntl_gbigint"
{
    f.kill();
	SetCoeff(f,3);
	SetCoeff(f,1,a);
	SetCoeff(f,0,b);
}

void EC_p::init(const ZZ& a, const ZZ& b)
// y^2 = x^3 + ax + b
{
    f.kill();
	SetCoeff(f,3);
	conv(f[1], a);
	conv(f[0], b);
}

void EC_p::init_random()
// random choice of (a,b)
{
    f.kill();
	SetCoeff(f,3);
	do {
        random(f[1]);
		random(f[0]);
	} while(IsSingular());
}

long EC_p::IsSingular()
// return (discriminant 4a^3+27b^2 == 0)
{
	ZZ_p u,v;
	power(u,a(),3);
	sqr(v,b());
	u *= 4;
	v *= 27;
	u += v;
    return IsZero(u);
}

std::ostream& operator<<(std::ostream& s, const EC_p& a) {
	s << '[' << a.x << ' ' << a.z << ']';
	return s;
}

void addh(EC_p& p, const EC_p& a, const EC_p& b, const EC_p& q) 
// set p=a+b; require q=|a-b|; assume &p!=&q
// reference: Crandall & Pomerance, eq(7.6) with C=0
{
    if(IsZero(a)) p=b;
	else if(IsZero(b)) p=a;
	else if(IsZero(q)) doubleh(p,a);
	else if(IsZero(q.x)) {
		ZZ_p s,t,u;
		mul(s, a.x, b.z);
		mul(t, a.z, b.x);
		mul(u, a.z, b.z);
		sub(p.z, s, t);
		sqr(p.z, p.z);
		mul(p.x, a.x, b.x);
		s += t;
		mul(t, u, EC_p::a());
		p.x += t;
		p.x *= s;
		p.x += p.x;
		sqr(u,u);
		u *= EC_p::b();
		u *= 4;
		p.x += u;		
	}
	else {
        ZZ_p s,t,u,v;
		mul(s, a.x, b.z);
		mul(t, a.z, b.x);
		sub(u, s, t);
		sqr(u, u);
		mul(v, a.z, b.z);
		mul(p.z, u, q.x);
		s += t;
		s *= v;
		s *= EC_p::b();
		s *= 4;
		mul(u, a.x, b.x);
		v *= EC_p::a();
		u -= v;
		sqr(u,u);
		u -= s;
		mul(p.x, u, q.z);
	}
}

void doubleh(EC_p& p, const EC_p& a)
// set p=a+a
// reference: Crandall & Pomerance, eq(7.7) with C=0
{
    if(IsZero(a)) clear(p);
	else {
        ZZ_p x2,az2,bz3,u;
		sqr(x2, a.x);
		sqr(az2, a.z);
		mul(bz3, az2, a.z);
		az2 *= EC_p::a();
		bz3 *= EC_p::b();
		add(u, x2, az2);
		u *= a.x;
		u += bz3;
		u *= 4;
		mul(p.z, a.z, u);
		x2 -= az2;
		sqr(x2,x2);
		mul(u, a.x, bz3);
		u *= 8;
		sub(p.x, x2, u);
	}	
}

void mul(EC_p& p, const EC_p& a, const ZZ& n)
// set p=n*a; assume n>0
// reference: Crandall & Pomerance, Algorithm 7.2.7
{
    if(IsZero(a) || IsZero(n)) clear(p);
	else if(IsOne(n)) p=a;
	else if(n==2) doubleh(p,a);
	else if(&p==&a) {
		EC_p b(a);
		mul(p,b,n);
	}
	else {
        EC_p q;
		p = a;
		doubleh(q,a);
		for(long i=NumBits(n)-2; i>=1; i--) {
			if(bit(n,i)) {
				addh(p,q,p,a);
				doubleh(q,q);
			}
			else {
                addh(q,q,p,a);
				doubleh(p,p);
			}
		}
		if(bit(n,0)) addh(p,q,p,a);
		else doubleh(p,p);
	}
}

void mul(EC_p& p, EC_p& q, const EC_p& a, const ZZ& n)
// set p=n*a, set q=(n+1)*a; assume n>0
{
    if(IsZero(a)) { clear(p); clear(q); }
	else if(IsZero(n)) { clear(p); q=a; } 
	else if(IsOne(n)) { p=a; doubleh(q,a); }
	else if(n==2) { doubleh(p,a); addh(q,p,a,a); }
	else if(&p==&a || &q==&a) {
		EC_p b(a);
		mul(p,q,b,n);
	}
	else {
        p = a;
		doubleh(q,a);
		for(long i=NumBits(n)-2; i>=0; i--) {
			if(bit(n,i)) {
				addh(p,q,p,a);
				doubleh(q,q);
			}
			else {
                addh(q,q,p,a);
				doubleh(p,p);
			}
		}
	}
}

void mul(EC_p& p, const EC_p& a, long n)
// set p=n*a; assume n>0
{
    if(IsZero(a) || n==0) clear(p);
	else if(n==1) p=a;
	else if(n==2) doubleh(p,a);
	else if(&p==&a) {
		EC_p b(a);
		mul(p,b,n);
	}
	else {
        unsigned long k(1L<<(NumBits(n)-1));
		EC_p q;
		p = a;
		doubleh(q,a);
		for(k>>=1; k>1; k>>=1) {// Lucas chain
			if(n&k) {
				addh(p,q,p,a);
				doubleh(q,q);
			}
			else {
                addh(q,q,p,a);
				doubleh(p,p);
			}
		}
		if(n&1) addh(p,q,p,a);
		else doubleh(p,p);
	}
}

long x_coord(ZZ_p& x, const EC_p& a)
// if a.z is invertible mod p, set x=a.x/a.z and return 0
// else set d=(divisor of p) and return 1
{
	ZZ s,t;
	XGCD(EC_p::d, s, t, rep(a.z), ZZ_p::modulus());
	if(!IsOne(EC_p::d)) return 1;
	if(sign(s) < 0) s += ZZ_p::modulus();
	conv(x,s);
	x *= a.x;
	return 0;
}

long z_coprime(const EC_p& a)
// if a.z is invertible mod p, return 1
// else set d=(divisor of p) and return 0
{
	GCD(EC_p::d, rep(a.z), ZZ_p::modulus());
	return IsOne(EC_p::d);
}

long IsOnEC(const EC_p& a)
// return 1 or 0 if a is on the curve or not
{
	if(IsZero(a)) return 1;
	ZZ_p x,y;
	if(x_coord(x,a)) return 0;
	eval(y, EC_p::f, x);
	return Jacobi(rep(y), ZZ_p::modulus())>=0;
}

void random(EC_p& p)
// set p = random point on the curve (except infty)
{
	do {
        random(p.x);
		eval(p.z, EC_p::f, p.x);
	} while(Jacobi(rep(p.z), ZZ_p::modulus()) < 0);
	set(p.z);
}

long CountOrder(long a, long b, long p)
// count the number of points on the curve
// reference: Crandall & Pomerance, eq(7.8)
{
	long x,s,t;
	for(x=s=0; x<p; x++) {
		t = MulMod(x,x,p);
		t = AddMod(t,a,p);
		t = MulMod(t,x,p);
		t = AddMod(t,b,p);
		s += Jacobi(t,p);
	}
	return p+1+s;
}