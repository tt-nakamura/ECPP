// uses NTL
//   http://www.shoup.net/ntl

#ifndef __EC_p_h__
#define __EC_p_h__

#include<NTL/ZZ_pX.h>
#include<NTL/ZZX.h>
using namespace NTL;

struct EC_p {// Elliptic Curve over finite field F_p
    ZZ_p x,z;// montgomery coordinates
    static ZZ d;// d=gcd(z,p) to check d==1
    static ZZ_pX f;// cubic polynomial x^3+ax+b
    inline static ZZ_p& a() { return f[1]; }
    inline static ZZ_p& b() { return f[0]; }
    static long IsSingular();
    static void init_random();
    static void init(const ZZ&, const ZZ&);
    static void init(const ZZ_p&, const ZZ_p&);
    static long schoof(ZZ&);// order counting
};

inline long IsZero(const EC_p& a) { return IsZero(a.z); }
inline void clear(EC_p& a) { clear(a.z); }
inline void set(EC_p& a, const ZZ_p& x1) { a.x=x1; set(a.z); }

void addh(EC_p&, const EC_p&, const EC_p&, const EC_p&);
void doubleh(EC_p&, const EC_p&);
void mul(EC_p&, const EC_p&, const ZZ&);
void mul(EC_p&, const EC_p&, long);
void mul(EC_p&, EC_p&, const EC_p&, const ZZ&);
long x_coord(ZZ&, const EC_p&);
long z_coprime(const EC_p&);
long IsOnEC(const EC_p&);
void random(EC_p& a);
std::ostream& operator<<(std::ostream&, const EC_p&);

long CountOrder(long, long, long);
long Jacobi(long, long);

struct EC_pPush {// save context
    ZZ_pX f;// cubic polynomial x^3+ax+b
    ZZ_pPush push;
    EC_pPush() : f(EC_p::f) {;}
    ~EC_pPush() { EC_p::f.kill(); EC_p::f = f; }
    // kill f to avoid "internal error: can't grow this _ntl_gbigint"
};


#endif // __EC_p_h__