// uses NTL
//   http://www.shoup.net/ntl

#ifndef __EC_pE_h__
#define __EC_pE_h__

#include<NTL/ZZ_pE.h>
#include "EC_p.h"
using namespace NTL;

struct EC_pE {// torsion points on EC over algebraic closure of F_p
	ZZ_pE x,y,z;// projective coordinates
	static Vec<ZZ_pE> ff;// ff[0]=x^3+ax+b, ff[1]=ff[0]^2
	static Vec<ZZ_pX> g;// division polynomial
	static void init(long);
	static void DivPoly(long);
	inline static const ZZ_pE& f() { return ff[0]; }
	inline static const ZZ_pE& f2() { return ff[1]; }
	inline static void kill() { ff.kill(); g.kill(); }
};

inline void clear(EC_pE& a) { clear(a.z); }
inline long IsZero(const EC_pE& a) { return IsZero(a.z); }
void set(EC_pE&);
void set(EC_pE&, long);
void negate(EC_pE&, const EC_pE&);
void add(EC_pE&, const EC_pE&, const EC_pE&);
void sub(EC_pE&, const EC_pE&, const EC_pE&);
void doubleh(EC_pE&, const EC_pE&);
void mul(EC_pE&, const EC_pE&, long);
long compare_x(const EC_pE&, const EC_pE&);
long compare_y(const EC_pE&, const EC_pE&);
void FrobConj(EC_pE&, const EC_pE&);

#endif // __EC_pE_h__