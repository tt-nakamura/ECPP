// uses NTL
//   http://www.shoup.net/ntl

#ifndef __CC_h__
#define __CC_h__

#include<complex>
#include<NTL/RR.h>
using namespace NTL;

// complex numbers in arbitrary-precision
typedef std::complex<RR> CC;

inline void operator+=(CC& z, double a) { real(z) += a; }
inline void operator-=(CC& z, double a) { real(z) -= a; }
inline void operator*=(CC& z, double a) { real(z) *= a, imag(z) *= a; }
inline void operator/=(CC& z, double a) { real(z) /= a, imag(z) /= a; }
inline void clear(CC& z) { clear(real(z)), clear(imag(z)); }
inline void set(CC& z) { set(real(z)), clear(imag(z)); }

void abs(RR&, const CC&);
void norm(RR&, const CC&);
void sqr(CC&, const CC&);
void pow(CC&, const CC&, long);
void exp(CC&, const CC&);

#endif // __CC_h__