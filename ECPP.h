// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ECPP_h__
#define __ECPP_h__

#include<NTL/ZZ.h>
#include "EC_p.h"
using namespace NTL;

struct CertElem {// element of certificate
    ZZ m,n;
    ZZ_p a,b;
    EC_p P;
};

typedef Vec<CertElem> ECPPCert;// certificate

void set(CertElem& e, const EC_p&, const ZZ&);

long brent_rho(ZZ&, const ZZ&, double);
long IsPrime(long);
long cornacchia(ZZ&, ZZ&, const ZZ&, const ZZ&);
long ClassPoly(ZZX&, long);
long ProbPrimeFactor(ZZ&, const ZZ&, const ZZ&);

long AtkinMorain(ECPPCert&, const ZZ&, long=0);
long AtkinMorain(const ZZ&, long=0);
long GoldKil(ECPPCert&, const ZZ&, long=0);
long GoldKil(const ZZ&, long=0);
long Certify(const ECPPCert&);

extern long ECPP_TRYDIV_MAX;
extern double ECPP_RHO_TIMEOUT;

#endif // __ECPP_h__