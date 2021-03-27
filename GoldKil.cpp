// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   R. Crandall and C. Pomerance
//    "Prime Numbers: A Computational Perspective"
//     2nd edition, Algorithm 7.6.2

#include "ECPP.h"

#define GK_NMIN (1<<30) // > ECPP_TRYDIV_MAX

long GoldKil_(ECPPCert& c, const ZZ& n, long verbose)
// primarity proving by Goldwasser-Kilian method
// return 1 or 0 if n is prime or not, respectively
// if n is prime, set c = certifcate of proof
// if verbose!=0, print recursion steps
// assume n is probable prime (by Miller-Rabin test)
{
    long i,p;
    ZZ t,q,m;
    EC_p P,Q;
    PrimeSeq ps;
    if(n <= GK_NMIN) {
        conv(p,n);
        if(!IsPrime(p)) return 0;
        c.SetLength(1);
        c[0].n = n;
        return 1;
    }
    SqrRoot(t,n); LeftShift(q,t,2);
    SqrRoot(q,q);
    t+=q; t++;// t = ([n^{1/4}] + 1)^2
    ZZ_p::init(n);
    for(;;) {
        EC_p::init_random();
        if(EC_p::schoof(m)) return 0;
        if(verbose)
            std::cout << "curve order : " << m << std::endl;
        if(ProbPrimeFactor(q,m,t)) continue;
        div(t,m,q);
        do {
            random(P);
            mul(Q,P,m);
            if(!IsZero(Q)) return 0;
            mul(Q,P,t);
        } while(IsZero(Q));
        if(!z_coprime(Q)) return 0;
        if(verbose)
            std::cout << "prime factor: " << q << std::endl;
        EC_pPush push;
        if(GoldKil(c,q,verbose)) break;
    }
    i = c.length();
    c.SetLength(i+1);
    set(c[i], P, m);
    return 1;
}

long GoldKil(ECPPCert& c, const ZZ& n, long verbose)
// perform Miller-Rabin test before proving primality
{
    if(!ProbPrime(n)) return 0;
    c.SetLength(0);
    return GoldKil_(c,n,verbose);
}

long GoldKil(const ZZ& n, long verbose)
// do not output certificate
{
    ECPPCert c;
    return GoldKil(c,n,verbose);
}