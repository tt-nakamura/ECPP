#include "ECPP.h"

main() {
    ZZ n;
    ECPPCert c;

    GenPrime(n,72);
    std::cout << "Goldwasser-Kilian primality prooving" << std::endl;
    std::cout << n << std::endl;
    std::cout << (GoldKil(c,n)>0 ? "succeeded" : "failed");
    std::cout << std::endl << "certificate is ";
    std::cout << (Certify(c) ? "verified" : "invalid");
    std::cout << std::endl << std::endl;

    GenPrime(n,1024);
    std::cout << "Atkin-Morain primality prooving" << std::endl;
    std::cout << n << std::endl;
    std::cout << (AtkinMorain(c,n)>0 ? "succeeded" : "failed");
    std::cout << std::endl <<  "certificate is ";
    std::cout << (Certify(c) ? "verified" : "invalid");
    std::cout << std::endl;
}