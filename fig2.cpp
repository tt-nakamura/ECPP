#include<fstream>
#include "ECPP.h"

void fig2(const char *fname, long l) {
    std::ofstream ofs(fname);
    long i,j,k,N(12),M(20);
    double T1(1e-4),T2(1),&T(ECPP_RHO_TIMEOUT);
    double dT(pow(T2/T1,1./N)),t,s;
    ZZ n;
    for(i=0; i<=N; i++) {
        T = T1*pow(dT,i);
        t = 0;
        for(j=k=0; j<M; j++) {
            GenPrime(n,l);
            s = GetTime();
            if(AtkinMorain(n)<=0) k++;
            t += GetTime() - s;
        }
        t /= (M-k);
        s = k/double(M);
        ofs << T << ' ';
        ofs << t << ' ';
        ofs << s << std::endl;
        std::cout << T << ' ';
        std::cout << t << ' ';
        std::cout << s << std::endl;
    }
}

main() {
    fig2("fig2a.txt", 256);
    fig2("fig2b.txt", 512);
    fig2("fig2c.txt", 1024);
}