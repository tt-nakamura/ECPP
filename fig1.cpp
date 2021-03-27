#include<fstream>
#include "ECPP.h"

void fig1(const char *fname, long l) {
    std::ofstream ofs(fname);
    long i,j,k,N(12),M(20),&B(ECPP_TRYDIV_MAX);
    double p1(1<<12),p2(1<<24);
    double dp(pow(p2/p1,1./N)),t,s;
    ZZ n;
    for(i=0; i<=N; i++) {
        B = long(p1*pow(dp,i));
        t = 0;
        for(j=k=0; j<M; j++) {
            GenPrime(n,l);
            s = GetTime();
            if(AtkinMorain(n)<=0) k++;
            t += GetTime() - s;
        }
        t /= (M-k);
        s = k/double(M);
        ofs << B << ' ';
        ofs << t << ' ';
        ofs << s << std::endl;
        std::cout << B << ' ';
        std::cout << t << ' ';
        std::cout << s << std::endl;
    }
}

main() {
    fig1("fig1a.txt", 256);
    fig1("fig1b.txt", 512);
    fig1("fig1c.txt", 1024);
}