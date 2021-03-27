#include<fstream>
#include "ECPP.h"

long apr(const ZZ&);

void fig3(const char *fname, long ptest, long l1, long l2, long N) {
    std::ofstream ofs(fname);
    long i,j,l,k,M(10),m(ptest==2 ? 1:M);
    double dl(pow(double(l2)/l1,1./N)),t,s;
    ZZ n;
    for(i=0; i<=N; i++) {
        l = long(l1*pow(dl,i));
        t = 0;
        for(j=0; j<m; j++) {
            GenPrime(n,l);
            s = GetTime();
            switch(ptest) {
                case 0: k = AtkinMorain(n); break;
                case 1: k = GoldKil(n); break;
                case 2: k = apr(n); break;
            }
            s = GetTime() - s;
            if(k<=0) j--;
            else t += s;
        }
        t /= m;
        ofs << l << ' ' << t << std::endl;
        std::cout << l << ' ' << t << std::endl;
    }
}

main() {
    fig3("fig3a.txt", 0, 60, 2000, 20);
    fig3("fig3b.txt", 1, 30, 150, 8);
    fig3("fig3c.txt", 2, 50, 1200, 10);
}