#include<iostream>

long conductor(long);
long ClassNum(long);

main()
// make table of class numbers below hmax
//  with |discriminant| <= dmax
//   for imaginary quadratic fields
{
    long dmax(10000),hmax(8),d,h;
    for(h=1; h<=hmax; h++) {
        std::cout << "Class Number ";
        std::cout << h << std::endl;
        for(d=-3; d>=-dmax; d--) {
            if((d&2)==2 || (d&2)==3) continue;
            if(conductor(d)!=1) continue;
            if(ClassNum(d) == h)
                std::cout << d << ", ";
        }
        std::cout << std::endl;
    }
}
