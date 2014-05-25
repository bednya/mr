#include <iostream>
#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Compare with:
      
      SMinput KVPhys(4.40, 80.385, 91.1876, 125.6, 173.5);
      

      alphaS   = 0.1184;

      // \mu = Mt
      alphaMt  = 0.00779305;
      alphaSMt   = 0.1079;

      // \mu = Mb
      long double alphaMb  = 0.00784257;
      long double alphaSMb   = 0.1905;

      

      AlphaS as;
      
      long double alphaMZ = 1./137.035999;


      // Yukawa top
      tt dMt  = tt(KVPhys, KVPhys.MMt());
      std::cout << "[ Top quark ]" << std::endl;
      std::cout << "Mh= " << KVPhys.MH()  << std::endl;
      std::cout << "as(MMt) = " << as(KVPhys.MMt()) << std::endl;          
      std::cout << "\t1-loop \\alpha         " << alphaMt/4./Pi*dMt.my10() << std::endl;
      std::cout << "\t1-loop \\alpha_S       " << alphaSMt/4./Pi*dMt.my01() << std::endl;
      std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMt/4./Pi*alphaSMt/4./Pi*dMt.my11() << std::endl;
      std::cout << "\t2-loop \\alpha^2       " << pow(alphaMt/4./Pi,2)*dMt.my20() << std::endl;


      // Yukawa bottom
      bb dMb  = bb(KVPhys, KVPhys.MMb());
      std::cout << "[ Bottom quark ]" << std::endl;
      std::cout << "\t1-loop \\alpha         " << alphaMb/4./Pi*dMb.my10() << std::endl;
      std::cout << "\t1-loop \\alpha_S       " << alphaSMb/4./Pi*dMb.my01() << std::endl;
      std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMb/4./Pi*alphaSMb/4./Pi*dMb.my11() << std::endl;
      std::cout << "\t2-loop \\alpha^2       " << pow(alphaMb/4./Pi,2)*dMb.my20() << std::endl;
      
      
      std::cout << "[ ratios ]" << std::endl;
      std::cout << "<10> yt/yb = " << dMb.my10()/dMt.my10() << " mt/mb = " << dMb.m10()/dMt.m10() << std::endl;
      std::cout << "<01> yt/yb = " << dMb.my01()/dMt.my01() << " mt/mb = " << dMb.m01()/dMt.m01() << std::endl;
      std::cout << "<11> yt/yb = " << dMb.my11()/dMt.my11() << " mt/mb = " << dMb.m11()/dMt.m11() << std::endl;
      std::cout << "<20> yt/yb = " << dMb.my20()/dMt.my20() << " mt/mb = " << dMb.m20()/dMt.m20() << std::endl;
      
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}
