// Example of Pole masses and Gf conversion to
// set of running couplings, running Higgs mass
// term and running vev in MS scheme

#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {

      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);

      // Set of all running parameters at scale Mt
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

      // Set of all running parameters at scale MZ
      P2MS pMSmZ(oi,pdg2014::Gf, as(oi.MZ()), oi.MZ(), order::all);

      
      // Input at mu=MZ for running as in [hep-ph]1208.3357
      std::cout << "alpha(1) = " << pMSmZ.a1()*4*Pi << std::endl
                << "alpha(2) = " << pMSmZ.a2()*4*Pi << std::endl
                << "alpha(3) = " << pMSmZ.as()*4*Pi << std::endl
                << "alpha(t) = " << pMSmZ.at()*4*Pi << std::endl
                << "alpha(b) = " << pMSmZ.ab()*4*Pi << std::endl
                << "4*Pi*lam = " << pMSmZ.alam()*pow(4*Pi,2) << std::endl;
        
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

