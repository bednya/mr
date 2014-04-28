#include <iostream>
#include "betaSM.hpp"
#include "SM3r2m.hpp"

BetaSMFull::BetaSMFull(int pocoa1_, int pocoa2_, int pocoas_, int pocoat_, int pocoab_, int pocoatau_, int pocolam_)
{
 // pocoa1(pocoa1_) pocoa2(pocoa2_) pocoas(pocoas_) pocoat(pocoat_) pocoab(pocoab_) pocoatau(pocoatau_) pocolam(pocolam_);
  // size_t powers[] = {pocoa1,pocoa2,pocoas,pocoat,pocoab,pocoatau,pocolam};
  // maxPower = *std::max_element(powers, powers + 7);
  
  for(size_t ia1 = 0; ia1 < maxPoco; ia1++)
    for(size_t ia2 = 0; ia2 < maxPoco - ia1; ia2++)
      for(size_t ias = 0; ias < maxPoco - ia1 - ia2; ias++)
        for(size_t iat = 0; iat < maxPoco - ia1 - ia2 - ias; iat++)
          for(size_t iab = 0; iab < maxPoco - ia1 - ia2 - ias - iat; iab++)
            for(size_t iatau = 0; iatau < maxPoco - ia1 - ia2 - ias - iat - iab; iatau++)
              for(size_t ilam = 0; ilam < maxPoco - ia1 - ia2 - ias - iat - iab - iatau; ilam++)
                {
                  be1[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be2[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be3[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be4[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be5[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be6[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                  be7[ia1][ia2][ias][iat][iab][iatau][ilam] = 0;
                }
    
  // 
  // Gauge beta-functions
  // 
    
  // * a1
  be1[2][0][0][0][0][0][0] = 4.1;
  be1[2][0][0][0][0][0][2] = -1.8;
  be1[2][0][0][0][0][1][0] = -1.5;
  be1[2][0][0][0][0][2][0] = 5.7375;
  be1[2][0][0][0][1][0][0] = -0.5;
  be1[2][0][0][0][1][1][0] = 7.85;
  be1[2][0][0][0][2][0][0] = 3.5625;
  be1[2][0][0][1][0][0][0] = -1.7;
  be1[2][0][0][1][0][1][0] = 9.95;
  be1[2][0][0][1][1][0][0] = 9.225;
  be1[2][0][0][2][0][0][0] = 11.8125;
  be1[2][0][1][0][0][0][0] = 8.8;
  be1[2][0][1][0][1][0][0] = -3.4;
  be1[2][0][1][1][0][0][0] = -5.8;
  be1[2][0][2][0][0][0][0] = 59.4;
  be1[2][1][0][0][0][0][0] = 2.7;
  be1[2][1][0][0][0][0][1] = 0.9;
  be1[2][1][0][0][0][1][0] = -10.18125;
  be1[2][1][0][0][1][0][0] = -8.19375;
  be1[2][1][0][1][0][0][0] = -14.71875;
  be1[2][1][1][0][0][0][0] = -0.6;
  be1[2][2][0][0][0][0][0] = 12.328125;
  be1[3][0][0][0][0][0][0] = 3.98;
  be1[3][0][0][0][0][0][1] = 0.54;
  be1[3][0][0][0][0][1][0] = -3.16125;
  be1[3][0][0][0][1][0][0] = -1.58375;
  be1[3][0][0][1][0][0][0] = -3.53375;
  be1[3][0][1][0][0][0][0] = -1.8266666666666666666666666666666666666666666666666666666666666666667;
  be1[3][1][0][0][0][0][0] = 0.76875;
  be1[4][0][0][0][0][0][0] = -16.192208333333333333333333333333333333333333333333333333333333333333;

  // * a2
  be2[0][2][0][0][0][0][0] = -3.1666666666666666666666666666666666666666666666666666666666666666667;
  be2[0][2][0][0][0][0][2] = -3.;
  be2[0][2][0][0][0][1][0] = -0.5;
  be2[0][2][0][0][0][2][0] = 1.8125;
  be2[0][2][0][0][1][0][0] = -1.5;
  be2[0][2][0][0][1][1][0] = 3.75;
  be2[0][2][0][0][2][0][0] = 9.1875;
  be2[0][2][0][1][0][0][0] = -1.5;
  be2[0][2][0][1][0][1][0] = 3.75;
  be2[0][2][0][1][1][0][0] = 14.625;
  be2[0][2][0][2][0][0][0] = 9.1875;
  be2[0][2][1][0][0][0][0] = 12.;
  be2[0][2][1][0][1][0][0] = -7.;
  be2[0][2][1][1][0][0][0] = -7.;
  be2[0][2][2][0][0][0][0] = 81.;
  be2[0][3][0][0][0][0][0] = 5.8333333333333333333333333333333333333333333333333333333333333333333;
  be2[0][3][0][0][0][0][1] = 1.5;
  be2[0][3][0][0][0][1][0] = -7.59375;
  be2[0][3][0][0][1][0][0] = -22.78125;
  be2[0][3][0][1][0][0][0] = -22.78125;
  be2[0][3][1][0][0][0][0] = 39.;
  be2[0][4][0][0][0][0][0] = 188.05150462962962962962962962962962962962962962962962962962962962963;
  be2[1][2][0][0][0][0][0] = 0.9;
  be2[1][2][0][0][0][0][1] = 0.3;
  be2[1][2][0][0][0][1][0] = -1.59375;
  be2[1][2][0][0][1][0][0] = -3.33125;
  be2[1][2][0][1][0][0][0] = -3.70625;
  be2[1][2][1][0][0][0][0] = -0.2;
  be2[1][3][0][0][0][0][0] = 5.45625;
  be2[2][2][0][0][0][0][0] = -3.498125;

  // * as
  be3[0][0][2][0][0][0][0] = -7.;
  be3[0][0][2][0][1][0][0] = -2.;
  be3[0][0][2][0][1][1][0] = 3.5;
  be3[0][0][2][0][2][0][0] = 15.;
  be3[0][0][2][1][0][0][0] = -2.;
  be3[0][0][2][1][0][1][0] = 3.5;
  be3[0][0][2][1][1][0][0] = 18.;
  be3[0][0][2][2][0][0][0] = 15.;
  be3[0][0][3][0][0][0][0] = -26.;
  be3[0][0][3][0][1][0][0] = -40.;
  be3[0][0][3][1][0][0][0] = -40.;
  be3[0][0][4][0][0][0][0] = 32.5;
  be3[0][1][2][0][0][0][0] = 4.5;
  be3[0][1][2][0][1][0][0] = -11.625;
  be3[0][1][2][1][0][0][0] = -11.625;
  be3[0][1][3][0][0][0][0] = 21.;
  be3[0][2][2][0][0][0][0] = 13.625;
  be3[1][0][2][0][0][0][0] = 1.1;
  be3[1][0][2][0][1][0][0] = -2.225;
  be3[1][0][2][1][0][0][0] = -2.525;
  be3[1][0][3][0][0][0][0] = 5.1333333333333333333333333333333333333333333333333333333333333333333;
  be3[1][1][2][0][0][0][0] = -0.075;
  be3[2][0][2][0][0][0][0] = -4.3583333333333333333333333333333333333333333333333333333333333333333;

  // 
  // Yukawa-couplings
  // 

  // * at
  be4[0][0][0][1][0][0][2] = 6.;
  be4[0][0][0][1][0][0][3] = -36.;
  be4[0][0][0][1][0][1][0] = 1.;
  be4[0][0][0][1][0][1][2] = -22.5;
  be4[0][0][0][1][0][2][0] = -2.25;
  be4[0][0][0][1][0][2][1] = 15.;
  be4[0][0][0][1][0][3][0] = 4.4375 + 3.*Zeta3;
  be4[0][0][0][1][1][0][0] = 1.5;
  be4[0][0][0][1][1][0][2] = -72.75;
  be4[0][0][0][1][1][1][0] = 1.25;
  be4[0][0][0][1][1][2][0] = 13.25;
  be4[0][0][0][1][2][0][0] = -0.25;
  be4[0][0][0][1][2][0][1] = 15.;
  be4[0][0][0][1][2][1][0] = 22.;
  be4[0][0][0][1][3][0][0] = 29.8125 + 4.5*Zeta3;
  be4[0][0][0][2][0][0][0] = 4.5;
  be4[0][0][0][2][0][0][1] = -12.;
  be4[0][0][0][2][0][0][2] = 3.75;
  be4[0][0][0][2][0][1][0] = -2.25;
  be4[0][0][0][2][0][1][1] = 30.;
  be4[0][0][0][2][0][2][0] = 25.875;
  be4[0][0][0][2][1][0][0] = -2.75;
  be4[0][0][0][2][1][0][1] = 93.;
  be4[0][0][0][2][1][1][0] = 3.5;
  be4[0][0][0][2][2][0][0] = 103.125 - 48.*Zeta3;
  be4[0][0][0][3][0][0][0] = -12.;
  be4[0][0][0][3][0][0][1] = 198.;
  be4[0][0][0][3][0][1][0] = 10.5;
  be4[0][0][0][3][1][0][0] = 46.1875;
  be4[0][0][0][4][0][0][0] = 42.375 + 13.5*Zeta3;
  be4[0][0][1][1][0][0][0] = -8.;
  be4[0][0][1][1][1][0][0] = 4.;
  be4[0][0][1][1][1][1][0] = -7.1666666666666666666666666666666666666666666666666666666666666666667;
  be4[0][0][1][1][2][0][0] = 82. - 64.*Zeta3;
  be4[0][0][1][2][0][0][0] = 36.;
  be4[0][0][1][2][0][0][1] = 16.;
  be4[0][0][1][2][0][1][0] = 2.5;
  be4[0][0][1][2][1][0][0] = 27. - 32.*Zeta3;
  be4[0][0][1][3][0][0][0] = -157.;
  be4[0][0][2][1][0][0][0] = -108.;
  be4[0][0][2][1][1][0][0] = -152.5 - 44.*Zeta3;
  be4[0][0][2][2][0][0][0] = 637.83333333333333333333333333333333333333333333333333333333333333333 - 228.*Zeta3;
  be4[0][0][3][1][0][0][0] = -1388.6666666666666666666666666666666666666666666666666666666666666667 + 640.*Zeta3;
  be4[0][1][0][1][0][0][0] = -2.25;
  be4[0][1][0][1][0][0][2] = 45.;
  be4[0][1][0][1][0][1][0] = 1.875;
  be4[0][1][0][1][0][2][0] = -19.6875 + 9.*Zeta3;
  be4[0][1][0][1][1][0][0] = 6.1875;
  be4[0][1][0][1][1][1][0] = -19.125 + 9.*Zeta3;
  be4[0][1][0][1][2][0][0] = -71.34375 + 31.5*Zeta3;
  be4[0][1][0][2][0][0][0] = 14.0625;
  be4[0][1][0][2][0][0][1] = -67.5;
  be4[0][1][0][2][0][1][0] = -20.25 - 9.*Zeta3;
  be4[0][1][0][2][1][0][0] = -72.09375 - 4.5*Zeta3;
  be4[0][1][0][3][0][0][0] = -99.5625;
  be4[0][1][1][1][0][0][0] = 9.;
  be4[0][1][1][1][1][0][0] = -13.5 - 108.*Zeta3;
  be4[0][1][1][2][0][0][0] = -168. + 180.*Zeta3;
  be4[0][1][2][1][0][0][0] = 246.75 - 144.*Zeta3;
  be4[0][2][0][1][0][0][0] = -5.75;
  be4[0][2][0][1][0][0][1] = -10.6875;
  be4[0][2][0][1][0][1][0] = 8.6953125 - 20.25*Zeta3;
  be4[0][2][0][1][1][0][0] = 40.39453125 - 28.125*Zeta3;
  be4[0][2][0][2][0][0][0] = 126.52734375 - 91.125*Zeta3;
  be4[0][2][1][1][0][0][0] = 108.75 - 108.*Zeta3;
  be4[0][3][0][1][0][0][0] = 0.78993055555555555555555555555555555555555555555555555555555555555556 + 140.625*Zeta3;
  be4[1][0][0][1][0][0][0] = -0.85;
  be4[1][0][0][1][0][0][2] = 9.;
  be4[1][0][0][1][0][1][0] = 1.875;
  be4[1][0][0][1][0][2][0] = -1.6875 - 5.4*Zeta3;
  be4[1][0][0][1][1][0][0] = 0.0875;
  be4[1][0][0][1][1][1][0] = 4.0916666666666666666666666666666666666666666666666666666666666666667 - 5.4*Zeta3;
  be4[1][0][0][1][2][0][0] = -5.99375 + 1.9*Zeta3;
  be4[1][0][0][2][0][0][0] = 4.9125;
  be4[1][0][0][2][0][0][1] = -12.7;
  be4[1][0][0][2][0][1][0] = -12.6 + 4.8*Zeta3;
  be4[1][0][0][2][1][0][0] = -8.64375 + 0.5*Zeta3;
  be4[1][0][0][3][0][0][0] = -30.4625;
  be4[1][0][1][1][0][0][0] = 1.2666666666666666666666666666666666666666666666666666666666666666667;
  be4[1][0][1][1][1][0][0] = -15.233333333333333333333333333333333333333333333333333333333333333333 - 5.6*Zeta3;
  be4[1][0][1][2][0][0][0] = -25.2 + 36.*Zeta3;
  be4[1][0][2][1][0][0][0] = 27.216666666666666666666666666666666666666666666666666666666666666667 - 35.2*Zeta3;
  be4[1][1][0][1][0][0][0] = -0.45;
  be4[1][1][0][1][0][0][1] = 2.925;
  be4[1][1][0][1][0][1][0] = -3.253125 - 1.8*Zeta3;
  be4[1][1][0][1][1][0][0] = 5.8359375 + 2.7*Zeta3;
  be4[1][1][0][2][0][0][0] = 12.6515625 + 18.45*Zeta3;
  be4[1][1][1][1][0][0][0] = -16.05;
  be4[1][2][0][1][0][0][0] = 2.559375 - 6.075*Zeta3;
  be4[2][0][0][1][0][0][0] = 1.9783333333333333333333333333333333333333333333333333333333333333333;
  be4[2][0][0][1][0][0][1] = -2.7225;
  be4[2][0][0][1][0][1][0] = -15.0921875 - 8.07*Zeta3;
  be4[2][0][0][1][1][0][0] = -2.1183854166666666666666666666666666666666666666666666666666666666667 - 0.995*Zeta3;
  be4[2][0][0][2][0][0][0] = -23.863489583333333333333333333333333333333333333333333333333333333333 - 0.465*Zeta3;
  be4[2][0][1][1][0][0][0] = 13.646666666666666666666666666666666666666666666666666666666666666667 - 29.92*Zeta3;
  be4[2][1][0][1][0][0][0] = 3.834375 - 6.885*Zeta3;
  be4[3][0][0][1][0][0][0] = 31.813458333333333333333333333333333333333333333333333333333333333333 - 13.073*Zeta3;

  // * ab
  be5[0][0][0][0][1][0][2] = 6.;
  be5[0][0][0][0][1][0][3] = -36.;
  be5[0][0][0][0][1][1][0] = 1.;
  be5[0][0][0][0][1][1][2] = -22.5;
  be5[0][0][0][0][1][2][0] = -2.25;
  be5[0][0][0][0][1][2][1] = 15.;
  be5[0][0][0][0][1][3][0] = 4.4375 + 3.*Zeta3;
  be5[0][0][0][0][2][0][0] = 4.5;
  be5[0][0][0][0][2][0][1] = -12.;
  be5[0][0][0][0][2][0][2] = 3.75;
  be5[0][0][0][0][2][1][0] = -2.25;
  be5[0][0][0][0][2][1][1] = 30.;
  be5[0][0][0][0][2][2][0] = 25.875;
  be5[0][0][0][0][3][0][0] = -12.;
  be5[0][0][0][0][3][0][1] = 198.;
  be5[0][0][0][0][3][1][0] = 10.5;
  be5[0][0][0][0][4][0][0] = 42.375 + 13.5*Zeta3;
  be5[0][0][0][1][1][0][0] = 1.5;
  be5[0][0][0][1][1][0][2] = -72.75;
  be5[0][0][0][1][1][1][0] = 1.25;
  be5[0][0][0][1][1][2][0] = 13.25;
  be5[0][0][0][1][2][0][0] = -2.75;
  be5[0][0][0][1][2][0][1] = 93.;
  be5[0][0][0][1][2][1][0] = 3.5;
  be5[0][0][0][1][3][0][0] = 46.1875;
  be5[0][0][0][2][1][0][0] = -0.25;
  be5[0][0][0][2][1][0][1] = 15.;
  be5[0][0][0][2][1][1][0] = 22.;
  be5[0][0][0][2][2][0][0] = 103.125 - 48.*Zeta3;
  be5[0][0][0][3][1][0][0] = 29.8125 + 4.5*Zeta3;
  be5[0][0][1][0][1][0][0] = -8.;
  be5[0][0][1][0][2][0][0] = 36.;
  be5[0][0][1][0][2][0][1] = 16.;
  be5[0][0][1][0][2][1][0] = 2.5;
  be5[0][0][1][0][3][0][0] = -157.;
  be5[0][0][1][1][1][0][0] = 4.;
  be5[0][0][1][1][1][1][0] = -7.1666666666666666666666666666666666666666666666666666666666666666667;
  be5[0][0][1][1][2][0][0] = 27. - 32.*Zeta3;
  be5[0][0][1][2][1][0][0] = 82. - 64.*Zeta3;
  be5[0][0][2][0][1][0][0] = -108.;
  be5[0][0][2][0][2][0][0] = 637.83333333333333333333333333333333333333333333333333333333333333333 - 228.*Zeta3;
  be5[0][0][2][1][1][0][0] = -152.5 - 44.*Zeta3;
  be5[0][0][3][0][1][0][0] = -1388.6666666666666666666666666666666666666666666666666666666666666667 + 640.*Zeta3;
  be5[0][1][0][0][1][0][0] = -2.25;
  be5[0][1][0][0][1][0][2] = 45.;
  be5[0][1][0][0][1][1][0] = 1.875;
  be5[0][1][0][0][1][2][0] = -19.6875 + 9.*Zeta3;
  be5[0][1][0][0][2][0][0] = 14.0625;
  be5[0][1][0][0][2][0][1] = -67.5;
  be5[0][1][0][0][2][1][0] = -20.25 - 9.*Zeta3;
  be5[0][1][0][0][3][0][0] = -99.5625;
  be5[0][1][0][1][1][0][0] = 6.1875;
  be5[0][1][0][1][1][1][0] = -19.125 + 9.*Zeta3;
  be5[0][1][0][1][2][0][0] = -72.09375 - 4.5*Zeta3;
  be5[0][1][0][2][1][0][0] = -71.34375 + 31.5*Zeta3;
  be5[0][1][1][0][1][0][0] = 9.;
  be5[0][1][1][0][2][0][0] = -168. + 180.*Zeta3;
  be5[0][1][1][1][1][0][0] = -13.5 - 108.*Zeta3;
  be5[0][1][2][0][1][0][0] = 246.75 - 144.*Zeta3;
  be5[0][2][0][0][1][0][0] = -5.75;
  be5[0][2][0][0][1][0][1] = -10.6875;
  be5[0][2][0][0][1][1][0] = 23.6953125 - 24.75*Zeta3;
  be5[0][2][0][0][2][0][0] = 126.52734375 - 91.125*Zeta3;
  be5[0][2][0][1][1][0][0] = 40.39453125 - 28.125*Zeta3;
  be5[0][2][1][0][1][0][0] = 108.75 - 108.*Zeta3;
  be5[0][3][0][0][1][0][0] = 0.78993055555555555555555555555555555555555555555555555555555555555556 + 140.625*Zeta3;
  be5[1][0][0][0][1][0][0] = -0.25;
  be5[1][0][0][0][1][0][2] = 9.;
  be5[1][0][0][0][1][1][0] = 1.875;
  be5[1][0][0][0][1][2][0] = -1.6875 - 5.4*Zeta3;
  be5[1][0][0][0][2][0][0] = 2.9625;
  be5[1][0][0][0][2][0][1] = -13.9;
  be5[1][0][0][0][2][1][0] = -13.7 + 4.2*Zeta3;
  be5[1][0][0][0][3][0][0] = -24.7625;
  be5[1][0][0][1][1][0][0] = 1.1375;
  be5[1][0][0][1][1][1][0] = 2.4416666666666666666666666666666666666666666666666666666666666666667 - 3.6*Zeta3;
  be5[1][0][0][1][2][0][0] = -26.26875 + 7.7*Zeta3;
  be5[1][0][0][2][1][0][0] = -2.26875 - 1.7*Zeta3;
  be5[1][0][1][0][1][0][0] = 2.0666666666666666666666666666666666666666666666666666666666666666667;
  be5[1][0][1][0][2][0][0] = -18. + 26.4*Zeta3;
  be5[1][0][1][1][1][0][0] = -26.833333333333333333333333333333333333333333333333333333333333333333 + 4.*Zeta3;
  be5[1][0][2][0][1][0][0] = 69.416666666666666666666666666666666666666666666666666666666666666667 - 35.2*Zeta3;
  be5[1][1][0][0][1][0][0] = -1.35;
  be5[1][1][0][0][1][0][1] = -0.675;
  be5[1][1][0][0][1][1][0] = -3.853125 + 10.8*Zeta3;
  be5[1][1][0][0][2][0][0] = 13.2703125 - 1.8*Zeta3;
  be5[1][1][0][1][1][0][0] = 5.1046875 + 9.45*Zeta3;
  be5[1][1][1][0][1][0][0] = -7.65;
  be5[1][2][0][0][1][0][0] = -1.978125 - 6.075*Zeta3;
  be5[2][0][0][0][1][0][0] = -0.21166666666666666666666666666666666666666666666666666666666666666667;
  be5[2][0][0][0][1][0][1] = -0.5625;
  be5[2][0][0][0][1][1][0] = -14.1371875 + 3.51*Zeta3;
  be5[2][0][0][0][2][0][0] = -10.919739583333333333333333333333333333333333333333333333333333333333 - 0.855*Zeta3;
  be5[2][0][0][1][1][0][0] = -5.4546354166666666666666666666666666666666666666666666666666666666667 - 0.065*Zeta3;
  be5[2][0][1][0][1][0][0] = -4.4933333333333333333333333333333333333333333333333333333333333333333 - 8.8*Zeta3;
  be5[2][1][0][0][1][0][0] = 6.684375 - 2.025*Zeta3;
  be5[3][0][0][0][1][0][0] = 11.655125 - 3.845*Zeta3;

  // * atau
  be6[0][0][0][0][0][1][2] = 6.;
  be6[0][0][0][0][0][1][3] = -36.;
  be6[0][0][0][0][0][2][0] = 2.5;
  be6[0][0][0][0][0][2][1] = -12.;
  be6[0][0][0][0][0][2][2] = 48.75;
  be6[0][0][0][0][0][3][0] = -3.;
  be6[0][0][0][0][0][3][1] = 108.;
  be6[0][0][0][0][0][4][0] = -10. + 7.5*Zeta3;
  be6[0][0][0][0][1][1][0] = 3.;
  be6[0][0][0][0][1][1][2] = -67.5;
  be6[0][0][0][0][1][2][0] = -6.75;
  be6[0][0][0][0][1][2][1] = 90.;
  be6[0][0][0][0][1][3][0] = 9.;
  be6[0][0][0][0][2][1][0] = -6.75;
  be6[0][0][0][0][2][1][1] = 45.;
  be6[0][0][0][0][2][2][0] = 34.875;
  be6[0][0][0][0][3][1][0] = 49.3125 + 9.*Zeta3;
  be6[0][0][0][1][0][1][0] = 3.;
  be6[0][0][0][1][0][1][2] = -67.5;
  be6[0][0][0][1][0][2][0] = -6.75;
  be6[0][0][0][1][0][2][1] = 90.;
  be6[0][0][0][1][0][3][0] = 9.;
  be6[0][0][0][1][1][1][0] = 1.5;
  be6[0][0][0][1][1][2][0] = -21.75;
  be6[0][0][0][1][2][1][0] = 51.9375;
  be6[0][0][0][2][0][1][0] = -6.75;
  be6[0][0][0][2][0][1][1] = 45.;
  be6[0][0][0][2][0][2][0] = 34.875;
  be6[0][0][0][2][1][1][0] = 51.9375;
  be6[0][0][0][3][0][1][0] = 49.3125 + 9.*Zeta3;
  be6[0][0][1][0][1][1][0] = 20.;
  be6[0][0][1][0][1][2][0] = -96. + 72.*Zeta3;
  be6[0][0][1][0][2][1][0] = 7.5 - 72.*Zeta3;
  be6[0][0][1][1][0][1][0] = 20.;
  be6[0][0][1][1][0][2][0] = -96. + 72.*Zeta3;
  be6[0][0][1][1][1][1][0] = 57. - 48.*Zeta3;
  be6[0][0][1][2][0][1][0] = 7.5 - 72.*Zeta3;
  be6[0][0][2][0][1][1][0] = 207.33333333333333333333333333333333333333333333333333333333333333333 - 24.*Zeta3;
  be6[0][0][2][1][0][1][0] = 207.33333333333333333333333333333333333333333333333333333333333333333 - 24.*Zeta3;
  be6[0][1][0][0][0][1][0] = -2.25;
  be6[0][1][0][0][0][1][2] = 45.;
  be6[0][1][0][0][0][2][0] = 10.3125;
  be6[0][1][0][0][0][2][1] = -67.5;
  be6[0][1][0][0][0][3][0] = -33.1875;
  be6[0][1][0][0][1][1][0] = 5.625;
  be6[0][1][0][0][1][2][0] = -33.75 - 27.*Zeta3;
  be6[0][1][0][0][2][1][0] = -72.5625 + 27.*Zeta3;
  be6[0][1][0][1][0][1][0] = 5.625;
  be6[0][1][0][1][0][2][0] = -33.75 - 27.*Zeta3;
  be6[0][1][0][1][1][1][0] = -48.375;
  be6[0][1][0][2][0][1][0] = -72.5625 + 27.*Zeta3;
  be6[0][1][1][0][1][1][0] = -122.25 + 108.*Zeta3;
  be6[0][1][1][1][0][1][0] = -122.25 + 108.*Zeta3;
  be6[0][2][0][0][0][1][0] = -5.75;
  be6[0][2][0][0][0][1][1] = -10.6875;
  be6[0][2][0][0][0][2][0] = 79.13671875 - 41.625*Zeta3;
  be6[0][2][0][0][1][1][0] = 71.0859375 - 74.25*Zeta3;
  be6[0][2][0][1][0][1][0] = 26.0859375 - 60.75*Zeta3;
  be6[0][2][1][0][0][1][0] = 87.75 - 108.*Zeta3;
  be6[0][3][0][0][0][1][0] = 0.78993055555555555555555555555555555555555555555555555555555555555556 + 140.625*Zeta3;
  be6[1][0][0][0][0][1][0] = -2.25;
  be6[1][0][0][0][0][1][2] = 9.;
  be6[1][0][0][0][0][2][0] = 6.7125;
  be6[1][0][0][0][0][2][1] = -9.9;
  be6[1][0][0][0][0][3][0] = -28.6875;
  be6[1][0][0][0][1][1][0] = 0.625;
  be6[1][0][0][0][1][2][0] = -4.35 - 1.8*Zeta3;
  be6[1][0][0][0][2][1][0] = -15.4125 + 5.4*Zeta3;
  be6[1][0][0][1][0][1][0] = 2.125;
  be6[1][0][0][1][0][2][0] = -11.55 + 3.6*Zeta3;
  be6[1][0][0][1][1][1][0] = -10.425 + 4.8*Zeta3;
  be6[1][0][0][2][0][1][0] = -11.9625 - 1.8*Zeta3;
  be6[1][0][1][0][1][1][0] = -16.516666666666666666666666666666666666666666666666666666666666666667 + 12.*Zeta3;
  be6[1][0][1][1][0][1][0] = -40.316666666666666666666666666666666666666666666666666666666666666667 + 40.8*Zeta3;
  be6[1][1][0][0][0][1][0] = 1.35;
  be6[1][1][0][0][0][1][1] = 6.525;
  be6[1][1][0][0][0][2][0] = -4.2609375 + 35.1*Zeta3;
  be6[1][1][0][0][1][1][0] = 17.090625 - 5.4*Zeta3;
  be6[1][1][0][1][0][1][0] = 10.678125 - 29.7*Zeta3;
  be6[1][2][0][0][0][1][0] = 9.196875 - 6.075*Zeta3;
  be6[2][0][0][0][0][1][0] = 6.855;
  be6[2][0][0][0][0][1][1] = -7.7625;
  be6[2][0][0][0][0][2][0] = -26.44453125 - 3.375*Zeta3;
  be6[2][0][0][0][1][1][0] = -5.6782291666666666666666666666666666666666666666666666666666666666667 - 0.87*Zeta3;
  be6[2][0][0][1][0][1][0] = -17.480729166666666666666666666666666666666666666666666666666666666667 - 34.71*Zeta3;
  be6[2][0][1][0][0][1][0] = 72.27 - 79.2*Zeta3;
  be6[2][1][0][0][0][1][0] = 4.336875 - 18.225*Zeta3;
  be6[3][0][0][0][0][1][0] = 75.907625 - 34.605*Zeta3;

  // 
  // Higgs self-coupling
  // 

  // * lam
  be7[0][0][0][0][0][0][2] = 12.;
  be7[0][0][0][0][0][0][3] = -156.;
  be7[0][0][0][0][0][0][4] = 3588. + 2016.*Zeta3;
  be7[0][0][0][0][0][1][1] = 2.;
  be7[0][0][0][0][0][1][2] = -24.;
  be7[0][0][0][0][0][1][3] = 291.;
  be7[0][0][0][0][0][2][0] = -1.;
  be7[0][0][0][0][0][2][1] = -0.5;
  be7[0][0][0][0][0][2][2] = 358.5 + 252.*Zeta3;
  be7[0][0][0][0][0][3][0] = 5.;
  be7[0][0][0][0][0][3][1] = -155.125 - 66.*Zeta3;
  be7[0][0][0][0][0][4][0] = -17.875 - 12.*Zeta3;
  be7[0][0][0][0][1][0][1] = 6.;
  be7[0][0][0][0][1][0][2] = -72.;
  be7[0][0][0][0][1][0][3] = 873.;
  be7[0][0][0][0][1][1][2] = -216.;
  be7[0][0][0][0][1][2][1] = 240.;
  be7[0][0][0][0][1][3][0] = -37.125;
  be7[0][0][0][0][2][0][0] = -3.;
  be7[0][0][0][0][2][0][1] = -1.5;
  be7[0][0][0][0][2][0][2] = 859.5 + 756.*Zeta3;
  be7[0][0][0][0][2][1][1] = 240.;
  be7[0][0][0][0][2][2][0] = -72.;
  be7[0][0][0][0][3][0][0] = 15.;
  be7[0][0][0][0][3][0][1] = 14.625 - 198.*Zeta3;
  be7[0][0][0][0][3][1][0] = -37.125;
  be7[0][0][0][0][4][0][0] = -199.875 - 36.*Zeta3;
  be7[0][0][0][1][0][0][1] = 6.;
  be7[0][0][0][1][0][0][2] = -72.;
  be7[0][0][0][1][0][0][3] = 873.;
  be7[0][0][0][1][0][1][2] = -216.;
  be7[0][0][0][1][0][2][1] = 240.;
  be7[0][0][0][1][0][3][0] = -37.125;
  be7[0][0][0][1][1][0][1] = -21.;
  be7[0][0][0][1][1][0][2] = 117. - 864.*Zeta3;
  be7[0][0][0][1][1][1][1] = 21.;
  be7[0][0][0][1][1][2][0] = 12.;
  be7[0][0][0][1][2][0][0] = -3.;
  be7[0][0][0][1][2][0][1] = 799.875 + 144.*Zeta3;
  be7[0][0][0][1][2][1][0] = 5.625;
  be7[0][0][0][1][3][0][0] = -89.625 - 36.*Zeta3;
  be7[0][0][0][2][0][0][0] = -3.;
  be7[0][0][0][2][0][0][1] = -1.5;
  be7[0][0][0][2][0][0][2] = 859.5 + 756.*Zeta3;
  be7[0][0][0][2][0][1][1] = 240.;
  be7[0][0][0][2][0][2][0] = -72.;
  be7[0][0][0][2][1][0][0] = -3.;
  be7[0][0][0][2][1][0][1] = 799.875 + 144.*Zeta3;
  be7[0][0][0][2][1][1][0] = 5.625;
  be7[0][0][0][2][2][0][0] = 72.*Zeta3;
  be7[0][0][0][3][0][0][0] = 15.;
  be7[0][0][0][3][0][0][1] = 14.625 - 198.*Zeta3;
  be7[0][0][0][3][0][1][0] = -37.125;
  be7[0][0][0][3][1][0][0] = -89.625 - 36.*Zeta3;
  be7[0][0][0][4][0][0][0] = -199.875 - 36.*Zeta3;
  be7[0][0][1][0][1][0][1] = 40.;
  be7[0][0][1][0][1][0][2] = -1224. + 1152.*Zeta3;
  be7[0][0][1][0][2][0][0] = -16.;
  be7[0][0][1][0][2][0][1] = 895. - 1296.*Zeta3;
  be7[0][0][1][0][3][0][0] = -38. + 240.*Zeta3;
  be7[0][0][1][1][0][0][1] = 40.;
  be7[0][0][1][1][0][0][2] = -1224. + 1152.*Zeta3;
  be7[0][0][1][1][1][0][1] = 82. - 96.*Zeta3;
  be7[0][0][1][1][2][0][0] = -2. - 48.*Zeta3;
  be7[0][0][1][2][0][0][0] = -16.;
  be7[0][0][1][2][0][0][1] = 895. - 1296.*Zeta3;
  be7[0][0][1][2][1][0][0] = -2. - 48.*Zeta3;
  be7[0][0][1][3][0][0][0] = -38. + 240.*Zeta3;
  be7[0][0][2][0][1][0][1] = 414.66666666666666666666666666666666666666666666666666666666666666667 - 48.*Zeta3;
  be7[0][0][2][0][2][0][0] = -88.666666666666666666666666666666666666666666666666666666666666666667 + 32.*Zeta3;
  be7[0][0][2][1][0][0][1] = 414.66666666666666666666666666666666666666666666666666666666666666667 - 48.*Zeta3;
  be7[0][0][2][1][1][0][0] = 192.;
  be7[0][0][2][2][0][0][0] = -88.666666666666666666666666666666666666666666666666666666666666666667 + 32.*Zeta3;
  be7[0][1][0][0][0][0][1] = -4.5;
  be7[0][1][0][0][0][0][2] = 54.;
  be7[0][1][0][0][0][0][3] = -474. + 72.*Zeta3;
  be7[0][1][0][0][0][1][1] = 3.75;
  be7[0][1][0][0][0][1][2] = 53.25 - 144.*Zeta3;
  be7[0][1][0][0][0][2][1] = -198.375 + 171.*Zeta3;
  be7[0][1][0][0][0][3][0] = 35.53125 - 9.*Zeta3;
  be7[0][1][0][0][1][0][1] = 11.25;
  be7[0][1][0][0][1][0][2] = 159.75 - 432.*Zeta3;
  be7[0][1][0][0][1][1][1] = -27.;
  be7[0][1][0][0][2][0][1] = -622.125 + 513.*Zeta3;
  be7[0][1][0][0][3][0][0] = 106.59375 - 27.*Zeta3;
  be7[0][1][0][1][0][0][1] = 11.25;
  be7[0][1][0][1][0][0][2] = 159.75 - 432.*Zeta3;
  be7[0][1][0][1][0][1][1] = -27.;
  be7[0][1][0][1][1][0][1] = -132.75 + 54.*Zeta3;
  be7[0][1][0][1][2][0][0] = 14.90625;
  be7[0][1][0][2][0][0][1] = -622.125 + 513.*Zeta3;
  be7[0][1][0][2][1][0][0] = 14.90625;
  be7[0][1][0][3][0][0][0] = 106.59375 - 27.*Zeta3;
  be7[0][1][1][0][1][0][1] = -244.5 + 216.*Zeta3;
  be7[0][1][1][0][2][0][0] = -15.5 + 24.*Zeta3;
  be7[0][1][1][1][0][0][1] = -244.5 + 216.*Zeta3;
  be7[0][1][1][1][1][0][0] = -8. + 96.*Zeta3;
  be7[0][1][1][2][0][0][0] = -15.5 + 24.*Zeta3;
  be7[0][2][0][0][0][0][0] = 0.5625;
  be7[0][2][0][0][0][0][1] = -4.5625;
  be7[0][2][0][0][0][0][2] = -173.625 - 513.*Zeta3;
  be7[0][2][0][0][0][1][0] = -0.375;
  be7[0][2][0][0][0][1][1] = -36.234375 - 58.5*Zeta3;
  be7[0][2][0][0][0][2][0] = 25.4296875 - 17.0625*Zeta3;
  be7[0][2][0][0][1][0][0] = -1.125;
  be7[0][2][0][0][1][0][1] = -108.703125 - 175.5*Zeta3;
  be7[0][2][0][0][1][1][0] = 1.125;
  be7[0][2][0][0][2][0][0] = 77.4140625 - 51.1875*Zeta3;
  be7[0][2][0][1][0][0][0] = -1.125;
  be7[0][2][0][1][0][0][1] = -108.703125 - 175.5*Zeta3;
  be7[0][2][0][1][0][1][0] = 1.125;
  be7[0][2][0][1][1][0][0] = -41.484375 + 58.5*Zeta3;
  be7[0][2][0][2][0][0][0] = 77.4140625 - 51.1875*Zeta3;
  be7[0][2][1][0][0][0][1] = 202.5 - 216.*Zeta3;
  be7[0][2][1][0][1][0][0] = 81.375 - 54.*Zeta3;
  be7[0][2][1][1][0][0][0] = 81.375 - 54.*Zeta3;
  be7[0][3][0][0][0][0][0] = 9.53125;
  be7[0][3][0][0][0][0][1] = 201.49652777777777777777777777777777777777777777777777777777777777778 + 552.375*Zeta3;
  be7[0][3][0][0][0][1][0] = -8.91796875 + 24.75*Zeta3;
  be7[0][3][0][0][1][0][0] = -26.75390625 + 74.25*Zeta3;
  be7[0][3][0][1][0][0][0] = -26.75390625 + 74.25*Zeta3;
  be7[0][3][1][0][0][0][0] = -57.375 + 54.*Zeta3;
  be7[0][4][0][0][0][0][0] = 74.303059895833333333333333333333333333333333333333333333333333333333 - 156.7265625*Zeta3;
  be7[1][0][0][0][0][0][1] = -0.9;
  be7[1][0][0][0][0][0][2] = 10.8;
  be7[1][0][0][0][0][0][3] = -94.8 + 14.4*Zeta3;
  be7[1][0][0][0][0][1][1] = 3.75;
  be7[1][0][0][0][0][1][2] = -81.15 + 57.6*Zeta3;
  be7[1][0][0][0][0][2][0] = -1.2;
  be7[1][0][0][0][0][2][1] = 38.025 - 70.2*Zeta3;
  be7[1][0][0][0][0][3][0] = 2.53125 + 19.8*Zeta3;
  be7[1][0][0][0][1][0][1] = 1.25;
  be7[1][0][0][0][1][0][2] = 62.55 - 115.2*Zeta3;
  be7[1][0][0][0][1][1][1] = -5.4;
  be7[1][0][0][0][2][0][0] = 0.4;
  be7[1][0][0][0][2][0][1] = -143.425 + 149.4*Zeta3;
  be7[1][0][0][0][3][0][0] = 31.94375 - 15.*Zeta3;
  be7[1][0][0][1][0][0][1] = 4.25;
  be7[1][0][0][1][0][0][2] = -29.25 - 28.8*Zeta3;
  be7[1][0][0][1][0][1][1] = -5.4;
  be7[1][0][0][1][1][0][1] = -46.45 - 1.2*Zeta3;
  be7[1][0][0][1][2][0][0] = -14.36875 + 15.6*Zeta3;
  be7[1][0][0][2][0][0][0] = -0.8;
  be7[1][0][0][2][0][0][1] = -62.125 + 34.2*Zeta3;
  be7[1][0][0][2][1][0][0] = 8.35625 - 16.8*Zeta3;
  be7[1][0][0][3][0][0][0] = 21.66875 + 10.2*Zeta3;
  be7[1][0][1][0][1][0][1] = -33.033333333333333333333333333333333333333333333333333333333333333333 + 24.*Zeta3;
  be7[1][0][1][0][2][0][0] = -21.366666666666666666666666666666666666666666666666666666666666666667 + 27.2*Zeta3;
  be7[1][0][1][1][0][0][1] = -80.633333333333333333333333333333333333333333333333333333333333333333 + 81.6*Zeta3;
  be7[1][0][1][2][0][0][0] = 31.033333333333333333333333333333333333333333333333333333333333333333 - 11.2*Zeta3;
  be7[1][1][0][0][0][0][0] = 0.225;
  be7[1][1][0][0][0][0][1] = 2.925;
  be7[1][1][0][0][0][0][2] = -199.8 - 97.2*Zeta3;
  be7[1][1][0][0][0][1][0] = 1.65;
  be7[1][1][0][0][0][1][1] = -70.70625 + 75.6*Zeta3;
  be7[1][1][0][0][0][2][0] = -0.140625 - 28.575*Zeta3;
  be7[1][1][0][0][1][0][0] = 1.35;
  be7[1][1][0][0][1][0][1] = -56.41875 + 7.2*Zeta3;
  be7[1][1][0][0][1][1][0] = -0.75;
  be7[1][1][0][0][2][0][0] = -10.121875 - 23.325*Zeta3;
  be7[1][1][0][1][0][0][0] = 3.15;
  be7[1][1][0][1][0][0][1] = -122.04375 + 106.2*Zeta3;
  be7[1][1][0][1][0][1][0] = 4.35;
  be7[1][1][0][1][1][0][0] = 6.25625 + 9.3*Zeta3;
  be7[1][1][0][2][0][0][0] = -3.371875 - 55.725*Zeta3;
  be7[1][1][1][0][1][0][0] = 34.95 - 21.6*Zeta3;
  be7[1][1][1][1][0][0][0] = 37.35 - 21.6*Zeta3;
  be7[1][2][0][0][0][0][0] = -1.80625;
  be7[1][2][0][0][0][0][1] = 115.06875 - 29.475*Zeta3;
  be7[1][2][0][0][0][1][0] = 3.39609375 - 0.9*Zeta3;
  be7[1][2][0][0][1][0][0] = 14.29453125 + 5.4*Zeta3;
  be7[1][2][0][1][0][0][0] = 8.17265625 + 4.05*Zeta3;
  be7[1][2][1][0][0][0][0] = -11.475 + 10.8*Zeta3;
  be7[1][3][0][0][0][0][0] = -28.761284722222222222222222222222222222222222222222222222222222222222 - 7.59375*Zeta3;
  be7[2][0][0][0][0][0][0] = 0.0675;
  be7[2][0][0][0][0][0][1] = 4.7175;
  be7[2][0][0][0][0][0][2] = -150.48 - 29.16*Zeta3;
  be7[2][0][0][0][0][1][0] = -1.125;
  be7[2][0][0][0][0][1][1] = -27.579375 - 22.14*Zeta3;
  be7[2][0][0][0][0][2][0] = 21.8728125 + 8.4375*Zeta3;
  be7[2][0][0][0][1][0][0] = 0.225;
  be7[2][0][0][0][1][0][1] = -31.171458333333333333333333333333333333333333333333333333333333333333 - 2.82*Zeta3;
  be7[2][0][0][0][1][1][0] = 0.615;
  be7[2][0][0][0][2][0][0] = -10.873229166666666666666666666666666666666666666666666666666666666667 - 5.0875*Zeta3;
  be7[2][0][0][1][0][0][0] = -0.855;
  be7[2][0][0][1][0][0][1] = -42.476458333333333333333333333333333333333333333333333333333333333333 - 26.94*Zeta3;
  be7[2][0][0][1][0][1][0] = 10.515;
  be7[2][0][0][1][1][0][0] = -3.988125 - 0.36*Zeta3;
  be7[2][0][0][2][0][0][0] = 7.0617708333333333333333333333333333333333333333333333333333333333333 + 7.3925*Zeta3;
  be7[2][0][1][0][0][0][1] = 29.7 - 31.68*Zeta3;
  be7[2][0][1][0][1][0][0] = 10.245 - 6.48*Zeta3;
  be7[2][0][1][1][0][0][0] = 8.805 - 6.48*Zeta3;
  be7[2][1][0][0][0][0][0] = -2.09625;
  be7[2][1][0][0][0][0][1] = 69.705 - 6.615*Zeta3;
  be7[2][1][0][0][0][1][0] = 8.46140625 - 2.7*Zeta3;
  be7[2][1][0][0][1][0][0] = 10.69171875 + 1.62*Zeta3;
  be7[2][1][0][1][0][0][0] = 11.92546875 - 1.08*Zeta3;
  be7[2][1][1][0][0][0][0] = -8.415 + 7.92*Zeta3;
  be7[2][2][0][0][0][0][0] = -8.4905208333333333333333333333333333333333333333333333333333333333333 + 12.470625*Zeta3;
  be7[3][0][0][0][0][0][0] = -0.85275;
  be7[3][0][0][0][0][0][1] = 44.3195 - 13.437*Zeta3;
  be7[3][0][0][0][0][1][0] = 9.25509375 - 0.81*Zeta3;
  be7[3][0][0][0][1][0][0] = 4.54903125 + 0.27*Zeta3;
  be7[3][0][0][1][0][0][0] = 11.76590625 - 0.54*Zeta3;
  be7[3][0][1][0][0][0][0] = -5.049 + 4.752*Zeta3;
  be7[3][1][0][0][0][0][0] = -7.43084375 + 4.89825*Zeta3;
  be7[4][0][0][0][0][0][0] = -6.5616328125 + 4.2042375*Zeta3;
}


void BetaSMFull::operator() (const state_type &a, state_type &dadt, const double t)
{
  // 
  //    - couplings: a1->a[0], a2->a[1], as->a[2], at->a[3], ab->a[4], atau->a[5], lam->a[6]
  // 
  //    - parameter: t = Log(mu/mu0)
  // 
    
  dadt[0] = 0;                // beta a1
  dadt[1] = 0;                // beta a2
  dadt[2] = 0;                // beta as
  dadt[3] = 0;                // beta at
  dadt[4] = 0;                // beta ab
  dadt[5] = 0;                // beta atau
  dadt[6] = 0;                // beta lam
    
  long double a1   = pocoa1   < 0  ? 0 : a[0];
  long double a2   = pocoa2   < 0  ? 0 : a[1];
  long double as   = pocoas   < 0  ? 0 : a[2];
  long double at   = pocoat   < 0  ? 0 : a[3];
  long double ab   = pocoab   < 0  ? 0 : a[4];
  long double atau = pocoatau < 0  ? 0 : a[5];
  long double lam  = pocolam  < 0  ? 0 : a[6];

  for(size_t ia1 = 0; ia1 < maxPoco; ia1++)
    for(size_t ia2 = 0; ia2 < maxPoco - ia1; ia2++)
      for(size_t ias = 0; ias < maxPoco - ia1 - ia2; ias++)
        for(size_t iat = 0; iat < maxPoco - ia1 - ia2 - ias; iat++)
          for(size_t iab = 0; iab < maxPoco - ia1 - ia2 - ias - iat; iab++)
            for(size_t iatau = 0; iatau < maxPoco - ia1 - ia2 - ias - iat - iab; iatau++)
              for(size_t ilam = 0; ilam < maxPoco - ia1 - ia2 - ias - iat - iab - iatau; ilam++)
                {
                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoa1 + 2) 
                    dadt[0] += be1[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);

                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoa2 + 2) 
                    dadt[1] += be2[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);
                    
                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoas + 2) 
                    dadt[2] += be3[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);

                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoat + 2) 
                    dadt[3] += be4[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);

                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoab + 2) 
                    dadt[4] += be5[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);

                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocoatau + 2) 
                    dadt[5] += be6[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);
                    
                  if(ia1+ia2+ias+iat+iab+iatau+ilam <= pocolam + 2) 
                    dadt[6] += be7[ia1][ia2][ias][iat][iab][iatau][ilam]
                      *pow(a1,ia1)*pow(a2,ia2)*pow(as,ias)
                      *pow(at,iat)*pow(ab,iab)*pow(atau,iatau)
                      *pow(lam,ilam);
                }
}

long double BetaSMFull::betaQCD(long double as)
{
  long double bQCD = 0;
  for(size_t i = 0; i < 5; i++)
    bQCD += be3[0][0][i][0][0][0][0]*pow(as,i);
  
  std::cout << "bQCD: " << bQCD << std::endl;
    
}


