:Evaluate:  Print["MR - Matching and Running Mathmatica interface"];
:Evaluate:  Print["To see available functions use Names[\"mr`*\"]"];
:Evaluate:  Print["Andrey Pikelner <pikelner@theor.jinr.ru>"];

:Evaluate:  BeginPackage["mr`"]

:Evaluate:  RunQCD::usage  = "RunQCD[oscale,asMZ, MZscale,L=4, mtth] evaluate alphas at specified oscale given asMZ at MZscale (nf=5) with L-loop RGE with top threshold at mtth"
:Evaluate:  RunSM::usage  = "RunSM[gp,g,gs,yt,lam,m,iscale,oscale] return running parameters at specified oscale given the values at specified iscale"


:Evaluate:  MW::usage  = "MW[gp,g,gs,yb,yt,lam,m,scale] returns pole W-boson mass MW given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MWp::usage  = "MWp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole W-boson mass MW given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MZp::usage  = "MZp[gp,g,gs,yb,yt,lam,m,scale] returns pole Z-boson mass MZ given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MZ::usage  = "MZ[gp,g,gs,yb,yt,lam,m,scale] returns pole Z-boson mass MZ given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MH::usage  = "MH[gp,g,gs,yb,yt,lam,m,scale] returns pole H-boson mass MH given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MHp::usage  = "MHp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole H-boson mass MH given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MT::usage  = "MT[gp,g,gs,yb,yt,lam,m,scale] returns pole t-quark mass MT given  MSbar parameters at specified scale at 2-loop level + 3-loop QCD"
:Evaluate:  MTp::usage  = "MTp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole t-quark mass MT given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  GF::usage  = "GF[gp,g,gs,yb,yt,lam,m,scale] returns Fermi constant GF given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  GFp::usage  = "GFp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns Fermi constant GF given  MSbar parameters at specified scale at L-loop level"


:Evaluate:  XMMW::usage  = "XMMW[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole W-boson mass MW^2 given  MSbar parameters at specified scale"
:Evaluate:  XMMZ::usage  = "XMMZ[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole Z-boson mass MZ^2 given  MSbar parameters at specified scale"
:Evaluate:  XMMH::usage  = "XMMH[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole H-boson mass MH^2 given  MSbar parameters at specified scale"
:Evaluate:  XdRbar::usage  = "XdRbar[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the Fermi constant GF - running vev relation given  MSbar parameters at specified scale"

:Evaluate:  XMT::usage  = "XMT[gp,g,gs,yb,yt,lam,m,scale] returns electroweak contributions to the pole top quark mass MT given  MSbar parameters at specified scale"
:Evaluate:  XMTQCD::usage  = "XMT[gp,g,gs,yb,yt,lam,m,scale] returns pure QCD contributions to the pole top quark mass MT given  MSbar parameters at specified scale"


:Evaluate:  Xb::usage  = "Xb[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XW::usage  = "XW[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XZ::usage  = "XZ[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XH::usage  = "XH[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  Xt::usage  = "Xt[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  dROS::usage  = "dROS[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XbQCD::usage  = "XbQCD[Mb,MW,MZ,MH,Mt,scale,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 
:Evaluate:  XtQCD::usage  = "XtQCD[Mb,MW,MZ,MH,Mt,scale,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 
:Evaluate:  mmWMMW::usage  = "mmWMMW[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mW^2/MW^2, full correction to relation between MS-bar mass mW and pole MW"
:Evaluate:  mmZMMZ::usage  = "mmZMMZ[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mZ^2/MZ^2, full correction to relation between MS-bar mass mZ and pole MZ"
:Evaluate:  mmHMMH::usage  = "mmHMMH[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mH^2/MH^2, full correction to relation between MS-bar mass mH and pole MH"
:Evaluate:  mtMt::usage  = "mtMt[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mt/Mt, full correction to relation between MS-bar mass mt and pole Mt"
:Evaluate:  mbMb::usage  = "mbMb[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mb/Mb, full correction to relation between MS-bar mass mb and pole Mb"


:Evaluate:  a1::usage  = "a1[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns a1=5/3(g1/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2, we use GUT normalization"
:Evaluate:  a2::usage  = "a2[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns a2=(g2/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  at::usage  = "a1[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns at=(yt/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  ab::usage  = "ab[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns ab=(yb/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  alam::usage  = "alam[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] returns alam=\[Lambda]/(4\[Pi])^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  vev::usage  = "vev[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns running vev[scale] as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"

:Evaluate:  dalphaGF::usage  = "dalphaGF[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] - calculates loop contributions to the relation between GF and running alpha[scale]"
:Evaluate:  alphaGF::usage  = "alphaGF[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] - calculates running alpha[scale] in terms GF and running alpha_s[scale]" 

:Evaluate:  aEW::usage  = "\[Alpha]/(4\[Pi]), running EW constant"
:Evaluate:  aQCD::usage  = "\[Alpha]_S/(4\[Pi]), running QCD constant"

:Evaluate:  Gf::usage  = "G_F Fermi constant"
:Evaluate:  h::usage  = "h - loop counter"

:Evaluate:  xW::usage  = "xW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mW and the pole mass MW"
:Evaluate:  xZ::usage  = "xZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mZ and the pole mass MZ"
:Evaluate:  xH::usage  = "xH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mH and the pole mass MH"
:Evaluate:  xt::usage  = "xt[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mt and the pole mass Mt"
:Evaluate:  xb::usage  = "xb[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mb and the pole mass Mb"

:Evaluate:  yW::usage  = "yW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to W-boson and the pole mass MW"
:Evaluate:  yZ::usage  = "yZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to Z-boson and the pole mass MZ"
:Evaluate:  yH::usage  = "yH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to H-boson and the pole mass MH"
:Evaluate:  yt::usage  = "yt[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to t-quark and the pole mass Mt"
:Evaluate:  yb::usage  = "yb[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to b-quark and the pole mass Mb"


:Evaluate:  xMMW::usage  = "xMW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MW^2 and running parameters"
:Evaluate:  xMMZ::usage  = "xMZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MZ^2 and running parameters"
:Evaluate:  xMMH::usage  = "xMH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MH^2 and running parameters" 
:Evaluate:  xMT::usage  = "xMT[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MT and running parameters"
:Evaluate:  xMTQCD::usage  = "xMTQCD[0,b] represents a coefficient of pure QCD contribution aQCD^b in the relation between the pole mass MT and running parameters"
:Evaluate:  xdRbar::usage  = "xdRbar[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between GF and running parameters (vev, etc)"



:Evaluate:  dr::usage  = "dr[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs vev and the Fermi constant GF"

:Evaluate:  daGF::usage  = "daGF[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running electromagnetic alpha and Fermi constant GF. Note that aEW should be again expressed in terms of Fermi constant"

:Evaluate:   g1::usage  = "running U(1) coupling"
:Evaluate:   g2::usage  = "running SU(2) coupling"
:Evaluate:   gs::usage  = "running SU(3) strong coupling"
:Evaluate:   yt::usage  = "running top Yukawa coupling"
:Evaluate:   yb::usage  = "running bottom Yukawa coupling"
:Evaluate:   lam::usage  = "running higgs self-coupling"
:Evaluate:   m::usage    = "running higgs mass parameter"
:Evaluate:   scale::usage  = "renormalization scale"

:Evaluate:  Protect[g1,g2,gs,yt,yb,lam,m,scale];

:Evaluate:  Begin["`Private`"]

// Mathematica part

// Masses
:Evaluate:  mmWMMW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[scale]*xW[1,0]+aEW[scale]*aQCD[scale]*xW[1,1]+aEW[scale]^2*xW[2,0])/.XW[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  mmZMMZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[scale]*xZ[1,0]+aEW[scale]*aQCD[scale]*xZ[1,1]+aEW[scale]^2*xZ[2,0])/.XZ[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  mmHMMH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[scale]*xH[1,0]+aEW[scale]*aQCD[scale]*xH[1,1]+aEW[scale]^2*xH[2,0])/.XH[mb,mW,mZ,mH,mt,scale,nL,nH];
                                                                                                                             
:Evaluate:  mtMt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[scale]*xt[1,0]+aEW[scale]*aQCD[scale]*xt[1,1]+aEW[scale]^2*xt[2,0])/.Xt[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  mbMb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[scale]*xb[1,0]+aEW[scale]*aQCD[scale]*xb[1,1]+aEW[scale]^2*xb[2,0])/.Xb[mb,mW,mZ,mH,mt,scale,nL,nH];

// Couplings

:Evaluate:  a1[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 5/3*2^(5/2)*Gf/(4*Pi)^2*((mZ^2*(1+aEW[scale]*yZ[1,0]+aEW[scale]*aQCD[scale]*yZ[1,1]+aEW[scale]^2*yZ[2,0])/.XZ[mb,mW,mZ,mH,mt,scale,nL,nH]) - (mW^2*(1+aEW[scale]*yW[1,0]+aEW[scale]*aQCD[scale]*yW[1,1]+aEW[scale]^2*yW[2,0])/.XW[mb,mW,mZ,mH,mt,scale,nL,nH]));

:Evaluate:  a2[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(5/2)*Gf/(4*Pi)^2*mW^2*(1+aEW[scale]*yW[1,0]+aEW[scale]*aQCD[scale]*yW[1,1]+aEW[scale]^2*yW[2,0])/.XW[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  at[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(3/4)*Sqrt[Gf]/(4*Pi)^2*mt*(1+aQCD[scale]*xt[0,1] + aEW[scale]*yt[1,0]+aEW[scale]*aQCD[scale]*yt[1,1]+aEW[scale]^2*yt[2,0] + aQCD[scale]^2*xt[0,2] + aQCD[scale]^3*xt[0,3])/.Xt[mb,mW,mZ,mH,mt,scale,nL,nH] /. XtQCD[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  ab[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(3/4)*Sqrt[Gf]/(4*Pi)^2*mb*(1+aQCD[scale]*xb[0,1]+aEW[scale]*yb[1,0]+aEW[scale]*aQCD[scale]*yb[1,1]+aEW[scale]^2*yb[2,0] + aQCD[scale]^2*xb[0,2] + aQCD[scale]^3*xb[0,3])/.Xb[mb,mW,mZ,mH,mt,scale,nL,nH] /. XbQCD[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  alam[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(-1/2)*Gf/(4*Pi)^2*mH^2*(1+aEW[scale]*yH[1,0]+aEW[scale]*aQCD[scale]*yH[1,1]+aEW[scale]^2*yH[2,0])/.XH[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  alphaGF[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := Module[{alGF = 2^(1/2)*Gf*mW^2*(1-mW^2/mZ^2)/Pi},alGF*(1 + h alGF/(4 Pi)*daGF[1,0] + h^2 alGF/(4 Pi)*aQCD[scale]*daGF[1,1]+ h^2 alGF^2/(4 Pi)^2*daGF[2,0])/.dalphaGF[mb,mW,mZ,mH,mt,scale,nL,nH]];


:Evaluate:  vev[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(-1/4)/Sqrt[Gf]*Sqrt[(1+aEW[scale]*dr[1,0]+aEW[scale]*aQCD[scale]*dr[1,1]+aEW[scale]^2*dr[2,0])] /.dROS[mb,mW,mZ,mH,mt,scale,nL,nH];

:Evaluate:  MW[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, looptag_:1] := Block[{vev,lc,aEW,aQCD,mW},
			(* loop corrections *)	lc = XMMW[ Sequence @@ pars];
			vev = Sqrt[M^2/LAM/2]; mW = G2 * vev/2.0;
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			Return[ mW * Sqrt[1 + looptag * aEW * xMMW[1, 0] + looptag^2 ( aEW^2 * xMMW[2,0] + aEW*aQCD * xMMW[1,1])] /. lc ]]	


:Evaluate:  MW[runpars_List, looptag_:1]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MW[ Sequence @@ pars, looptag ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];


:Evaluate:  MWold[runpars_List]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, vev, lc, aEW, aQCD, mW}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			(* loop corrections *)	lc = XMMW[ Sequence @@ pars];
			vev = Sqrt[m^2/lam/2] /. runpars;
			mW = g2 * vev/2.0 /. runpars;
			{aEW, aQCD} = {g1^2*g2^2/(g1^2 + g2^2), gs^2}/(16 Pi^2) /. runpars; 
			Return[ mW * Sqrt[1 + h * aEW * xMMW[1, 0] + h^2 ( aEW^2 * xMMW[2,0] + aEW*aQCD * xMMW[1,1])] /. lc ],	
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];

:Evaluate:  MZ[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, looptag_:1] := Block[{vev,lc,aEW,aQCD,mZ},
			(* loop corrections *)	lc = XMMZ[ Sequence @@ pars];
			vev = Sqrt[M^2/LAM/2]; mZ = Sqrt[G1^2 + G2^2] * vev/2.0;
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			Return[ mZ * Sqrt[1 + looptag * aEW * xMMZ[1, 0] + looptag^2 ( aEW^2 * xMMZ[2,0] + aEW*aQCD * xMMZ[1,1])] /. lc ]]	
:Evaluate:  MZ[runpars_List, looptag_:1]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MZ[ Sequence @@ pars, looptag ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];

:Evaluate:  MZold[runpars_List]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, vev, lc, aEW, aQCD, mZ}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			(* loop corrections *)	lc = XMMZ[ Sequence @@ pars];
			vev = Sqrt[m^2/lam/2] /. runpars;
			mZ = Sqrt[g1^2 + g2^2] * vev/2.0 /. runpars;
			{aEW, aQCD} = {g1^2*g2^2/(g1^2 + g2^2), gs^2}/(16 Pi^2) /. runpars; 
			Return[ mZ * Sqrt[1 + h * aEW * xMMZ[1, 0] + h^2 ( aEW^2 * xMMZ[2,0] + aEW*aQCD * xMMZ[1,1])] /. lc ],	
			(* else *) Print[" Not All parameters specified  ", pars, " from ", runpars]]];


:Evaluate:  MH[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, looptag_:1] := Block[{vev,lc,aEW,aQCD},
			(* loop corrections *)	lc = XMMH[ Sequence @@ pars];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			Return[ M * Sqrt[1 + looptag * aEW * xMMH[1, 0] + looptag^2 ( aEW^2 * xMMH[2,0] + aEW*aQCD * xMMH[1,1])] /. lc ]]	

:Evaluate:  MH[runpars_List, looptag_:1]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MH[ Sequence @@ pars, looptag ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];


:Evaluate:  MHold[runpars_List]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, vev, lc, aEW, aQCD, mH}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			(* loop corrections *)	lc = XMMH[ Sequence @@ pars];
			mH = m /. runpars; (* notation *)
			{aEW, aQCD} = {g1^2*g2^2/(g1^2 + g2^2), gs^2}/(16 Pi^2) /. runpars; 
			Return[ mH * Sqrt[1 + h * aEW * xMMH[1, 0] + h^2 ( aEW^2 * xMMH[2,0] + aEW*aQCD * xMMH[1,1])] /. lc ],	
			(* else *) Print[" Not All parameters specified  ", pars, " from ", runpars]]];


:Evaluate:  MT[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, looptag_:1] := Block[{vev,lc,aEW,aQCD,mt},
			(* loop corrections *)	lc = Join[XMT[ Sequence @@ pars], XMTQCD[ Sequence @@ pars]];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			vev = Sqrt[M^2/LAM/2]; mt = YT * vev / Sqrt[2];
			Return[ mt * (1 + looptag * (aQCD * xMTQCD[0,1] + aEW * xMT[1, 0]) + looptag^2 ( aQCD^2 * xMTQCD[0,2] + aEW^2 * xMT[2,0] + aEW*aQCD * xMT[1,1] ) + looptag^3 * aQCD^3 * xMTQCD[0,3]) /. lc ]]	

:Evaluate:  MT[runpars_List, looptag_:1]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MT[ Sequence @@ pars, looptag ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];


:Evaluate:  MTold[runpars_List]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, vev, lc, aEW, aQCD, mt}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			(* loop corrections *)	lc = Join[ XMT[ Sequence @@ pars], XMTQCD[ Sequence @@ pars]];
			vev = Sqrt[m^2/lam/2] /. runpars;
			mt = yt * vev/ Sqrt[2] /. runpars; 
			{aEW, aQCD} = {g1^2*g2^2/(g1^2 + g2^2), gs^2}/(16 Pi^2) /. runpars; 
			Return[ mt * (1 + h * (aQCD * xMTQCD[0,1] + aEW * xMT[1, 0]) + h^2 ( aQCD^2 * xMTQCD[0,2] + aEW^2 * xMT[2,0] + aEW*aQCD * xMT[1,1] ) + h^3 * aQCD^3 * xMTQCD[0,3]) /. lc ],	
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars ]]];

:Evaluate:  GF[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, looptag_:1] := Block[{vv,lc,aEW,aQCD,gf},
			(* loop corrections *)	lc = XdRbar[ Sequence @@ pars];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			vv = M^2/LAM/2; gf = 1/Sqrt[2]/vv; 
			Return[ gf * (1 + looptag * aEW * xdRbar[1, 0] + looptag^2 (  aEW^2 * xdRbar[2,0] + aEW*aQCD * xdRbar[1,1] ) ) /. lc ]];	
				
:Evaluate:  GF[runpars_List, looptag_:1]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ GF[ Sequence @@ pars, looptag ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];

:Evaluate:  GFold[runpars_List]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, vv, lc, aEW, aQCD, gf}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			(* loop corrections *)	lc = XdRbar[ Sequence @@ pars];
			vv = m^2/lam/2 /. runpars;
			gf = 1/Sqrt[2]/vv /. runpars; 
			{aEW, aQCD} = {g1^2*g2^2/(g1^2 + g2^2), gs^2}/(16 Pi^2) /. runpars; 
			Return[ gf * (1 + h * aEW * xdRbar[1, 0] + h^2 (  aEW^2 * xdRbar[2,0] + aEW*aQCD * xdRbar[1,1] ) ) /. lc ],	
			(* else *) Print[" Not All parameters specified " , pars, " from ", runpars ]]];



// C++ part
:Begin:
:Function: XMMW
:Pattern: XMMW[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: XMMZ
:Pattern: XMMZ[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XMMH
:Pattern: XMMH[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XdRbar
:Pattern: XdRbar[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: XMT
:Pattern: XMT[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XMTQCD
:Pattern: XMTQCD[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: MWp
:Pattern: MWp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:


:Begin:
:Function: MZp
:Pattern: MZp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MHp
:Pattern: MHp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:


:Begin:
:Function: MTp
:Pattern: MTp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: GFp
:Pattern: GFp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: RunQCD
:Pattern: RunQCD[oscale_?NumericQ,asMZ_?NumericQ,MZscale_?NumericQ,nL_Integer,mtth_?NumericQ]
:Arguments: {N[oscale],N[asMZ],N[MZscale],nL,N[mtth]}
:ArgumentTypes: {Real128,Real128,Real128,Integer,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: RunSM
:Pattern: RunSM[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,iscale_?NumericQ,oscale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[iscale],N[oscale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XW
:Pattern: XW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XZ
:Pattern: XZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XH
:Pattern: XH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xt
:Pattern: Xt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: dROS
:Pattern: dROS[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: dalphaGF
:Pattern: dalphaGF[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xb
:Pattern: Xb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XtQCD
:Pattern: XtQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XbQCD
:Pattern: XbQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:


:Evaluate: End[]
:Evaluate: EndPackage[]
