:Evaluate:  Print["MR - Matching and Running Mathmatica interface"];
:Evaluate:  Print["To see available functions use Names[\"mr`*\"]"];
:Evaluate:  Print["Andrey Pikelner <pikelner@theor.jinr.ru>"];

:Evaluate:  BeginPackage["mr`"]

:Evaluate:  RunQCD::usage  = "RunQCD[oscale,asMZ, MZscale,L=4, mtth] evaluate alphas at specified oscale given asMZ at MZscale (nf=5) with L-loop RGE with top threshold at mtth"
:Evaluate:  RunSM::usage  = "RunSM[gp,g,gs,yt,lam,m,iscale,oscale] return running parameters at specified oscale given the values at specified iscale"


:Evaluate:  MW::usage  = "MW[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole W-boson mass MW given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MW::usage  = "MW[runpars,L=2] returns pole W-boson mass MW given  MSbar parameters (runpars = { scale -> ..., g1 -> .., etc }) at L-loop level"
:Evaluate:  MWp::usage  = "MWp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole W-boson mass MW given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MZp::usage  = "MZp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole Z-boson mass MZ given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MZ::usage  = "MZ[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole Z-boson mass MZ given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MH::usage  = "MH[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole H-boson mass MH given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MHp::usage  = "MHp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole H-boson mass MH given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MT::usage  = "MT[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole t-quark mass MT given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  MTp::usage  = "MTp[gp,g,gs,yb,yt,lam,m,scale,L=2] returns pole t-quark mass MT given  MSbar parameters at specified scale at L-loop level"
:Evaluate:  GF::usage  = "GF[gp,g,gs,yb,yt,lam,m,scale,L=2] returns Fermi constant GF given  MSbar parameters at specified scale at L-loop level"


:Evaluate:  XMW::usage  = "XMW[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole W-boson mass MW given  MSbar parameters at specified scale"
:Evaluate:  XMZ::usage  = "XMZ[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole Z-boson mass MZ given  MSbar parameters at specified scale"
:Evaluate:  XMH::usage  = "XMH[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole H-boson mass MH given  MSbar parameters at specified scale"
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

:Evaluate:  dr::usage  = "dr[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs vev and the Fermi constant GF"

:Evaluate:  daGF::usage  = "daGF[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running electromagnetic alpha and Fermi constant GF. Note that aEW should be again expressed in terms of Fermi constant"

:Evaluate:   g1::usage  = "running U(1) coupling"
:Evaluate:   g2::usage  = "running SU(2) coupling"
:Evaluate:   gs::usage  = "running SU(3) strong coupling"
:Evaluate:   yt::usage  = "running top Yukawa coupling"
:Evaluate:   yb::usage  = "running bottom Yukawa coupling"
:Evaluate:   lam::usage  = "running higgs self-coupling"
:Evaluate:   m::usage    = "running higgs mass parameter"


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

:Evaluate:  MW[runpars_List,L_?NumericQ]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
							Return[ Sequence @@ pars, L], 
			(* else *) Print[" Not All parameters specified "]]];

// C++ part
:Begin:
:Function: XMW
:Pattern: XMW[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: XMZ
:Pattern: XMZ[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XMH
:Pattern: XMH[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
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
:Function: MW
:Pattern: MW[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
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
:Function: MZ
:Pattern: MZ[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
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
:Function: MH
:Pattern: MH[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MT
:Pattern: MT[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
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
:Function: GF
:Pattern: GF[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ,L_Integer:2]
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
