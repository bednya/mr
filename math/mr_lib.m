link=Install["mr"];

(* PDG 2014 data *)

$PDGYear = 2014;
$PDGdata[2014] = { 
		  ASqcdNf5atMZ -> 0.1185, d[ASqcdNf5atMZ] -> 0.0006, 
		  invAEatMZ -> 127.940, d[invAEatMZ] -> 0.014, (* from Review on EW physics *)
		  SW2NDatMZ -> 0.23144, (* from Review on EW physics *) 
		  SW2OS -> 0.22333,(* from Review on EW physics *)
                  MTpole -> 173.21,  d[MTpole] -> Sqrt[ 0.51^2 + 0.71^2], (* GeV *)
		  MBpole ->  4.78, d[MBpole] ->  0.06, (* GeV *)
		  MWpole -> 80.385,  d[MWpole] -> 0.015, (* GeV *)
		  MZpole -> 91.1876, d[MZpole] -> 0.0021, (* GeV *)
		  MHpole -> 125.7,   d[MHpole] -> 0.4,(* GeV *)
		  GF -> 0.000011663787, d[GF] -> 0.000000000006 (* GeV^-2 *) 
		 };

$DegrassiEtAlValues = {
		  ASqcdNf5atMZ -> 0.1184, d[ASqcdNf5atMZ] -> 0.0007, 
                  MTpole -> 173.34,  d[MTpole] -> Sqrt[ 0.3^2 + 0.76^2], (* GeV *)
		  MBpole ->  4.78, d[MBpole] ->  0.06, (* GeV *)
		  MWpole -> 80.384,  d[MWpole] -> 0.014, (* GeV *)
		  MZpole -> 91.1876, d[MZpole] -> 0.0021, (* GeV *)
		  MHpole -> 125.15,   d[MHpole] -> 0.24,(* GeV *)
		  GF -> 0.000011663787, d[GF] -> 0.000000000006 (* GeV^-2 *) 
			};

$SMdata = $PDGdata[$PDGYear];


SetPDGValues[smdata_:$PDGdata[$PDGYear]] := (
			PDG`MW  =   MWpole  /. smdata /. $PDGdata[$PDGYear];
			PDG`dMW = d[MWpole] /. smdata /. $PDGdata[$PDGYear];
			PDG`MZ  =   MZpole  /. smdata /. $PDGdata[$PDGYear];
			PDG`dMZ = d[MZpole] /. smdata /. $PDGdata[$PDGYear];
			PDG`MH  =   MHpole  /. smdata /. $PDGdata[$PDGYear];
			PDG`dMH = d[MHpole] /. smdata /. $PDGdata[$PDGYear];
			PDG`MT  =   MTpole  /. smdata /. $PDGdata[$PDGYear];
			PDG`dMT = d[MTpole] /. smdata /. $PDGdata[$PDGYear];
			PDG`MB  =   MBpole  /. smdata /. $PDGdata[$PDGYear];
			PDG`dMB = d[MBpole] /. smdata /. $PDGdata[$PDGYear];
			PDG`GF  =   GF      /. smdata /. $PDGdata[$PDGYear];
			PDG`dGF = d[GF]     /. smdata /. $PDGdata[$PDGYear];
		        PDG`asQCD  =   ASqcdNf5atMZ  /. smdata /. $PDGdata[$PDGYear];
		        PDG`dasQCD = d[ASqcdNf5atMZ] /. smdata /. $PDGdata[$PDGYear];
		  	)
		

(* Default "PDG" values *)
SetPDGValues[]

ReportPDGValues[] := (
		      Print["MW = ", PDG`MW, "+/-",PDG`dMW," GeV"];
		      Print["MZ = ", PDG`MZ, "+/-",PDG`dMZ," GeV"];
		      Print["MH = ", PDG`MH, "+/-",PDG`dMH," GeV"];
		      Print["MT = ", PDG`MT, "+/-",PDG`dMT," GeV"];
		      Print["MB = ", PDG`MB, "+/-",PDG`dMB," GeV"];
		      Print["GF = (", PDG`GF*1^5, "+/-",PDG`dGF*10^5,")*10^-5 GeV^-2"];
		      Print["as(MZ) = ", PDG`asQCD, "+/-",PDG`dasQCD];
		     );


(* 3 sigma range [x - 3 s, x + 3 s] divided into n subintervals *)

ThreeSigmaRange[x_,s_,n_:2] := Module[{step = 6 s /n}, Table[ x - 3 s + step * (i-1), {i,1,n+1}]] 


(* find running parameters given (pseudo)observables, and a renormalization scale *)

sol[sc_][mzp_:PDG`MZ,mwp_:PDG`MW,mtp_:PDG`MT,mhp_:PDG`MH,gfp_:PDG`GF,asp_:PDG`asQCD] := Join[FindRoot[ {
	   MZ[g1, g2, gs, 0, yt, lam, m, sc] == mzp (* 91.1876 *),
	   MW[g1, g2, gs, 0, yt, lam, m, sc] == mwp (* 80.385 *),
	   MT[g1, g2, gs, 0, yt, lam, m, sc] == mtp (* sc *),
	   MH[g1, g2, gs, 0, yt, lam, m, sc] == mhp (* 125.7 *),
	   GF[g1, g2, gs, 0, yt, lam, m, sc] == gfp (* 0.000011663787 *),
	   gs^2/(4*Pi)==RunQCD[sc, asp,PDG`MZ,4, PDG`MT] 
	  },{
		{g1, 0.357561},
		{g2 , 0.64822},
		{gs , 1.1666},
		{yt,0.93558},
		{lam,0.12711},
		{m,132.03}
	  }],{yb -> 0,scale -> sc}];

(* evaluate chi2 function for given set of running parameters at the given scale *)

chiSquare[sol_]:=(   
		(MZ[g1, g2, gs, 0, yt, lam, m, scale] - PDG`MZ)^2/(PDG`dMZ)^2 + 
	 	(MW[g1, g2, gs, 0, yt, lam, m, scale] - PDG`MW)^2/(PDG`dMW)^2 + 
		(MT[g1, g2, gs, 0, yt, lam, m, scale] - PDG`MT)^2/(PDG`dMT)^2 + 
		(MH[g1, g2, gs, 0, yt, lam, m, scale] - PDG`MH)^2/(PDG`dMH)^2 + 
		(GF[g1, g2, gs, 0, yt, lam, m, scale] - PDG`GF)^2/(PDG`dGF)^2 (* check this *)	 + 
		(gs^2/(4 Pi) -  RunQCD[scale, PDG`asQCD,PDG`MZ,4, PDG`MT])^2/(0.0001)^2
			) /. sol ;

(* theoretical error for observables *)
(*
EstimateTheorUncertaintyInMatchingDegrassi[scale_, scalefactor_:10] := Module[ {
	  		ref = sol[scale][], 
			refe,
	  		ppp = sol[scale scalefactor][],
			pppe,
	  		mmm = sol[scale/scalefactor][],
			mmme
	 		 },
			refe = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. ref, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
			pppe = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. ppp, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
			mmme = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. mmm, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
Print["refe=", refe];
Print["pppe=", pppe];
Print["mmme=", mmme];
			ppp = RunSM[ Sequence @@ (Last /@ ppp), scale]; (* scale*scalefactor -> scale *)
			mmm = RunSM[ Sequence @@ (Last /@ mmm), scale]; (* scale/scalefactor -> scale *)
			MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]	

			] 
*)
EstimateTheorUncertaintyInMatchingDegrassi[sc_, scfactor_:10] := Module[ {
	  		ref = sol[sc][], 
			refe,
	  		ppp = sol[sc scfactor][],
			pppe,
	  		mmm = sol[sc/scfactor][],
			mmme
	 		 },
			refe = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. ref, sc];
			pppe = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. ppp, sc scfactor];
			pppeR = RunSM[ pppe, sc]; 
			mmme = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. mmm, sc / scfactor];
			mmmeR = RunSM[ mmme, sc]; 
(*
			Print[" check ppp := ", RunSM[ Sequence @@ ({g1,g2,gs,yb,yt,lam,m,scale} /. ppp), sc]];
			Print[" check mmm := ", RunSM[ Sequence @@ ({g1,g2,gs,yb,yt,lam,m,scale} /. mmm), sc]];
*)
			ppp = RunSM[ ppp, sc]; (* sc*scfactor -> sc *)
			mmm = RunSM[ mmm, sc]; (* sc/scfactor -> sc *)
Print["ref=", ref];
Print["refe=", refe];

Print["ppp=", ppp];
Print["pppe=", pppe];
Print["pppeR=", pppeR];

Print["mmm=", mmm];
Print["mmme=", mmme];
Print["mmmeR=", mmmeR];

			(*MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]*)	
	Return[0];
			] 

EstimateTheorUncertainty[ gp_, g_, gs_, yt_, lam_, mu0_,scale_, loop_:2, scalefactor_:2] := Module[ { 
			(*
			refMZ = MZ[gp, g, gs, yt, lam, mu0, scale,loop],
			refMW = MW[gp, g, gs, yt, lam, mu0, scale,loop],
			refMT = MT[gp, g, gs, yt, lam, mu0, scale,loop],
			refMH = MH[gp, g, gs, yt, lam, mu0, scale,loop],
			refGF = GF[gp, g, gs, yt, lam, mu0, scale,loop],
			*)
			parsP = Flatten[{RunSM[gp, g, gs, yt, lam, mu0, scale, scalefactor  * scale],  loop}],
			parsM = Flatten[{RunSM[gp, g, gs, yt, lam, mu0, scale, 1/scalefactor * scale],  loop}],
			ppp, mmm, ref
			},
			Print[ parsP];
			ref = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[{gp,g,gs,yt,lam,mu0,scale,loop}]);
			ppp = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[parsP]);
			mmm = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[parsM]);
			Return[MapThread[List[##] &,{ref,ppp/ref-1,mmm/ref-1}]]
			];


(*
ScaleDependence[O_][runpars_,loop_:2] := O[ Sequence @@ (RunSM[ Sequence @@ Evaluate[{gp, g, gs, yt, lam, mu0, mu} /. runpars], #]),loop] & 
*)

ScaleDependence[O_][runpars_] := O[ Sequence @@ (RunSM[ Sequence @@ Evaluate[{g1, g2, gs, yb,yt, lam, m, scale} /. runpars], #])] & 


(*
RunParsFromPoleMassesAndAs[mz_:PDG`MZ,mw_:PDG`MW,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF,asmu_,smu_] := Module[{aew,seq,aa1,aa2, impliciteq},
		aa1 = (4 Pi) 3/5 * a1[0,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[0,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
		(* Couplings *)
		{ gp->Sqrt[ (4 Pi)^2 3/5 a1[0, mw,mz,mh,mt,smu]], 
		  g -> Sqrt[ (4 Pi)^2 a2[0, mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yt -> Sqrt[ (4 Pi)^2 at[0, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[0, mw,mz,mh,mt,smu], 
		  mu -> smu} /. {aQCD[smu]->asmu/(4 Pi), Gf->gf, aEW[smu]->aew/(4 Pi)} /. solaew
		  ];
*)
(* NB: one needs to supply running as! *)

RunParsFromPoleMassesAndAsWithMb[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, asmu_, smu_ ] := Module[{aew,seq,aa1,aa2, impliciteq,res,solaew},
		aa1 = (4 Pi) 3/5 * a1[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
		(* Couplings *)
		res = { g1->Sqrt[ (4 Pi)^2 3/5 a1[mb, mw,mz,mh,mt,smu]], 
		  g2 -> Sqrt[ (4 Pi)^2 a2[mb,mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yb -> Sqrt[ (4 Pi)^2 ab[mb, mw,mz,mh,mt,smu]],
		  yt -> Sqrt[ (4 Pi)^2 at[mb, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[mb,mw,mz,mh,mt,smu], 
		(* ms(mu) = mh(mu) ? *)
		    m -> mh * Sqrt[mmHMMH[mb,mw,mz,mh,mt,smu]],
		    vev -> vev[mb,mw,mz,mh,mt,smu],
		    scale -> smu};
		(*		Print[res, ":::",2^(-1/2)*Gf/(4*Pi)^2*mh^2*(1+aEW[mu]*yH[1,0]+aEW[mu]*aQCD[mu]*yH[1,1]+aEW[mu]^2*yH[2,0]) /. mu->smu/. XH[mb,mw,mz,mh,mt,smu,1,2], "::",alam[mb, mw,mz,mh,mt,smu], ":::", XH[mb,mw,mz,mh,mt,smu,2,1]];*)
	       	res = res /. {aQCD[smu]-> looptag "asmu" asmu/(4 Pi), Gf->gf, aEW[smu]-> looptag "aew" aew/(4 Pi)} /. solaew;
		(*Print["1/a=", 1/aew /. solaew];
		Print["as=", asmu,":", res, "->",res /. h->1];*)
		If[ NumericQ[looptag], Return[ res /. a_String -> 1 (* remove string tags *)], 
		(* else *) 
		(* Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. lootag->1)])] *)
		Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. lootag->1 /. a_String->1)/(b /. looptag->0)])]
		  ];
		];


alpha[gp_,g_] := g^2 gp^2/(g^2 + gp^2)/(4 Pi)
alpha[runpars_List] := (g1^2 g2^2/(g1^2 + g2^2)/(4 Pi)) /. runpars

