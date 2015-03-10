link=Install["mr"];

(* PDG 2014 data *)

$PDGYear = 2014;
$PDGdata[2014] = { 
		  ASqcdNf5atMZ -> 0.1185, d[ASqcdNf5atMZ] -> 0.0006, 
		  invAEatMZ -> 127.940, d[invAEatMZ] -> 0.014, (* from Review on EW physics *)
		  SW2NDatMZ -> 0.23144, (* from Review on EW physics *) 
		  SW2OS -> 0.22333,(* from Review on EW physics *)
                  MTpole -> 173.21,  d[MTpole] -> Sqrt[ 0.51^2 + 0.71^2], (* GeV *)
		  MBpole ->  4.93482 (* PDG 2-loop quoted 4.78 *), d[MBpole] ->  0.06, (* GeV *)
		  MWpole -> 80.385,  d[MWpole] -> 0.015, (* GeV *)
		  MZpole -> 91.1876, d[MZpole] -> 0.0021, (* GeV *)
		  MHpole -> 125.7,   d[MHpole] -> 0.4,(* GeV *)
		  GF -> 0.000011663787, d[GF] -> 0.000000000006 (* GeV^-2 *) 
		 };

$DegrassiEtAlValues = {
		  ASqcdNf5atMZ -> 0.1184, d[ASqcdNf5atMZ] -> 0.0007, 
                  MTpole -> 173.34,  d[MTpole] -> Sqrt[ 0.3^2 + 0.76^2], (* GeV *)
		  (* MBpole ->  4.78, d[MBpole] ->  0.06, (* GeV *) *)
		  MBpole ->  0, d[MBpole] ->  0, (* GeV *)
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

ReportPDGValues[] := (
		      Print["MW = ", PDG`MW, "\[PlusMinus]",PDG`dMW," GeV"];
		      Print["MZ = ", PDG`MZ, "\[PlusMinus]",PDG`dMZ," GeV"];
		      Print["MH = ", PDG`MH, "\[PlusMinus]",PDG`dMH," GeV"];
		      Print["MT = ", PDG`MT, "\[PlusMinus]",PDG`dMT," GeV"];
		      Print["MB = ", PDG`MB, "\[PlusMinus]",PDG`dMB," GeV"];
		      Print["GF = (", PDG`GF*10^5, "\[PlusMinus]",PDG`dGF*10^5,")*10^-5 GeV^-2"];
		      Print["as(MZ) = ", PDG`asQCD, "\[PlusMinus]",PDG`dasQCD];
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
	   gs^2/(4*Pi)==RunQCDnf6[sc, asp,PDG`MZ,4, PDG`MT] 
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
		(gs^2/(4 Pi) -  RunQCDnf6[scale, PDG`asQCD,PDG`MZ,4, PDG`MT])^2/(0.0001)^2
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
			mmme,
			pars = {g1,g2,gs,yb,yt,lam,m,scale}
	 		 },

DebugPrint["Check QCD vs SM (MT):",    RunQCDnf6[sc, PDG`asQCD,PDG`MZ,4, PDG`MT], " vs ", gs^2/(4 Pi) /. ref];
DebugPrint["Check QCD vs SM (MT*10):", RunQCDnf6[sc*scfactor, PDG`asQCD,PDG`MZ,4, PDG`MT]," vs ", gs^2/(4 Pi) /. ppp];
DebugPrint["Check QCD vs SM (MT/10):", RunQCDnf6[sc/scfactor, PDG`asQCD,PDG`MZ,4, PDG`MT], " vs ", gs^2/(4 Pi) /. mmm];

			refe = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. ref, sc];
			refeAGF = RunParsFromPoleMassesAndAsAndGf[][gs^2/(4 Pi) /. ref, sc];

			pppe = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. ppp, sc scfactor];
			pppeR = RunSM[ pppe, sc]; 

			pppeAGF = RunParsFromPoleMassesAndAsAndGf[][gs^2/(4 Pi) /. ppp, sc scfactor];
			pppeAGFR = RunSM[ pppeAGF, sc]; 

			mmme = RunParsFromPoleMassesAndAsWithMb[][gs^2/(4 Pi) /. mmm, sc / scfactor];
			mmmeR = RunSM[ mmme, sc]; 

			mmmeAGF = RunParsFromPoleMassesAndAsAndGf[][gs^2/(4 Pi) /. mmm, sc / scfactor];
			mmmeAGFR = RunSM[ mmmeAGF, sc]; 
(*
			Print[" check ppp := ", RunSM[ Sequence @@ ({g1,g2,gs,yb,yt,lam,m,scale} /. ppp), sc]];
			Print[" check mmm := ", RunSM[ Sequence @@ ({g1,g2,gs,yb,yt,lam,m,scale} /. mmm), sc]];
*)
			ppps = RunSM[ ppp, sc]; (* sc*scfactor -> sc *)
			mmms = RunSM[ mmm, sc]; (* sc/scfactor -> sc *)
DebugPrint["ref=", ref];
DebugPrint["refe=", refe];
DebugPrint["refeAGF=", refeAGF];

DebugPrint["ppps=", ppps];
DebugPrint["mmms=", mmms];
DebugPrint["pppeR=", pppeR];
DebugPrint["pppeAGFR=", pppeAGFR];
DebugPrint["mmmeR=", mmmeR];
DebugPrint["mmmeAGFR=", mmmeAGFR];

DebugPrint["ppp=", ppp];
DebugPrint["pppe=", pppe];
DebugPrint["pppeAGF=", pppeAGF];

DebugPrint["mmm=", mmm];
DebugPrint["mmme=", mmme];
DebugPrint["mmmeAGF=", mmmeAGF];

			

			(*MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]*)	
	Return[Map[ (# -> { {# /. ref, (# /. ppps) - (# /. ref), (# /. mmms) - (# /. ref)},
			    {# /. refe, (# /. pppeR) - (# /. refe), (# /. mmmeR) - (# /. refe)},	
			    {# /. refeAGF, (# /. pppeAGFR) - (# /. refeAGF), (# /. mmmeAGFR) - (# /. refeAGF)}}	
				)&, pars]];
			] 

(*
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

*)

EstimateTheorUncertainty[ runpars_, scalefactor_:2] := Module[ { 
			(*
			refMZ = MZ[gp, g, gs, yt, lam, mu0, scale,loop],
			refMW = MW[gp, g, gs, yt, lam, mu0, scale,loop],
			refMT = MT[gp, g, gs, yt, lam, mu0, scale,loop],
			refMH = MH[gp, g, gs, yt, lam, mu0, scale,loop],
			refGF = GF[gp, g, gs, yt, lam, mu0, scale,loop],
			*)
			pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars,
			parsP, 
			parsM, 
			ppp, mmm, ref
			},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			parsP = RunSM[runpars,  scalefactor  * (scale /. runpars)];
			parsM = RunSM[runpars, 1/scalefactor * (scale /. runpars)];
			ref = ({
					MZ[#],
					MW[#],
					MT[#],
					MH[#],
					GF[#]} &[runpars]);
			ppp = ({
					MZ[#],
					MW[#],
					MT[#],
					MH[#],
					GF[#]} &[parsP]);
			mmm = ({
					MZ[#],
					MW[#],
					MT[#],
					MH[#],
					GF[#]} &[parsM]);
			Return[MapThread[List[##] &,{ref,ppp/ref-1,mmm/ref-1}]]
			, Print[" Not All parameters specified ", pars, " from ", runpars]]];


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

(* move important parrameters to the beginning *)
RunParsFromPoleMassesAndGfalpha[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			RunParsFromPoleMassesAndAsWithMb[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};
RunParsFromPoleMassesAndGfalphaList[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			RunParsFromPoleMassesAndAsWithMb[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};

RunParsFromPoleMassesAndAsWithMb[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, asmu_, smu_ ] := Module[{aew,seq,aa1,aa2, impliciteq,res,solaew},
DebugPrint["DEBUG: mb = ", mb];
		aa1 = (4 Pi) 3/5 * a1[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
DebugPrint["DEBUG: 1/aewSol = ", 1/aew /. solaew];
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

(*
RunParsFromPoleMassesAndGf[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,PDG`asQCD,PDG`MZ,4,PDG`MT]},
			RunParsFromPoleMassesAndAsAndGf[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]];
*)
RunParsFromPoleMassesAndGf[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			RunParsFromPoleMassesAndAsAndGf[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};
RunParsFromPoleMassesAndGfList[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			Last /@ RunParsFromPoleMassesAndAsAndGf[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};

RunParsFromPoleMassesAndAsAndGf[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, asmu_, smu_ ] := Module[{aew,seq,aa1,aa2, impliciteq,res},
		aew = alphaGF[ mb, mw, mz, mh, mt, smu] /. aQCD[smu] -> asmu/(4 Pi) /. Gf -> gf; 
DebugPrint["DEBUG: 1/aewGF = ", 1/aew];
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
	       	res = res /. {aQCD[smu]-> looptag "asmu" asmu/(4 Pi), Gf->gf, aEW[smu]-> looptag "aew" aew/(4 Pi)};
		(*Print["1/a=", 1/aew /. solaew];
		Print["as=", asmu,":", res, "->",res /. h->1];*)
		If[ NumericQ[looptag], Return[ res /. xx_String -> 1 (* remove string tags *)], 
		(* else *) 
		(* Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. lootag->1)])] *)
		Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. looptag->1 /. xx_String->1)/(b /. looptag->0)])]
		  ];
		];


(* NB: this should be after all definitions, or default values will not change after update ofr PDG values *)
SetPDGValues[]

alpha[gp_,g_] := g^2 gp^2/(g^2 + gp^2)/(4 Pi)
alpha[runpars_List] := (g1^2 g2^2/(g1^2 + g2^2)/(4 Pi)) /. runpars
