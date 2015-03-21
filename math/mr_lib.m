
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


error[x_][0.]=0;
error[x_][0.,0.]=0;
error[x_][0]=0;
error[x_][0,0]=0;

(*
error1[0]=0; (* error due to different procedures to extract alpha from Gf - solution of implictit equation, or explicit reexpansion *)
error1[0.]=0; 

error10[0]=0; (* error due to different procedures to extract running parameters *)
error10[0.]=0; 

error12[0., 0.]=0; (* error due to matching scale variation, running parameters are obtained via solving pole masses and Gf wrt MS pars *)

error2[0.,0.]=0; (* error due to matching scale variation, alpha is extracted from gf via numerical solution of implicit equation *)
error3[0.,0.]=0; (* error due to matching scale variation, alpha is extracted from gf via explict formula *) 
*)

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

sol[sc_,qcdth_:PDG`MT][mzp_:PDG`MZ,mwp_:PDG`MW,mtp_:PDG`MT,mhp_:PDG`MH,gfp_:PDG`GF,asp_:PDG`asQCD] := Join[FindRoot[ {
	   MZ[g1, g2, gs, 0, yt, lam, m, sc] == mzp (* 91.1876 *),
	   MW[g1, g2, gs, 0, yt, lam, m, sc] == mwp (* 80.385 *),
	   MT[g1, g2, gs, 0, yt, lam, m, sc] == mtp (* sc *),
	   MH[g1, g2, gs, 0, yt, lam, m, sc] == mhp (* 125.7 *),
	   GF[g1, g2, gs, 0, yt, lam, m, sc] == gfp (* 0.000011663787 *),
	   gs^2/(4*Pi)==RunQCDnf6[sc, asp, mzp,4, mzp, qcdth]  (* threshold *)
	  },{
		{g1, 0.357561},
		{g2 , 0.64822},
		{gs , 1.1666},
		{yt,0.93558},
		{lam,0.12711},
		{m,132.03}
	  }],{yb -> 0,scale -> sc}];

sol2loopMT[sc_, qcdth_:PDG`MT][mzp_:PDG`MZ,mwp_:PDG`MW,mtp_:PDG`MT,mhp_:PDG`MH,gfp_:PDG`GF,asp_:PDG`asQCD] := Join[FindRoot[ {
	   MZ[g1, g2, gs, 0, yt, lam, m, sc] == mzp (* 91.1876 *),
	   MW[g1, g2, gs, 0, yt, lam, m, sc] == mwp (* 80.385 *),
	   MT[g1, g2, gs, 0, yt, lam, m, sc, 1, ExcludeQCD3loop -> True] == mtp (* sc *),
	   MH[g1, g2, gs, 0, yt, lam, m, sc] == mhp (* 125.7 *),
	   GF[g1, g2, gs, 0, yt, lam, m, sc] == gfp (* 0.000011663787 *),
	   gs^2/(4*Pi)==RunQCDnf6[sc, asp,mzp, 4, mtp, qcdth] 
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


Options[ EstimateTheorUncertaintyInMatchingDegrassi] = {"UseImplicitExtraction" -> False}

EstimateTheorUncertaintyInMatchingDegrassi[sc_?NumberQ, scfactor_Integer:10, OptionsPattern[]] := Module[ {
	  		ref, 
			refe,
	  		ppp,
			pppe,
	  		mmm,
			mmme,
			asLmatching = RunQCDnf6[sc,PDG`asQCD,PDG`MZ,4,PDG`MT,PDG`MT/scfactor], (* alpha_s from nf=5 QCD to nf=6 QCD matching at the scale Mt/2 *)
			asCmatching = RunQCDnf6[sc,PDG`asQCD,PDG`MZ,4,PDG`MT,PDG`MT], (* alpha_s from nf=5 QCD to nf=6 QCD matching at the scale Mt/2 *)
			asHmatching = RunQCDnf6[sc,PDG`asQCD,PDG`MZ,4,PDG`MT,scfactor*PDG`MT], (* alpha_s from nf=5 QCD to nf=6 QCD matching at the scale 2*Mt *)
			pars = {g1,g2,gs,yb,yt,lam,m,scale},
			DoImplicit = TrueQ[ OptionValue[ "UseImplicitExtraction"]]
	 		 },

If[ DoImplicit,
	  		ref = sol[sc][]; 
DebugPrint["ref=", ref];
Print["IMPLICIT!!!"];
 ];

Print["SCALE FACTOR: ", scfactor];
			(*Print["AS MATCING UNCERTAINTY at ", sc," GeV:", asCmatching, 100*(1 - asLmatching/asCmatching), "|", 100*(1-asHmatching/asCmatching)]; *)

		        refAsH = RunParsImplAlpha[][asHmatching,sc];	
		        refAsL = RunParsImplAlpha[][asLmatching,sc];	


			refe = RunParsFromGfAndPoleMasses[sc,"AlphaFromGfImplicit"->True];
			refeAGF = RunParsFromGfAndPoleMasses[sc,"AlphaFromGfImplicit"->False];

			pppe = RunParsFromGfAndPoleMasses[sc*scfactor,"AlphaFromGfImplicit"->True];
			pppeR = RunSM[ pppe, sc]; 

			pppeAGF = RunParsFromGfAndPoleMasses[sc*scfactor,"AlphaFromGfImplicit"->False];
			pppeAGFR = RunSM[ pppeAGF, sc]; 

			mmme = RunParsFromGfAndPoleMasses[sc/scfactor,"AlphaFromGfImplicit"->True];
			mmmeR = RunSM[ mmme, sc]; 

			mmmeAGF = RunParsFromGfAndPoleMasses[sc/scfactor,"AlphaFromGfImplicit"->False];
			mmmeAGFR = RunSM[ mmmeAGF, sc]; 

(*DebugPrint = Print;*)
DebugPrint["refe=", refe];
DebugPrint["refeAGF=", refeAGF];

DebugPrint["pppeR=", pppeR];
DebugPrint["pppeAGFR=", pppeAGFR];
DebugPrint["mmmeR=", mmmeR];
DebugPrint["mmmeAGFR=", mmmeAGFR];

DebugPrint["pppe=", pppe];
DebugPrint["pppeAGF=", pppeAGF];

DebugPrint["mmme=", mmme];
DebugPrint["mmmeAGF=", mmmeAGF];
Clear[DebugPrint];
			

			(*MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]*)	
If[ DoImplicit,
	Return[Map[ (# -> (* leave only errors  (# /. refe)  *)
				       + error[1][ (# /. refe) - (# /. refeAGF)] 
				       + error[2][(# /. pppeR) - (# /. refe), (# /. mmmeR) - (# /. refe)] 
				       + error[3][(# /. pppeAGFR) - (# /. refeAGF), (# /. mmmeAGFR) - (# /. refeAGF)]
				       + error[5][ (# /. refAsH) - (# /. refe), (# /. refAsL) - (# /. refe)] 
				       + error[10][ (# /. refe) - (# /. ref)] 
				)&, pars]],
(* else *)
	Return[Map[ (# ->  (* leave only errors *) 0 * (# /. refe) 
				       + error[1][ (# /. refe) - (# /. refeAGF)]  
				       + error[2][(# /. pppeR) - (# /. refe), (# /. mmmeR) - (# /. refe)] 
				       + error[5][ (# /. refAsH) - (# /. refe), (# /. refAsL) - (# /. refe)] 
				       + error[3][(# /. pppeAGFR) - (# /. refeAGF), (# /. mmmeAGFR) - (# /. refeAGF)]
				)&, pars]];
			] 

	];

(*
EstimateTheorUncertaintyInMatchingDegrassi[sc_?NumberQ, scfactor_Integer:10, OptionsPattern[]] := Module[ {
	  		ref, 
			refe,
	  		ppp,
			pppe,
	  		mmm,
			mmme,
			pars = {g1,g2,gs,yb,yt,lam,m,scale},
			asref = RunQCDnf6[sc, PDG`asQCD,PDG`MZ,4, PDG`MT],
			asppp = RunQCDnf6[sc*scfactor, PDG`asQCD,PDG`MZ,4, PDG`MT],
			asmmm = RunQCDnf6[sc/scfactor, PDG`asQCD,PDG`MZ,4, PDG`MT],
			DoImplicit = TrueQ[ OptionValue[ "UseImplicitExtraction"]]
	 		 },

If[ DoImplicit,
	  		ref = sol[sc][]; 
DebugPrint["ref=", ref];
Print["IMPLICIT!!!"];
 ];

			(*	
			Print["ASREF:", asref]; asref = gs^2/(4 Pi) /. sol[sc][];	Print["ASREF:", asref]; 
			Print["ASPPP:", asppp]; asref = gs^2/(4 Pi) /. sol[sc*scfactor][];	Print["ASPPP:", asref];	
			Print["ASMMM:", asmmm]; asref = gs^2/(4 Pi) /. sol[sc/scfactor][];	Print["ASMMM:", asref];		
			*)

			refe = RunParsFromPoleMassesAndAsWithMb[][asref, sc];
			refeAGF = RunParsFromPoleMassesAndAsAndGf[][asref, sc];

			pppe = RunParsFromPoleMassesAndAsWithMb[][asppp, sc scfactor];
			pppeR = RunSM[ pppe, sc]; 

			pppeAGF = RunParsFromPoleMassesAndAsAndGf[][asppp, sc scfactor];
			pppeAGFR = RunSM[ pppeAGF, sc]; 

			mmme = RunParsFromPoleMassesAndAsWithMb[][asmmm, sc / scfactor];
			mmmeR = RunSM[ mmme, sc]; 

			mmmeAGF = RunParsFromPoleMassesAndAsAndGf[][asmmm, sc / scfactor];
			mmmeAGFR = RunSM[ mmmeAGF, sc]; 

DebugPrint["refe=", refe];
DebugPrint["refeAGF=", refeAGF];

DebugPrint["pppeR=", pppeR];
DebugPrint["pppeAGFR=", pppeAGFR];
DebugPrint["mmmeR=", mmmeR];
DebugPrint["mmmeAGFR=", mmmeAGFR];

DebugPrint["pppe=", pppe];
DebugPrint["pppeAGF=", pppeAGF];

DebugPrint["mmme=", mmme];
DebugPrint["mmmeAGF=", mmmeAGF];
			

			(*MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]*)	
If[ DoImplicit,
	Return[Map[ (# -> (* leave only errors  (# /. refe)  *)
				       + error[1][ (# /. refe) - (# /. refeAGF)] 
				       + error[2][(# /. pppeR) - (# /. refe), (# /. mmmeR) - (# /. refe)] 
				       + error[3][(# /. pppeAGFR) - (# /. refeAGF), (# /. mmmeAGFR) - (# /. refeAGF)]
				       + error[10][ (# /. refe) - (# /. ref)] 
				)&, pars]],
(* else *)
	Return[Map[ (# ->  (* leave only errors *) 0 * (# /. refe) 
				       + error[1][ (# /. refe) - (# /. refeAGF)]  
				       + error[2][(# /. pppeR) - (# /. refe), (# /. mmmeR) - (# /. refe)] 
				       + error[3][(# /. pppeAGFR) - (# /. refeAGF), (# /. mmmeAGFR) - (# /. refeAGF)]
				)&, pars]];
			] 

	];

*)
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





ScaleDependence[O_][runpars_] := O[ Sequence @@ (RunSM[ Sequence @@ Evaluate[{g1, g2, gs, yb,yt, lam, m, scale} /. runpars], #])] & 


(*
(*RunParsFromPoleMassesAndAs[mz_:PDG`MZ,mw_:PDG`MW,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF,asmu_,smu_] := Module[{aew,seq,aa1,aa2, impliciteq},*)
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

Options[ RunParsFromGfAndPoleMasses ] = {
						"LoopTag" ->1,
						"ReturnList" -> False,
						"AlphaFromGfImplicit" -> True 
					};


RunParsFromGfAndPoleMasses[mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_?NumericQ, OptionsPattern[]]  := Module[ { 
					looptag = OptionValue["LoopTag"]},
(*Print["mt = ", mt];*)
(*Print["mh = ", mh];*)
(*Print["mw = ", mw];*)
(*Print["mb = ", mb];*)
(*Print["mz = ", mz];*)
					If[ TrueQ[ OptionValue["AlphaFromGfImplicit"] ],
						If [ TrueQ [ OptionValue [ "ReturnList"] ], Return[ RunParsImplicitAlphaList[looptag][ mt, mh, asMZ, mw, mb, mz, gf, smu] ],
						(* else *)				      Return[ RunParsImplicitAlpha[looptag][ mt, mh, asMZ, mw, mb, mz, gf, smu]]],
					(* else *)
						If [ TrueQ [ OptionValue [ "ReturnList"] ], Return[ RunParsExplicitAlphaList[looptag][ mt, mh, asMZ, mw, mb, mz, gf, smu]],
						(* else *)				      Return[ RunParsExplicitAlpha[looptag][ mt, mh, asMZ, mw, mb, mz, gf, smu]]]
					  ];
];
				 


RunParsImplicitAlpha[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			RunParsImplAlpha[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};
RunParsImplicitAlphaList[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			Last /@ RunParsImplAlpha[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};


(* 
	alpha_ew is found via numerical solution of implicit equation relating Gf and running alpha_ew (given alpha_s)
*)

RunParsImplAlpha[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, asmu_?NumericQ, smu_?NumericQ] := Module[{aew,seq,aa1,aa2, impliciteq,res,solaew,mmW = mw^2, mmZ=mz^2, dyZ, dyW},
(* old procedure - left for debug *)
		aa1 = (4 Pi) 3/5 * a1[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
Print["DEBUG (old) : 1/aewSol = ", 1/aew /. solaew];
(* new procedure *)
		dyZ = 1 + looptag*"aew"*aEW[smu]*yZ[1,0]
			+ looptag^2 ("aew"*"as"*aEW[smu]*aQCD[smu]*yZ[1,1]+ "aew"^2 aEW[smu]^2*yZ[2,0])/.XZ[mb,mw,mz,mh,mt,smu] /. x_String->1 /. aEW[smu]->aew/(4 Pi) /. aQCD[smu]-> asmu/(4 Pi);
		dyW = 1 + looptag*"aew"*aEW[smu]*yW[1,0]
			+ looptag^2 ("aew"*"as"*aEW[smu]*aQCD[smu]*yW[1,1]+ "aew"^2 aEW[smu]^2*yW[2,0])/.XW[mb,mw,mz,mh,mt,smu] /. x_String->1 /. aEW[smu]->aew/(4 Pi) /. aQCD[smu]-> asmu/(4 Pi) ;

		solaew = FindRoot[Evaluate[gf] == aew*Pi/Sqrt[2]/mmW/dyW/(1 - mmW/mmZ * dyW/dyZ), {aew,1/126.94}];
Print["DEBUG (new) : 1/aewSol = ", 1/aew /. solaew];

		(* Couplings *)
		res = { g1->Sqrt[ (4 Pi)^2 3/5 a1[mb, mw,mz,mh,mt,smu]], 
		  g2 -> Sqrt[ (4 Pi)^2 a2[mb,mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yb -> Sqrt[ (4 Pi)^2 ab[mb, mw,mz,mh,mt,smu]],
		  yt -> Sqrt[ (4 Pi)^2 at[mb, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[mb,mw,mz,mh,mt,smu], 
		(* ms(mu) = mh(mu) ? *)
		    mcheck -> mh * Sqrt[mmHMMH[mb,mw,mz,mh,mt,smu]],
	       m -> (4 Pi) Sqrt[2 alam[mb,mw,mz,mh,mt,smu]] * vev[mb, mw,mz,mh,mt,smu],
		    vev -> vev[mb,mw,mz,mh,mt,smu],
		    scale -> smu};
		(*		Print[res, ":::",2^(-1/2)*Gf/(4*Pi)^2*mh^2*(1+aEW[mu]*yH[1,0]+aEW[mu]*aQCD[mu]*yH[1,1]+aEW[mu]^2*yH[2,0]) /. mu->smu/. XH[mb,mw,mz,mh,mt,smu,1,2], "::",alam[mb, mw,mz,mh,mt,smu], ":::", XH[mb,mw,mz,mh,mt,smu,2,1]];*)
	       	res = res /. {aQCD[smu]-> looptag "asmu" asmu/(4 Pi), Gf->gf, aEW[smu]-> looptag "aew" aew/(4 Pi)} /. solaew;
		(*Print["1/a=", 1/aew /. solaew];
		DebugPrint["as=", asmu,":", res, "->",res /. h->1];*)
		If[ NumericQ[looptag], Return[ res /. a_String -> 1 (* remove string tags *)], 
		(* else *) 
		(* Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. lootag->1)])] *)
		Return[ res /. (a_->b_):>(a->(b/. looptag->0)*corr[Collect[1/(b/.looptag->0)*(Normal[Series[b,{looptag,0,2}]]),looptag,Expand], (b /. lootag->1 /. a_String->1)/(b /. looptag->0)])]
		  ];
		];


(* 
	running parameters are extracted from observables, alpha_ew is found from Gf via explicit expression
*)

RunParsExplicitAlpha[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			RunParsExplAlpha[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};
RunParsExplicitAlphaList[looptag_:1][mt_:PDG`MT,mh_:PDG`MH,asMZ_:PDG`asQCD,mw_:PDG`MW,mb_:PDG`MB,mz_:PDG`MZ,gf_:PDG`GF, smu_ ] := Module[{asmu = RunQCDnf6[smu,asMZ,mz,4,mt]},
			Last /@ RunParsExplAlpha[looptag][mb,mw,mz,mh,mt,gf,asmu,smu]] /; And @@ NumericQ /@ {mt,mh,asMZ,mw,mb,mz,gf,smu};

RunParsExplAlpha[looptag_:1][mb_:PDG`MB,mw_:PDG`MW,mz_:PDG`MZ,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF, asmu_, smu_ ] := Module[{aew,seq,aa1,aa2, impliciteq,res, aewtree = Sqrt[2] gf/Pi * mw^2 * (1 - mw^2/mz^2), A10, cw2 = mw^2/mz^2, aEWtree},
(* old procedure, buggy? *)
		aew = alphaGF[ mb, mw, mz, mh, mt, smu] /. aQCD[smu] -> asmu/(4 Pi) /. Gf -> gf; 

Print["DEBUG: 1/aewtree = ", 1/aewtree];
Print["DEBUG: 1/aewGF = ", 1/aew];

(* new routine, explicit *)
{dZ1, dZ2, dW1, dW2}  =   aEWtree*{ yZ[1,0], 
			        aQCD[smu]*yZ[1,1] + aEWtree*yZ[2,0],
				          yW[1,0], 
			        aQCD[smu]*yW[1,1] + aEWtree*yW[2,0]} /. XW[mb,mw,mz,mh,mt,smu] /. XZ[mb,mw,mz,mh,mt,smu] ;
DebugPrint["DEBUG:", {dZ1, dZ2, dW1, dW2}];
		A10 = (dW1 * (1 - 2 * cw2) + cw2 * dZ1)/(1 - cw2);
		aew = aewtree * ( 1 + A10 
				    + ( (dW2 * (1 - 2 * cw2) + cw2 * ( dZ2 + 2 dW1 dW2 - dW1^2 - dW2^2))/(1 - cw2)
				    + A10^2 (* reexpansion *)
				      )
				) /. aQCD[smu]-> asmu/(4 Pi) /. aEWtree -> aewtree/(4 Pi);
Print["DEBUG: 1/aewGF = ", 1/aew] ;

		(* Couplings *)
		res = { g1->Sqrt[ (4 Pi)^2 3/5 a1[mb, mw,mz,mh,mt,smu]], 
		  g2 -> Sqrt[ (4 Pi)^2 a2[mb,mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yb -> Sqrt[ (4 Pi)^2 ab[mb, mw,mz,mh,mt,smu]],
		  yt -> Sqrt[ (4 Pi)^2 at[mb, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[mb,mw,mz,mh,mt,smu], 
		(* ms(mu) = mh(mu) ? *)
		    mcheck -> mh * Sqrt[mmHMMH[mb,mw,mz,mh,mt,smu]],
	       m -> (4 Pi) Sqrt[2 alam[mb,mw,mz,mh,mt,smu]] * vev[mb, mw,mz,mh,mt,smu],
		    vev -> vev[mb,mw,mz,mh,mt,smu],
		    scale -> smu};
		(*		Print[res, ":::",2^(-1/2)*Gf/(4*Pi)^2*mh^2*(1+aEW[mu]*yH[1,0]+aEW[mu]*aQCD[mu]*yH[1,1]+aEW[mu]^2*yH[2,0]) /. mu->smu/. XH[mb,mw,mz,mh,mt,smu,1,2], "::",alam[mb, mw,mz,mh,mt,smu], ":::", XH[mb,mw,mz,mh,mt,smu,2,1]];*)
	       	res = res /. {aQCD[smu]-> looptag "asmu" asmu/(4 Pi), Gf->gf, aEW[smu]-> looptag "aew" aew/(4 Pi)};
		(*Print["1/a=", 1/aew /. solaew];
		DebugPrint["as=", asmu,":", res, "->",res /. h->1];*)
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


