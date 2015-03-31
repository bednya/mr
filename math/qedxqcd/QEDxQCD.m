BeginPackage["QEDxQCD`"]


RunQEDxQCD::usage = "RunQEDxQCD[as0, ae0, mu0, muf, L, options] - Return {as(muf), ae(muf)} calculated with L-loop evolution from as(mu0) = as0, ae(mu0) = ae0, options are \"nL\"->3, \"nU\"->2, and \"nD\"->3 specify number of active charged leptons, up-type quarks and down-type quarks (default is nf=5 QCDxQED)"

Begin["Private`"]


Charges = {Qu -> 2/3, Qd -> -1/3, Qe -> -1}
Casimirs = {cF -> 4/3, cA -> 3 , Tf -> 1/2, dF -> 3}

AsBetaCoefficients = {
		      {2, 1} -> 8 cA Nd Qd^2 Tf - 4 cF Nd Qd^2 Tf + 8 cA Nu Qu^2 Tf - 4 cF Nu Qu^2 Tf, 
		      {1, 2} -> -2 Nd Qd^4 Tf - 44/9 dF Nd^2 Qd^4 Tf - 44/9 Nd Nl Qd^2 Qe^2 Tf - 88/9 dF Nd Nu Qd^2 Qu^2 Tf - 44/9 Nl Nu Qe^2 Qu^2 Tf - 2 Nu Qu^4 Tf - 44/9 dF Nu^2 Qu^4 Tf, 
		      {1, 1} -> 4 Nd Qd^2 Tf + 4 Nu Qu^2 Tf, 
(* pure QCD, 4-loop from RunDec *)
		      {4, 0} -> -(149753/6) - 1093/729 (Nd + Nu)^3 + (Nd + Nu)^2 (-(50065/162) - (6472 Zeta[3])/ 81) - 3564 Zeta[3] + (Nd + Nu) (1078361/162 + (6508 Zeta[3])/27),
		      {3, 0} -> -((2857 cA^3)/54) + 1415/27 cA^2 Nd Tf + 205/9 cA cF Nd Tf - 2 cF^2 Nd Tf + 1415/27 cA^2 Nu Tf + 205/9 cA cF Nu Tf - 2 cF^2 Nu Tf - 158/27 cA Nd^2 Tf^2 - 44/9 cF Nd^2 Tf^2 - 316/27 cA Nd Nu Tf^2 - 88/9 cF Nd Nu Tf^2 - 158/27 cA Nu^2 Tf^2 - 44/9 cF Nu^2 Tf^2,  
		      {2, 0} -> -((34 cA^2)/3) + (20 cA Nd Tf)/3 + 4 cF Nd Tf + (20 cA Nu Tf)/3 + 4 cF Nu Tf, 
                      {1, 0} -> -((11 cA)/3) + ( 4 Nd Tf)/3 + (4 Nu Tf)/3
		      } /. Casimirs /. Charges;


AlBetaCoefficients = { 
		{2, 1} -> 133/9 cA cF dF Nd Qd^2 - 2 cF^2 dF Nd Qd^2 + 133/9 cA cF dF Nu Qu^2 - 2 cF^2 dF Nu Qu^2 - 44/9 cF dF Nd^2 Qd^2 Tf - 44/9 cF dF Nd Nu Qd^2 Tf - 44/9 cF dF Nd Nu Qu^2 Tf - 44/9 cF dF Nu^2 Qu^2 Tf, 
	        {1, 2} -> -4 cF dF Nd Qd^4 - 4 cF dF Nu Qu^4, 
		{1, 1} -> 4 cF dF Nd Qd^2 + 4 cF dF Nu Qu^2, 
(* pure qed *)
		{0, 3} -> -2 dF Nd Qd^6 - 44/9 dF^2 Nd^2 Qd^6 - 44/9 dF Nd Nl Qd^4 Qe^2 - 44/9 dF Nd Nl Qd^2 Qe^4 - 2 Nl Qe^6 - (44 Nl^2 Qe^6)/9 - 44/9 dF^2 Nd Nu Qd^4 Qu^2 - 44/9 dF Nl Nu Qe^4 Qu^2 - 44/9 dF^2 Nd Nu Qd^2 Qu^4 - 44/9 dF Nl Nu Qe^2 Qu^4 - 2 dF Nu Qu^6 - 44/9 dF^2 Nu^2 Qu^6, 
		{0, 2} -> 4 dF Nd Qd^4 + 4 Nl Qe^4 + 4 dF Nu Qu^4, 
		{0, 1} -> 4/3 dF Nd Qd^2 + (4 Nl Qe^2)/3 + 4/3 dF Nu Qu^2	
		     } /. Casimirs /. Charges;

Options[ RunQEDxQCD ] = { 
				"nL" -> 3, (* charged leptons *)
				"nU" -> 2, (* up-type quarks *)
				"nD" -> 3 (* down-type quarks *)
			}

RunQEDxQCD[as0_?NumericQ,ae0_?NumericQ, muStart_?NumericQ, muEnd_?NumericQ, L_Integer:4,OptionsPattern[]] := Module[{nL = OptionValue["nL"],nU = OptionValue["nU"], nD = OptionValue["nD"]}, 	
(*				
				Print[ runAA[as0,ae0,2*Log[ muEnd/muStart], L, nU,nD,nL]];
*)
(*Print["L=",L," nU=", nU," nD=",nD, " nL=",nL];*)
				Return[ run[as0,ae0,2*Log[ muEnd/muStart], L, nU,nD,nL]]]
(* t = Log[MuF^2/Muo^2] *)
run[as0_?NumericQ,ae0_?NumericQ, tf_?NumericQ, L_Integer:3,nu_:2,nd_:3,nl_:3] :=  Block[{aS,aL,t, baSc,baS, baLc, baL, eqs, sol}, 
			baSc = Select[ AsBetaCoefficients, (Plus @@ First[#] <= L)&]; (* restrict to nL-loop coefficients *) 
			baLc = Select[ AlBetaCoefficients, (Plus @@ First[#] <= L)&]; (* restrict to nL-loop coefficients *) 
			baS = FromCoefficientRules[ baSc, {aS[t]/(4 Pi), aL[t]/(4 Pi)}] /. {Nu -> nu, Nd->nd, Nl->nl} ; (* build beta-function for aS *)
			baL = FromCoefficientRules[ baLc, {aS[t]/(4 Pi), aL[t]/(4 Pi)}] /. {Nu -> nu, Nd->nd, Nl->nl}  ; (* build beta-function for aL *) 
			eqs = { D[aS[t],t] == aS[t] * baS, D[aL[t],t] == aL[t] * baL, aS[0] == as0, aL[0]==ae0}; (* RGEs *)
			sol = NDSolve[ eqs, {aS,aL},{t,0, tf}];	
			Return[Flatten[{aS[tf], aL[tf]} /. sol]]
]


End[]
EndPackage[]
