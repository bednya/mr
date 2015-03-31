BeginPackage["ASdecQCDandEW`"];

Print[" A Mathematica package for two-loop matching of the strong coupling constant"];
Print["   from the SM to 5-flavour QCD x QED effective theory by A.V. Bednyakov." ];
Print[" Please, see Information on DecAS5up6OS, DecAS6down5OS, and $SMdata."];
Print[" Three-loop QCD corrections are taken from: "]; 
Print["   K.G. Chetyrkin, Johann H. Kuhn, M. Steinhauser, Comput.Phys.Commun. 133 (2000) 43-65"];

DEC::WrongNumberOfLoops  = "Wrong number of loops for running nl = `1`. Should be 1 < nl < 4";
DEC::NonNumericParameters  = "The parameter `1` should be a number";
OnlyQCD::usage = "An option to DecAS5up6OS to output only QCD contribution";
SMdata::usage = " An option to specify SM parameters. A list of rules should be given, e.g., { MTpole -> 172.5, etc.}";
$SMdata::usage = "Default values of SM parameters (PDG 2014)";
ContribTagging::usage = "An option to DecAS5up6OS/DecAS6down5OS to tag all the contributions (NB: effectively, the tags are equal to one) ";

DecAS5up6OS::usage = "
	DecAS5up6OS[ asQCD, mu, nl, options]
		given \[Alpha]_s(mu) in QCD x QED  returns \[Alpha]_s(mu) in the SM	

	asQCD -  running \[Alpha]_s(mu) in nf = 5 QCD, 
	mu - matching scale
	nl - number of loops: 
	     nl = 1 (tree-level decoupling), 
	     nl = 2 (1-loop decoupling), 
	     nl = 3 (2-loop decoupling),
	     nl = 4 (3-loop decoupling)

	options - {
		   OnlyQCD -> True/False, 
		   ContribTagging -> True/False,
		   SMdata -> List of rules with SM parameters (including pole top-quark and W/Z/H boson masses)
		  }
	";

DecAS6down5OS::usage = "
	DecAS6down5OS[ asSM, mu, nl, options]
		given \[Alpha]_s(mu) in the SM returns \[Alpha]_s(mu) in QCD x QED	
	asSM -  running \[Alpha]_s(mu) in the SM
	mu - matching scale
	nl - number of loops: 
	     nl = 1 (tree-level decoupling), 
	     nl = 2 (1-loop decoupling), 
	     nl = 3 (2-loop decoupling),
	     nl = 4 (3-loop decoupling)
	options - {
		   OnlyQCD -> True/False, 
		   ContribTagging -> True/False,
		   SMdata -> List of rules with SM parameters (including pole top-quark and W/Z/H boson masses) 
		  }
	";

ASqcdNf5atMZ::usage = "\*SubscriptBox[\[Alpha], s](\*SubscriptBox[M,Z]) in \*SubscriptBox[N,f] = 5 QCD"; 
invAEatMZ::usage = "1/\[Alpha](\*SubscriptBox[M,Z]) in \*SubscriptBox[N,f] = 5 QCD"; 
SW2NDatMZ::usage = "Running (non-decoupled) value of \*SuperscriptBox[Sin,2](\*SubscriptBox[\[Theta],W]) at \*SubscriptBox[M,Z]";
SW2OS::usage = "ON-shell  value of \*SuperscriptBox[Sin,2](\*SubscriptBox[\[Theta],W])";
MTpole::usage = "Top quark pole mass";
MWpole::usage = "W-boson pole mass";
MZpole::usage = "Z-boson pole mass";
MHpole::usage = "Higgs boson pole mass";

delta::usage = "delta[ SM parameter ] represents experimental uncertainty of the parameter";

Begin["Private`"]

(* DebugPrint = Print; *)

(*
Phi[z_?NumericQ] := 
 Sqrt[z/(z - 1)]*(-4 PolyLog[2, 1/2*(1 - Sqrt[1 - 1/z])] + 
     2 Log[1/2*(1 - Sqrt[1 - 1/z])]^2 - Log[4 z]^2 + Pi^2/3)
*)

(* Phi-function entering two-loop bubble with two equal masses *)

ClausenCl[n_Integer?EvenQ,x_] := (PolyLog[n, E^(-I x)] - PolyLog[n, E^(I x)])I/2
ClausenCl[n_Integer?OddQ,x_] :=  (PolyLog[n, E^(-I x)] + PolyLog[n, E^(I x)])/2

Phi[z_?NumberQ] := N[4 * Sqrt[ z/(1-z) ] * ClausenCl[2, 2 * ArcSin[Sqrt[z]] ]];

$PDGdata[2014] = { 
		  ASqcdNf5atMZ -> 0.1185, delta[ASqcdNf5atMZ] -> 0.0006, 
		  invAEatMZ -> 127.940, delta[invAEatMZ] -> 0.014, (* from Review on EW physics *)
		  SW2NDatMZ -> 0.23144, (* from Review on EW physics *) 
		  SW2OS -> 0.22333,(* from Review on EW physics *)
                  MTpole -> 173.21,  delta[MTpole] -> Sqrt[ 0.51^2 + 0.71^2], (* GeV *)
		  MWpole -> 80.385,  delta[MWpole] -> 0.015, (* GeV *)
		  MZpole -> 91.1876, delta[MZpole] -> 0.0021, (* GeV *)
		  MHpole -> 125.7,   delta[MHpole] -> 0.4 (* GeV *)
		 };

$SMdata = $PDGdata[2014];

(* 
   Expansion of 1/xiOS^2 and xiOS^2 in powers of 
	{AS,AEW,LT}
   where 
	(4 Pi) AS = GS^2/(4 Pi) - Running QCD coupling defined in nf=5 QCDxQED
	(4 Pi) AEW = EL^2/(4 Pi)/SW^2*mmt/mmW
	LT = Log[ mmt / mmu ], mmu - decoupling scale squared

  FromCoefficientRules[ oneloop, {AS, AEW, LT} ] restores original expression
  
  
  Xht = mmh/mmt
  Xwt = mmW/mmt
  Xzt = mmZ/mmt

 *)


(* ASQCD = ASSM * xiASOS *)
XiAsOSoneloop = {{1, 0, 1} -> 2/3};
XiAsOStwoloop = {
	{2, 0, 2} -> 4/9, 
	{2, 0, 1} -> 38/3, 
	{2, 0, 0} -> -14/3, 
 	{1, 1, 1} -> -1 + (2*Xwt)/9 + (22*Xwt^2)/(9*Xzt) + (11*Xzt)/6, 
 	{1, 1, 0} -> (-4*(-105 + 9*Xht + 25*Xwt*(3 + 2*Xwt)) + 
     (24*(5 - 2*Xwt)^2)/(-4 + Xzt) - (108*Xwt^2)/Xzt + (-291 + 160*Xwt)*Xzt - 
     68*Xzt^2)/216 + ((-4 + Xht)*Sqrt[-((-4 + Xht)*Xht)]*ArcCos[Sqrt[Xht]/2])/
    6 + (Sqrt[4 - Xzt]*(32*Xwt^2*(2 + Xzt) - 40*Xwt*Xzt*(2 + Xzt) + 
      Xzt^2*(7 + 17*Xzt))*ArcCos[Sqrt[Xzt]/2])/(54*Sqrt[Xzt]) + 
   ((6 + Xht*(27 + (-10 + Xht)*Xht))*Log[Xht])/(12*(-4 + Xht)) - 
   ((-1 + Xwt)^2*(1 + 2*Xwt)*Log[Abs[-1 + Xwt]])/6 + 
   (Xwt*(18 + 6/(-1 + Xwt) + Xwt*(-3 + 2*Xwt))*Log[Xwt])/6 + 
   ((8*Xwt^2*(384 + Xzt*(-216 + Xzt*(85 + 4*(-8 + Xzt)*Xzt))) - 
      4*Xwt*Xzt*(744 + Xzt*(-432 + Xzt*(199 + 10*(-8 + Xzt)*Xzt))) + 
      Xzt^2*(2226 + Xzt*(-1701 + Xzt*(638 + Xzt*(-163 + 17*Xzt)))))*Log[Xzt])/
    (108*(-4 + Xzt)^2*Xzt) - ((-1 + Xht)*Phi[Xht/4])/(2*(-4 + Xht)*Xht) + 
   ((32*Xwt^2*(-2 + Xzt) - 40*Xwt*(-2 + Xzt)*Xzt + Xzt^2*(-7 + 8*Xzt))*
     Phi[Xzt/4])/(6*(-4 + Xzt)^2*Xzt^2)};


(* Chetyrkin (RunDec manual) eq. 22, NB: L -> -LT  *)

nl = 5; (* 5-light flavours *)

XiAsOSthreeloop = { 
	{3,0,3} -> -4^3 ( -1/216	),
	{3,0,2} ->  4^3 ( -131/576      ),
	{3,0,1} -> -4^3 ( -8521/1728 + nl 409/1728 ),
	{3,0,0} ->  4^3 ( -58933/124416 - 2/3 Zeta[2] (1 + 1/3 Log[2]) - 80507/27648 Zeta[3] + nl ( 2479/31104 + Zeta[2]/9 ) )
};

(* ASSM = ASQCD * invxiASOS *)

InvXiAsOSoneloop = {{1, 0, 1} -> -2/3};
 
InvXiAsOStwoloop = {
     {2, 0, 2} -> 4/9, 
     {2, 0, 1} -> -38/3, 
     {2, 0, 0} -> 14/3, 
     {1, 1, 1} -> 1 - (2*Xwt)/9 - (22*Xwt^2)/(9*Xzt) - (11*Xzt)/6, 
     {1, 1, 0} -> (36*Xht + 20*(-21 + 5*Xwt*(3 + 2*Xwt)) - 
         (24*(5 - 2*Xwt)^2)/(-4 + Xzt) + (108*Xwt^2)/Xzt + 
         (291 - 160*Xwt)*Xzt + 68*Xzt^2)/216 - 
       ((-4 + Xht)*Sqrt[-((-4 + Xht)*Xht)]*ArcCos[Sqrt[Xht]/2])/6 + 
       ((-4 + Xzt)*(32*Xwt^2*(2 + Xzt) - 40*Xwt*Xzt*(2 + Xzt) + 
          Xzt^2*(7 + 17*Xzt))*ArcCos[Sqrt[Xzt]/2])/
        (54*Sqrt[-((-4 + Xzt)*Xzt)]) - ((6 + Xht*(27 + (-10 + Xht)*Xht))*
         Log[Xht])/(12*(-4 + Xht)) + ((-1 + Xwt)^2*(1 + 2*Xwt)*Log[(* NB *) Abs[-1 + Xwt]])/
        6 + (Xwt*(12 + Xwt*(-21 + (5 - 2*Xwt)*Xwt))*Log[Xwt])/
        (6*(-1 + Xwt)) + 
       ((Xzt^2*(-2226 + Xzt*(1701 + Xzt*(-638 + (163 - 17*Xzt)*Xzt))) - 
          8*Xwt^2*(384 + Xzt*(-216 + Xzt*(85 + 4*(-8 + Xzt)*Xzt))) + 
          4*Xwt*Xzt*(744 + Xzt*(-432 + Xzt*(199 + 10*(-8 + Xzt)*Xzt))))*
         Log[Xzt])/(108*(-4 + Xzt)^2*Xzt) + ((-1 + Xht)*Phi[Xht/4])/
        (2*(-4 + Xht)*Xht) + ((-32*Xwt^2*(-2 + Xzt) + 40*Xwt*(-2 + Xzt)*Xzt + 
          (7 - 8*Xzt)*Xzt^2)*Phi[Xzt/4])/(6*(-4 + Xzt)^2*Xzt^2)
     };

(* Chetyrkin et al. (RunDec manual) eq. 25, NB: L -> -LT  *)

InvXiAsOSthreeloop = {
	{3,0,3} -> -4^3 ( +1/216	),
	{3,0,2} ->  4^3 ( +511/576      ),
	{3,0,1} -> -4^3 ( +8941/1728 - nl 409/1728 ),
	{3,0,0} ->  4^3 ( 58933/124416 + 2/3*Zeta[2]*(1 + 1/3 Log[2]) + 80507/27648 Zeta[3] + nl ( -2479/31104 - Zeta[2]/9 ) )
	};


Options[DecAS5up6OS] = {OnlyQCD -> False, 
			 ContribTagging -> False,
			 SMdata -> $SMdata};


DecAS5up6OS[asQCD_?NumberQ, al__?NumericQ, mu_, nl_?IntegerQ, OptionsPattern[] ] := 
		Module[{invXiAs = 1, 
			smdata = OptionValue[ SMdata ],
			AS = asQCD/(4 Pi), 
			AEW, (* dummy variable (hack, since FromCoefficientRules returns Indetermined if one of variables is zero) *) aew, 
			LT,  
			lt,
			xwt,
			xzt,
			xht},
			AEW = 1/(4 Pi) * al * MTpole^2/MWpole^2/(1-MWpole^2/MZpole^2) /. smdata;
			Print["HELLLLLO:", AEW];
			LT = Log[ MTpole^2/mu^2] /. smdata;
			If [ Not[NumericQ[ AS  ]] , Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)",*) AS]; Return[$Failed]];
			If [ Not[NumericQ[ AEW ]],  Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)", *)AEW]; Return[$Failed]];
			If [ (nl<1) ||  (nl > 4), Message[DEC::WrongNumberOfLoops, nl]; Return[$Failed]]; 
			If [OptionValue[ ContribTagging ] === True,
				AS = AS * "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)";
				AEW = AEW * "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)";
			];
		 	If [ (nl>1), invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSoneloop, {AS, AEW, LT} ]];
			If[ OptionValue[ OnlyQCD ] === True, 
				AEW = aew,
				(* else *)
				xwt = MWpole^2/MTpole^2 /. smdata;
				xzt = MZpole^2/MTpole^2 /. smdata;
				xht = MHpole^2/MTpole^2 /. smdata;
			];
		 	If [ (nl>2), 
				DebugPrint["Coupling values: as/(4pi) = ", AS , " aew/(4pi) = ",AEW , ", LT = ", LT];
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOStwoloop, {AS, AEW, lt} ] /. lt->LT /. aew -> 0 /. {Xwt -> xwt, Xzt -> xzt, Xht -> xht};
			   ];
			If [ (nl>3), (* QCD ONLY *)
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSthreeloop, {AS,AEW, lt} ] /. lt->LT
				];
			DebugPrint[ {asQCD, Chop[invXiAs]} ];
			Return[ Times @@ {asQCD, Chop[invXiAs]} ];
		];	

DecAS5up6OS[asQCD_?NumberQ, al_?NumericQ, mu_, nl_?IntegerQ, mt_?NumericQ,mw_?NumericQ,mz_?NumericQ, mh_?NumericQ, OptionsPattern[] ] := 
		Module[{invXiAs = 1, 
			smdata = {MTpole -> mt, MWpole -> mw, MZpole -> mz, MHpole -> mh},
			AS = asQCD/(4 Pi), 
			AEW, (* dummy variable (hack, since FromCoefficientRules returns Indetermined if one of variables is zero) *) aew, 
			LT,  
			lt,
			xwt,
			xzt,
			xht},
			AEW = 1/(4 Pi) * al * MTpole^2/MWpole^2/(1-MWpole^2/MZpole^2) /. smdata;
			(*Print["HELLLLLO:", AEW];*)
			LT = Log[ MTpole^2/mu^2] /. smdata;
			If [ Not[NumericQ[ AS  ]] , Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)",*) AS]; Return[$Failed]];
			If [ Not[NumericQ[ AEW ]],  Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)", *)AEW]; Return[$Failed]];
			If [ (nl<1) ||  (nl > 4), Message[DEC::WrongNumberOfLoops, nl]; Return[$Failed]]; 
			If [OptionValue[ ContribTagging ] === True,
				AS = AS * "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)";
				AEW = AEW * "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)";
			];
		 	If [ (nl>1), invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSoneloop, {AS, AEW, LT} ]];
			If[ OptionValue[ OnlyQCD ] === True, 
				AEW = aew,
				(* else *)
				xwt = MWpole^2/MTpole^2 /. smdata;
				xzt = MZpole^2/MTpole^2 /. smdata;
				xht = MHpole^2/MTpole^2 /. smdata;
			];
		 	If [ (nl>2), 
				DebugPrint["Coupling values: as/(4pi) = ", AS , " aew/(4pi) = ",AEW , ", LT = ", LT];
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOStwoloop, {AS, AEW, lt} ] /. lt->LT /. aew -> 0 /. {Xwt -> xwt, Xzt -> xzt, Xht -> xht};
			   ];
			If [ (nl>3), (* QCD ONLY *)
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSthreeloop, {AS,AEW, lt} ] /. lt->LT
				];
			DebugPrint[ {asQCD, Chop[invXiAs]} ];
			Return[ Times @@ {asQCD, Chop[invXiAs]} ];
		];	

DecAS5up6OS[asQCD_, mu_, nl_?IntegerQ, OptionsPattern[] ] := 
		Module[{invXiAs = 1, 
			smdata = OptionValue[ SMdata ],
			AS = asQCD/(4 Pi), 
			AEW, (* dummy variable (hack, since FromCoefficientRules returns Indetermined if one of variables is zero) *) aew, 
			LT,  
			lt,
			xwt,
			xzt,
			xht},
			AEW = 1/(4 Pi) * (1/invAEatMZ * MTpole^2/MWpole^2/SW2NDatMZ) /. smdata;
			Print["HELLLLLO:", AEW];
			Print["Check:", 1/(4 Pi) * (1/invAEatMZ * MTpole^2/MWpole^2/(1-MWpole^2/MZpole^2)) /. smdata];
			LT = Log[ MTpole^2/mu^2] /. smdata;
			If [ Not[NumericQ[ AS  ]] , Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)",*) AS]; Return[$Failed]];
			If [ Not[NumericQ[ AEW ]],  Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)", *)AEW]; Return[$Failed]];
			If [ (nl<1) ||  (nl > 4), Message[DEC::WrongNumberOfLoops, nl]; Return[$Failed]]; 
			If [OptionValue[ ContribTagging ] === True,
				AS = AS * "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)";
				AEW = AEW * "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)";
			];
		 	If [ (nl>1), invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSoneloop, {AS, AEW, LT} ]];
			If[ OptionValue[ OnlyQCD ] === True, 
				AEW = aew,
				(* else *)
				xwt = MWpole^2/MTpole^2 /. smdata;
				xzt = MZpole^2/MTpole^2 /. smdata;
				xht = MHpole^2/MTpole^2 /. smdata;
			];
		 	If [ (nl>2), 
				DebugPrint["Coupling values: as/(4pi) = ", AS , " aew/(4pi) = ",AEW , ", LT = ", LT];
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOStwoloop, {AS, AEW, lt} ] /. lt->LT /. aew -> 0 /. {Xwt -> xwt, Xzt -> xzt, Xht -> xht};
			   ];
			If [ (nl>3), (* QCD ONLY *)
				invXiAs = invXiAs + FromCoefficientRules[ InvXiAsOSthreeloop, {AS,AEW, lt} ] /. lt->LT
				];
			DebugPrint[ {asQCD, Chop[invXiAs]} ];
			Return[ Times @@ {asQCD, Chop[invXiAs]} ];
		];	


Options[DecAS6down5OS] = {OnlyQCD -> False, 
			 ContribTagging -> False,
			 SMdata -> $SMdata};

DecAS6down5OS[asQCD_, mu_, nl_?IntegerQ, OptionsPattern[] ] := 
		Module[{XiAs = 1, 
			smdata = OptionValue[ SMdata ],
			AS = asQCD/(4 Pi), 
			AEW, (* dummy variable (hack, since FromCoefficientRules returns Indetermined if one of variables is zero) *) aew, 
			lt,
			LT,  
			xwt,
			xzt,
			xht},
			AEW = 1/(4 Pi) * (1/invAEatMZ * MTpole^2/MWpole^2/SW2NDatMZ) /. smdata;
			LT = Log[ MTpole^2/mu^2] /. smdata;
			If [ Not[NumericQ[ AS  ]] , Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)",*) AS]; Return[$Failed]];
			If [ Not[NumericQ[ AEW ]],  Message[DEC:NonNumericParameters,(* "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)", *)AEW]; Return[$Failed]];
			If [ (nl<1) ||  (nl > 4), Message[DEC::WrongNumberOfLoops, nl]; Return[$Failed]]; 
			If [OptionValue[ ContribTagging ] === True,
				AS = AS * "\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)";
				AEW = AEW * "\!\(\*SubscriptBox[\(\[Alpha]\), \(ew\)]\)";
			];
		 	If [ (nl>1), XiAs =XiAs + FromCoefficientRules[ XiAsOSoneloop, {AS, AEW, LT} ]];
			If[ OptionValue[ OnlyQCD ] === True, 
				AEW = aew,
				(* else *)
				xwt = MWpole^2/MTpole^2 /. smdata;
				xzt = MZpole^2/MTpole^2 /. smdata;
				xht = MHpole^2/MTpole^2 /. smdata;
			];
		 	If [ (nl>2), 
				DebugPrint["Coupling values: as/(4pi) = ", AS , " aew/(4pi) = ",AEW , ", LT = ", LT];
				XiAs = XiAs + FromCoefficientRules[ XiAsOStwoloop, {AS, AEW, lt} ] /. lt->LT /. aew -> 0 /. {Xwt -> xwt, Xzt -> xzt, Xht -> xht};
			   ];
			If [ (nl>3), (* QCD ONLY *)
				XiAs = XiAs + FromCoefficientRules[ XiAsOSthreeloop, {AS,AEW,lt} ] /. lt->LT;
				];
			DebugPrint[ {asQCD, Chop[XiAs]} ];
			Return[ Times @@ {asQCD, Chop[XiAs]} ];
		];	


End[]
EndPackage[]
