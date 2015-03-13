<<"mr_lib.m";

$CharacterEncoding = "UTF8";

SetOptions[$Output, PageWidth->200]
(*

sol = sol[ PDG`MT][]
chiSquare[sol]

*)

(* SetPDGValues[$DegrassiEtAlValues]   *)


ReportPDGValues[]



Needs["NumericalCalculus`"]

Do[ 




inipars = RunParsFromPoleMassesAndGfalpha[1][matchscale];
iniparsn = Last /@ inipars;


Print["MT"];
mtderiv = ND[ RunParsFromPoleMassesAndGfalphaList[1][mt,matchscale]*PDG`dMT, mt, PDG`MT];
Print["MH"];
mhderiv = ND[ RunParsFromPoleMassesAndGfalphaList[1][PDG`MT,mh,matchscale]*PDG`dMH, mh, PDG`MH];
Print["MW"];
mwderiv = ND[ RunParsFromPoleMassesAndGfalphaList[1][PDG`MT,PDG`MH,PDG`asQCD,mw,matchscale]*PDG`dMW, mw, PDG`MW];
Print["as"];
asderiv = ND[ RunParsFromPoleMassesAndGfalphaList[][PDG`MT,PDG`MH,as,matchscale]*PDG`dasQCD, as, PDG`asQCD, Scale -> 0.1 ]; (* scale prevents complex derivative of yb *)



uncert10 = EstimateTheorUncertaintyInMatchingDegrassi[matchscale,10];
uncert2 = EstimateTheorUncertaintyInMatchingDegrassi[matchscale,2];


ress = MapThread[First[#1] -> #2 &, {inipars,iniparsn + df["MT" - PDG`MT,PDG`dMT] * mtderiv 
						      + df["MH" - PDG`MH,PDG`dMH] * mhderiv 
						      + df["MW" - PDG`MW,PDG`dMW] * mwderiv  
						      + df["as" - PDG`asQCD,PDG`dasQCD] * asderiv}];


ress = ress /. (x_->y_):> (x-> Chop[y + "fc[2]"*(x /. uncert2 /. x->0) + "fc[10]"*(x /. uncert10 /. x->0)] );

PutAppend[ress, "pdg2014.out"],
{ matchscale, {50,60,70,80,PDG`MZ,100,110,120,130,140,150,160,PDG`MT,190,200,210,220,230,250, 260,270}}]; 
(* { matchscale, {PDG`MZ,PDG`MH,PDG`MT,2*PDG`MT,5*PDG`MT}}]; *)
(* { matchscale, {200,300}}]; *)

Quit[]


XtQCD[PDG`MB,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,PDG`MT]


(*
pars0 = sol[PDG`MT][]
parsp10 = RunParsFromPoleMassesAndAsWithMb[][ as, PDG`MT*10] 
parsm10 = RunParsFromPoleMassesAndAsWithMb[][ as, PDG`MT/10] 

pars1 = RunParsFromPoleMassesAndAsAndGf[][ as, PDG`MT] 

Quit[];
RunSM[pars, PDG`MZ]

checkMZ = ScaleDependence[MZ][pars]

PDG`MZ

*)



(* 
pars /. h->1 /. x_String->1 /. corr[x_,y_]:> y*corr[y-x]
*)
(*
solNEW[ PDG`MZ][]

*)
