<<"mr_lib.m";

$CharacterEncoding = "UTF8";

SetOptions[$Output, PageWidth->200]
(*

sol = sol[ PDG`MT][]
chiSquare[sol]

*)

(* SetPDGValues[$DegrassiEtAlValues] *)


ReportPDGValues[]


inipars = RunParsFromPoleMassesAndGf[1][PDG`MT]
iniparsn = Last /@ RunParsFromPoleMassesAndGf[1][PDG`MT]

Needs["NumericalCalculus`"]


Print["MT"];
mtderiv = ND[ RunParsFromPoleMassesAndGfList[1][mt,PDG`MT], mt, PDG`MT];
Print["MH"];
mhderiv = ND[ RunParsFromPoleMassesAndGfList[1][PDG`MT,mh,PDG`MT], mh, PDG`MH];
Print["MW"];
mwderiv = ND[ RunParsFromPoleMassesAndGfList[1][PDG`MT,PDG`MH,PDG`asQCD,mw,PDG`MT]*PDG`dMW, mw, PDG`MW];
Print["as"]
asderiv = ND[ RunParsFromPoleMassesAndGfList[][PDG`MT,PDG`MH,as,PDG`MT]*PDG`dasQCD, as, PDG`asQCD];

ress = MapThread[First[#1] -> #2 &, {inipars,iniparsn + df["MT" - PDG`MT] * mtderiv + df["MH" - PDG`MH] * mhderiv + df["MW" - PDG`MW,PDG`dMW] * mwderiv  + df["as" - PDG`asQCD,PDG`dasQCD] * asderiv}] 
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


EstimateTheorUncertaintyInMatchingDegrassi[PDG`MT,10]
EstimateTheorUncertaintyInMatchingDegrassi[PDG`MT,2]

(* 
pars /. h->1 /. x_String->1 /. corr[x_,y_]:> y*corr[y-x]
*)
(*
solNEW[ PDG`MZ][]

*)
