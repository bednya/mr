<<"mr_lib.m";

$CharacterEncoding = "UTF8";

SetOptions[$Output, PageWidth->200]
(*

sol = sol[ PDG`MT][]
chiSquare[sol]

*)

(* SetPDGValues[$DegrassiEtAlValues] *)


ReportPDGValues[]

matchscale = PDG`MT;

inipars = RunParsFromPoleMassesAndGf[1][matchscale]
iniparsn = Last /@ inipars;

Needs["NumericalCalculus`"]


Print["MT"];
mtderiv = ND[ RunParsFromPoleMassesAndGfList[1][mt,matchscale]*PDG`dMT, mt, PDG`MT];
Print["MH"];
mhderiv = ND[ RunParsFromPoleMassesAndGfList[1][PDG`MT,mh,matchscale]*PDG`dMH, mh, PDG`MH];
Print["MW"];
mwderiv = ND[ RunParsFromPoleMassesAndGfList[1][PDG`MT,PDG`MH,PDG`asQCD,mw,matchscale]*PDG`dMW, mw, PDG`MW];
Print["as"]
asderiv = ND[ RunParsFromPoleMassesAndGfList[][PDG`MT,PDG`MH,as,matchscale]*PDG`dasQCD, as, PDG`asQCD];

ress = MapThread[First[#1] -> #2 &, {inipars,iniparsn + df["MT" - PDG`MT,PDG`dMT] * mtderiv + df["MH" - PDG`MH,PDG`dMH] * mhderiv + df["MW" - PDG`MW,PDG`dMW] * mwderiv  + df["as" - PDG`asQCD,PDG`dasQCD] * asderiv}] 

ress >> "pdg2014_runpars_MT"

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
