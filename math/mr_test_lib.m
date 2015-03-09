<<"mr_lib.m";

$CharacterEncoding = "UTF8";

SetOptions[$Output, PageWidth->200]
(*

sol = sol[ PDG`MT][]
chiSquare[sol]

*)

as = 0.108002; (* gs^2/(4 Pi) /. sol *)

SetPDGValues[$DegrassiEtAlValues]
ReportPDGValues[]

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

(* 
pars /. h->1 /. x_String->1 /. corr[x_,y_]:> y*corr[y-x]
*)
(*
solNEW[ PDG`MZ][]

*)
