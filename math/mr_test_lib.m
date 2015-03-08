<<"mr_lib.m";


(*

sol = sol[ PDG`MT][]
chiSquare[sol]

*)

as = 0.108002; (* gs^2/(4 Pi) /. sol *)

pars = RunParsFromPoleMassesAndAsWithMb[][ as, PDG`MT] 

RunSM[pars, PDG`MZ]

checkMZ = ScaleDependence[MZ][pars]

PDG`MZ

ReportPDGValues[]

EstimateTheorUncertaintyInMatchingDegrassi[PDG`MT,10]
(* 
pars /. h->1 /. x_String->1 /. corr[x_,y_]:> y*corr[y-x]
*)
(*
solNEW[ PDG`MZ][]

*)
