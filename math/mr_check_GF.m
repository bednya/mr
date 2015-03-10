Get["mr_lib.m"];


RunParsFromPoleMassesAndAsAndGf[][RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT],PDG`MT]

RunParsFromPoleMassesAndAsWithMb[][RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT],PDG`MT]


Quit[];
DiffGF[looptag_:1][gf0_,as_,al_,mu_,nL_:2,nH_:1]:= Module[{mmW = PDG`MW^2, mmZ = PDG`MZ^2, 
			xw = XW[PDG`MB,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,mu,nL,nH],  
			xz = XZ[PDG`MB,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,mu,nL,nH], 
			dyW,dyZ, aEW = al/(4 Pi), aQCD = as/(4 Pi)},
				dyZ = 1 + looptag*"aew"*aEW*yZ[1,0]+ looptag^2 ("aew"*"as"*aEW*aQCD*yZ[1,1]+ "aew"^2 aEW^2*yZ[2,0])/.xz;
				dyW = 1 + looptag*"aew"*aEW*yW[1,0]+ looptag^2 ("aew"*"as"*aEW*aQCD*yW[1,1]+ "aew"^2 aEW^2*yW[2,0])/.xw;
				res = 10^5 * (gf0 - xxx al*Pi/Sqrt[2]/mmW/dyW/(1 - mmW/mmZ * dyW/dyZ));
				If[ NumericQ[looptag], res = res /. x_String -> 1];
				Return[ res ]
				];

xxx =1;
DiffGF[1][ PDG`GF, RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT], 1/127.56, PDG`MT] //Expand
DiffGF[1][ PDG`GF, RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT], 1/127.98, PDG`MT] //Expand
DiffGF[1][ PDG`GF, RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT], 1/127.7562045, PDG`MT] //Expand
DiffGF[1][ PDG`GF, RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT], 1/127.72, PDG`MT] //Expand
RunQCDnf6[PDG`MT,PDG`asQCD,PDG`MZ,4,PDG`MT]

