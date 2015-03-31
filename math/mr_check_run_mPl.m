Install["mr"];


DegrassiPars = {	
	g1 -> 0.35830,
	g2 -> 0.64779,
	gs -> 1.1666,
	yb -> 0,
	yt -> 0.93690,
	lam -> 0.12604,
	m -> 131.55,
	scale -> 173.34};
DegrassiParsNom = {	
	g1 -> 0.35830,
	g2 -> 0.64779,
	gs -> 1.1666,
	yb -> 0,
	yt -> 0.93690,
	lam -> 0.12604,
	scale -> 173.34};

mPlank = 1.2209 * 10^(19);


mPlankPars = RunSM[ DegrassiPars, mPlank];
mPlankPars1 = RunSMcouplings[ DegrassiParsNom, mPlank, 3];
mPlankPars12 = RunSMcouplings[ DegrassiParsNom, mPlank, 2];
mPlankPars11 = RunSMcouplings[ DegrassiParsNom, mPlank, 1];


backtoMT = RunSM[ mPlankPars, 173.34];
backtoMT1 = RunSMcouplings[ mPlankPars1, 173.34,3];


Print[ DegrassiPars, backtoMT];

Print[DegrassiParsNom, backtoMT1];

Map[ (# /. backtoMT1) - (# /. DegrassiParsNom) &, {g1,g2,gs,yt,yb,lam}]
Map[ (# /. backtoMT) - (# /. DegrassiPars) &, {g1,g2,gs,yt,yb,lam, m}]

Map[ (# /. mPlankPars1) - (# /. mPlankPars12) &, {g1,g2,gs,yt,yb,lam}]
Map[ (# /. mPlankPars1) - (# /. mPlankPars11) &, {g1,g2,gs,yt,yb,lam}]

Quit[];


downandbackMT = RunSM[RunSM[DegrassiPars, 17.334], 173.34]


FixG1[pars_]:= pars /. (g1->x_):> (g1 -> x Sqrt[5/3]);

mPlankPars = FixG1[ mPlankPars];

CheckValue[x_,value_,tol_:0.0001][pars_] := If[ Abs[ (x /. pars) - value] < tol, Print[x, "=", x /. pars,"  checked against ", value], Print[  x, " checking failed: ", x /. pars, "  vs ", value]];

CheckValue[g1,0.6154][mPlankPars]
CheckValue[g2,0.5055][mPlankPars]
CheckValue[gs,0.4873][mPlankPars]
CheckValue[yt,0.3825][mPlankPars]
CheckValue[lam,-0.0143][mPlankPars]
CheckValue[m,129.4][mPlankPars]


DegrassiPars
backtoMT
downandbackMT
