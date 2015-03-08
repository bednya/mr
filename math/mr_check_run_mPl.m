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

mPlank = 1.2209 * 10^(19);


mPlankPars = RunSM[ DegrassiPars, mPlank];

FixG1[pars_]:= pars /. (g1->x_):> (g1 -> x Sqrt[5/3]);

mPlankPars = FixG1[ mPlankPars];

CheckValue[x_,value_,tol_:0.0001][pars_] := If[ Abs[ (x /. pars) - value] < tol, Print[x, " checked against ", value], Print[  x, " checking failed: ", x /. pars, "  vs ", value]];

CheckValue[g1,0.6154][mPlankPars]
CheckValue[g2,0.5055][mPlankPars]
CheckValue[gs,0.4873][mPlankPars]
CheckValue[yt,0.3825][mPlankPars]
CheckValue[lam,-0.0143][mPlankPars]
CheckValue[m,129.4][mPlankPars]

