(* ::Package:: *)

AppendTo[$Path,NotebookDirectory[]];


BeginPackage["StateHomologyStable`",{"StateFunctorStable`","CechOpsStable`","BasicStable`"}];


(* ::Title:: *)
(*Descriptions of Public (Unhidden) Functions*)


comHomologyRk::usage="comHomologyRk[deg,rho,dimprim,cover] outputs the rank of the degree 'deg'
homology group associated to the state rho";


GNSHomologyRk::usage="GNSHomologyRk[deg,rho,dimprim,cover] outputs the rank of the degree 'deg'
homology group associated to the state rho";


comHomologyVect::usage="comHomologyVect[deg,rho,dimprim,cover] outputs a list of generators of 
the degree 'deg' homology group---identified with a vector space---associated to the state rho";


GNSHomologyVect::usage="GNSHomologyVect[deg,rho,dimprim,cover] outputs a list of generators of 
the degree 'deg' homology group---identified with a vector space---associated to the state rho";


comHomologyObj::usage="comHomologyObj[deg,rho,dimprim,cover] outputs a lists of generators of the degree
'deg' homology group---given by chains valued in elements of endFunctorObj---associated to the state rho.";


GNSHomologyObj::usage="GNSHomologyObj[deg,rho,dimprim,cover] outputs a lists of generators of the degree
'deg' homology group---given by chains valued in elements of GNSFunctorObj---associated to the state rho.";


comHomologyRkList::usage="comHomologyRkList[rho_,dimprim_,cover_] outputs a list of the ranks
of the homology groups (associated to the state rho) from degree -1 to degree Length[cover]-1.  This 
function is faster than sequentially using stateHomologyRk as it temporarily stores values
from endFunctorObj and endFunctorMor.";


GNSHomologyRkList::usage="GNSHomologyRkList[rho_,dimprim_,cover_] outputs a list of the ranks
of the homology groups (associated to the state rho) from degree -1 to degree Length[cover]-1.  This 
function is faster than sequentially using stateHomologyRk as it temporarily stores values
from GNSFunctorObj and GNSFunctorMor.";


(* ::Title:: *)
(*Function Definitions*)


(* ::Subsection:: *)
(*Homologically Graded*)


Begin["`Private`"];


comHomologyRk[rho_,dimprim_,cover_][deg_]:=homologyRank[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyRk[rho_,dimprim_,cover_][deg_]:=homologyRank[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


comHomologyVect[rho_,dimprim_,cover_][deg_]:=homologyVect[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyVect[rho_,dimprim_,cover_][deg_]:=homologyVect[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


comHomologyObj[rho_,dimprim_,cover_][deg_]:=homologyObj[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyObj[rho_,dimprim_,cover_][deg_]:=homologyObj[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


comHomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=endFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=endFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-1}]
];


GNSHomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=GNSFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=GNSFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-1}]
];


(* ::Subsection:: *)
(*Cohomologically Graded*)


fixDegree[deg_,cover_]:=Length@cover-deg-1;


comCohomologyRk[rho_,dimprim_,cover_][deg_]:=comHomologyRk[rho,dimprim,cover]@fixDegree;


GNSCohomologyRk[rho_,dimprim_,cover_][deg_]:=GNSHomologyRk[rho,dimprim,cover]@fixDegree;


comCohomologyVect[rho_,dimprim_,cover_][deg_]:=comHomologyVect[rho,dimprim,cover]@fixDegree;


GNSCohomologyVect[rho_,dimprim_,cover_][deg_]:=comHomologyVect[rho,dimprim,cover]@fixDegree;


comCohomologyObj[rho_,dimprim_,cover_][deg_]:=comCohomologyObj[rho,dimprim,cover]@fixDegree;


GNSCohomologyObj[rho_,dimprim_,cover_][deg_]:=comCohomologyObj[rho,dimprim,cover]@fixDegree;


comCohomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=endFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=endFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-1}]
];


GNSCohomologyRkList[rho_,dimprim_,cover_]:=Reverse@*GNSHomologyRkList;


(* ::Title:: *)
(*End Matter*)


End[];


EndPackage[]



