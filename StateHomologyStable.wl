(* ::Package:: *)

AppendTo[$Path,NotebookDirectory[]];


BeginPackage["StateHomologyStable`",{"StateFunctorStable`","CechOpsStable`","BasicStable`"}];


(* ::Title:: *)
(*Descriptions of Public (Unhidden) Functions*)


endHomologyRk::usage="endHomologyRk[deg,rho,dimprim,cover] outputs the rank of the degree 'deg'
homology group associated to the state rho";


GNSHomologyRk::usage="GNSHomologyRk[deg,rho,dimprim,cover] outputs the rank of the degree 'deg'
homology group associated to the state rho";


endHomologyVect::usage="endHomologyVect[deg,rho,dimprim,cover] outputs a list of generators of 
the degree 'deg' homology group---identified with a vector space---associated to the state rho";


GNSHomologyVect::usage="GNSHomologyVect[deg,rho,dimprim,cover] outputs a list of generators of 
the degree 'deg' homology group---identified with a vector space---associated to the state rho";


endHomologyObj::usage="endHomologyObj[deg,rho,dimprim,cover] outputs a lists of generators of the degree
'deg' homology group---given by chains valued in elements of endFunctorObj---associated to the state rho.";


GNSHomologyObj::usage="GNSHomologyObj[deg,rho,dimprim,cover] outputs a lists of generators of the degree
'deg' homology group---given by chains valued in elements of GNSFunctorObj---associated to the state rho.";


endHomologyRkList::usage="endHomologyRkList[rho_,dimprim_,cover_] outputs a list of the ranks
of the homology groups (associated to the state rho) from degree -1 to degree Length[cover]-1.  This 
function is faster than sequentially using stateHomologyRk as it temporarily stores values
from endFunctorObj and endFunctorMor.";


GNSHomologyRkList::usage="GNSHomologyRkList[rho_,dimprim_,cover_] outputs a list of the ranks
of the homology groups (associated to the state rho) from degree -1 to degree Length[cover]-1.  This 
function is faster than sequentially using stateHomologyRk as it temporarily stores values
from GNSFunctorObj and GNSFunctorMor.";


(* ::Title:: *)
(*Function Definitions*)


Begin["`Private`"];


endHomologyRk[deg_,rho_,dimprim_,cover_]:=homologyRank[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyRk[deg_,rho_,dimprim_,cover_]:=homologyRank[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


endHomologyVect[deg_,rho_,dimprim_,cover_]:=homologyVect[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyVect[deg_,rho_,dimprim_,cover_]:=homologyVect[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


endHomologyObj[deg_,rho_,dimprim_,cover_]:=homologyObj[deg, endFunctorObj[rho,dimprim], endFunctorMor[rho,dimprim], cover, localInProd];


GNSHomologyObj[deg_,rho_,dimprim_,cover_]:=homologyObj[deg, GNSFunctorObj[rho,dimprim], GNSFunctorMor[rho,dimprim], cover, localInProd];


endHomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=endFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=endFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-1}]
];


GNSHomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=GNSFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=GNSFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-1}]
];


(* ::Title:: *)
(*End Matter*)


End[];


EndPackage[]



