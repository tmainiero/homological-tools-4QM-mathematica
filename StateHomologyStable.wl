(* ::Package:: *)

AppendTo[$Path,NotebookDirectory[]];


BeginPackage["StateHomologyStable`",{"StateFunctorStable`","CechOpsStable`","BasicStable`"}];


(* ::Title:: *)
(*Descriptions of Public (Unhidden) Functions*)


comCohomologyRk::usage="Outputs a list of ranks of commutant cohomology groups of a multipartite state,
read left to right: the first component is the dimension of the degree 0 component
the second is the degree 1, ...,  the N-1th is the degree N-1 (where N is the number of tensor factors).

This function takes in several possible classes of inputs:
--comCohomologyRk[rho,dimprim][deg]: outputs the rank/dimension of the degree 'deg'  commutant
cohomology component associated to the multipartite density state 'rho' on a set of tensor
factors with Hilbert space dimensions given by the ordered list 'dimprim'.
**'deg' is an integer
**'rho' is a a  positive semidefinite matrix)
**'dimprim' is a list of integers (e.g. {2,2,3})

--comCohomologyRk[rho,dimprim,partition][deg]: outputs the dimension of the degree 'deg' commutant
cohomology component associated to a coarsening of the multipartite density state 'rho' defined by
the partition 'partition 'and living on the set of tensor factors 'dimprim' (before the coarsening).
Here, 'deg', 'rho', and 'dimprim' are as before, and
**'partition' is a list of collections of subsystems: e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem.

--comCohomologyRk[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comCohomologyRk[rho,{2,...,2}][deg].
**'numSys' is an integer

--comCohomologyRk[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comCohomologyRk[rho,{d,...,d}][deg].
**'d' is an integer.";


GNSCohomologyRk::usage="Outputs the rank of a fixed GNS cohomology component of a multipartite state.

This function takes in several possible classes of inputs:
--GNSCohomologyRk[rho,dimprim][deg]: outputs the rank/dimension of the degree 'deg' GNS
cohomology component associated to the multipartite density state 'rho' on a set of tensor
factors with Hilbert space dimensions given by the ordered list 'dimprim'.
**'deg' is an integer
**'rho' is a a  positive semidefinite matrix)
**'dimprim' is a list of integers (e.g. {2,2,3})

--GNSCohomologyRk[rho,dimprim,partition][deg]: outputs the dimension of the degree 'deg' GNS
cohomology component associated to a coarsening of the multipartite density state 'rho' defined by
the partition of tensor factors 'partition 'and living on the set of tensor factors 'dimprim' (before the coarsening).
Here, 'deg', 'rho', and 'dimprim' are as before, and
**'partition' is a list of collections of subsystems: e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem.

--GNSCohomologyRk[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSCohomologyRk[rho,{2,...,2}][deg].
**'numSys' is an integer

--GNSCohomologyRk[deg,rho,numSys,d]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSCohomologyRk[rho,{d,...,d}][deg].
**'d' is an integer.";


(* comCohomologyVect::usage="Outputs a list of generators of 
the commutant cohomology component of fixed degree---using the Frobenius inner product on matrices
this component is identified with a subspace of R^{d}, where d is the dimension of the cohomology component.
The output is a list of vectors whose span is this subspace.

This function takes in several possible classes of inputs:
--comCohomologyVect[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--comCohomologyVect[rho,dimprim,partition][deg]: computes the commutant cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite commutant cohomology.

--comCohomologyVect[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comCohomologyVect[rho,{2,...,2}][deg].
**'numSys' is an integer

--comCohomologyVect[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comCohomologyVect[rho,{d,...,d}][deg].
**'d' is an integer."; 

!!!FUNCTION COMMENTED DUE TO SIGN ISSUES!!!
*)


(* GNSCohomologyVect::usage="Outputs a list of generators of the GNS cohomology component of fixed degree.
Using the Frobenius inner product on matrices this component is identified with a subspace of R^{d}, where d
is the dimension of the cohomology component. The output is a list of vectors whose span is this subspace.

This function takes in several possible classes of inputs:
--GNSCohomologyVect[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--GNSCohomologyVect[rho,dimprim,partition][deg]: computes the GNS cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--GNSCohomologyVect[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSCohomologyVect[rho,{2,...,2}][deg].
**'numSys' is an integer

--GNSCohomologyVect[rho,numSys,d][0]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSCohomologyVect[rho,{d,...,d}][deg].
**'d' is an integer."; 

!!!FUNCTION COMMENTED DUE TO SIGN ISSUES!!!
*)


comCohomologyObj::usage="Outputs a list of generators of the GNS cohomology component of fixed degree.
Here the elements of the list are output as functions on the set of subsystems of size deg+1 (where deg is the degree under consideration).
The value of each of these functions is a representative living in the GNS module
associated to the reduced state on that subsystem.  For instance, setting 
R=comCohomologyObj[...][deg] then R is a list of size d, where d=comCohomologyRk[...][deg] is the
dimension of that cohomology component.
The kth generator of the component (where k = 1,...,d) is R[[k]], which is a function.
Suppose deg = 2, and we are working with a 4-partite state then R[[1]][{1,2,4}] should output a matrix/
operator living in the GNS module associated to the reduced state on the subsystem {1,2,4}.  R[[1]] in this
case takes in all size 3 subsets of {1,2,3,4} (ordering of the list elements does not matter).

This function takes in several possible classes of inputs:
--comCohomologyObj[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--comCohomologyObj[rho,dimprim,partition][deg]: computes the GNS cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite commutant cohomology.
!!!Note:  When working over a coarsening associated to a partition coarser than the finest partition, the output is a list of functions that
take in subsystems of elements of the partition: e.g. if the coarsening of a four-partite state is given by
the partition {{1,3},{4},{2}} (ordering intentionally meddled with) for the purposes of demonstration),
 then, if we are looking at the degree 1 component of cohomology, the syntax R[[1]][{1,2}] denotes the assignment
to the subsystem {{1,3},4}, where the element {1,3} is now treated as a single tensor factor/ primitive subsystem. ***

--comCohomologyObj[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comCohomologyObj[rho,{2,...,2}][deg].
**'numSys' is an integer

--comCohomologyObj[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comCohomologyObj[rho,{d,...,d}][deg].
**'d' is an integer.";


GNSCohomologyObj::usage="Outputs a list of generators of the GNS cohomology component of fixed degree.
Here the elements of the list are output as functions on the set of subsystems of size deg+1 (where deg is the degree under consideration).
The value of each of these functions is a representative living in the GNS module
associated to the reduced state on that subsystem.  For instance, setting 
R=GNSCohomologyObj[...][deg] then R is a list of size d, where d=GNSCohomologyRk[...][deg] is the
dimension of that cohomology component.
The kth generator of the component (where k = 1,...,d) is R[[k]], which is a function.
Suppose deg = 2, and we are working with a 4-partite state then R[[1]][{1,2,4}] should output a matrix/
operator living in the GNS module associated to the reduced state on the subsystem {1,2,4}.  R[[1]] in this
case takes in all size 3 subsets of {1,2,3,4} (ordering of the list elements does not matter).

This function takes in several possible classes of inputs:
--GNSCohomologyObj[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--GNSCohomologyObj[rho,dimprim,partition][deg]: computes the GNS cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.
!!!Note:  When working over a coarsening associated to a partition coarser than the finest partition, the output is a list of functions that
take in subsystems of elements of the partition: e.g. if the coarsening of a four-partite state is given by
the partition {{1,3},{4},{2}} (ordering intentionally meddled with) for the purposes of demonstration),
 then, if we are looking at the degree 1 component of cohomology, the syntax R[[1]][{1,2}] denotes the assignment
to the subsystem {{1,3},4}, where the element {1,3} is now treated as a single tensor factor/ primitive subsystem. ***

--GNSCohomologyObj[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSCohomologyObj[rho,{2,...,2}][deg].
**'numSys' is an integer

--GNSCohomologyObj[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSCohomologyObj[rho,{d,...,d}][deg].
**'d' is an integer.";


GNSCohomologyGens::usage="Outputs a comlpete list of generators (i.e. a basis) of the GNS cohomology component of fixed degree.
The output is a list, each of whose elements are lists of assignments of matrices to subsystems of size deg + 1,
where 'deg' is the degree under consideration.

This function takes in several possible classes of inputs:
--GNSCohomologyGens[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--GNSCohomologyGens[rho,dimprim,partition][deg]: computes the GNS cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--GNSCohomologyGens[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSCohomologyObj[rho,{2,...,2}][deg].
**'numSys' is an integer

--GNSCohomologyGens[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSCohomologyObj[rho,{d,...,d}][deg].
**'d' is an integer.";


comCohomologyGens::usage="Outputs a comlpete list of generators (i.e. a basis) of the commutant cohomology component of fixed degree.
The output is a list, each of whose elements are lists of assignments of matrices to subsystems of size deg + 1,
where 'deg' is the degree under consideration.

This function takes in several possible classes of inputs:
--comCohomologyGens[rho,dimprim][deg]:
**'deg' is an integer: the degree of the component under consideration
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--comCohomologyGens[rho,dimprim,partition][deg]: computes the commutant cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite commutant cohomology.

--comCohomologyGens[rho,numSys][deg]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comCohomologyObj[rho,{2,...,2}][deg].
**'numSys' is an integer

--comCohomologyGens[rho,numSys,d][deg]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comCohomologyObj[rho,{d,...,d}][deg].
**'d' is an integer.";


comCohomologyRkList::usage="Outputs a list of the ranks of the commutant cohomology groups (associated to the state rho) from degree 0 
to degree N-1 (where N is the number of tensor factors) read left to right.  This function is faster than sequentially using comCohomologyRk
as it temporarily stores values from endFunctorObj and endFunctorMor.

This function takes in several possible classes of inputs:
--comCohomologyRkList[rho,dimprim]:
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--comCohomologyRkList[rho,dimprim,partition]: computes the commutant cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--comCohomologyRkList[rho,numSys]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comCohomologyRkList[rho,{2,...,2}].
**'numSys' is an integer

--comCohomologyRkList[rho,numSys,d]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comCohomologyRkList[rho,{d,...,d}].
**'d' is an integer.";


GNSCohomologyRkList::usage="Outputs a list of the ranks of the GNS cohomology groups (associated to the state rho) from degree 0 
to degree N-1 (where N is the number of tensor factors) read left to right.  This function is faster than sequentially using GNSCohomologyRk
as it temporarily stores values from GNSFunctorObj and GNSFunctorMor.

This function takes in several possible classes of inputs:
--GNSCohomologyRkList[rho,dimprim]:
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--GNSCohomologyRkList[rho,dimprim,partition]: computes the GNS cohomology using a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--GNSCohomologyRkList[rho,numSys]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSCohomologyRkList[rho,{2,...,2}].
**'numSys' is an integer

--GNSCohomologyRkList[rho,numSys,d]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSCohomologyRkList[rho,{d,...,d}].
**'d' is an integer.";


comPoly::usage="comPoly[multiState][y] outputs the commutant Poincare polynomial, in variable y, of the multipartite density state given by the multiple
arguments 'multiState'.

This function takes in several possible classes of inputs:
--comPoly[rho,dimprim][y]:
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--comPoly[rho,dimprim,partition][y]: computes the GNS Poincare polynomial associated to a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--comPoly[rho,numSys][y]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates comPoly[rho,{2,...,2}].
**'numSys' is an integer

--comPoly[rho,numSys,d][y]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates comPoly[rho,{d,...,d}][y].
**'d' is an integer.";


GNSPoly::usage="GNSPoly[multiState][y] outputs the GNS Poincare polynomial, in variable y, of the multipartite density state given by the multiple
arguments 'multiState'.

This function takes in several possible classes of inputs:
--GNSPoly[rho,dimprim][y]:
**'rho' is a positive semidefinite matrix: the density state
**'dimprim' is a list of integers (e.g. {2,2,3}): the list of dimensions of Hilbert spaces at the tensor factors

--GNSPoly[rho,dimprim,partition][y]: computes the GNS Poincare polynomial associated to a coarsening of the state specified by the
partition of tensor factors 'partition'.
**'partition' is a list of list of integers: the partition that we wish to coarsen by, e.g. {{1,3},{2}} is a partition of three tensor
factors merging together the first and the third subsystem, and the function would output the component of the
resulting bipartite GNS cohomology.

--GNSPoly[rho,numSys][y]: assumes there are numSys tensor factors with Hilbert space
dimension 2 and calculates GNSPoly[rho,{2,...,2}].
**'numSys' is an integer

--GNSPoly[rho,numSys,d][y]: assumes there are numSys tensor factors with Hilbert space
dimension d and calculates GNSPoly[rho,{d,...,d}][y].
**'d' is an integer.";


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
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-2}]
];


GNSHomologyRkList[rho_,dimprim_,cover_]:=Module[{funMor,funObj},
funObj[srcobj_]:=funObj[srcobj]=GNSFunctorObj[rho,dimprim][srcobj];
funMor[src_,tgt_]:=funMor[src,tgt]=GNSFunctorMor[rho,dimprim][src,tgt];
Table[homologyRank[deg, funObj, funMor, cover, localInProd],{deg,-1,Length@cover-2}]
];


(* ::Subsection:: *)
(*Cohomologically Graded Wrappers*)


(**Because this package was originally built to compute cech homology of a pre-cosheaf, 
the functions here are mere wrappers to correct back to cohomological grading**)


fixDegree[deg_,cover_]:=Length@cover-deg-2;


multipartiteData[rho_,dimprim_?ListQ,partition_?ListQ]:={rho,dimprim,partitionToCover@partition};

multipartiteData[rho_,dimprim_?ListQ]:={rho,dimprim,complementaryCover[Length@dimprim]};

multipartiteData[rho_,numSys_?IntegerQ,d_?IntegerQ]:={rho,ConstantArray[d,numSys],complementaryCover[numSys]};

multipartiteData[rho_,numSys_?IntegerQ]:={rho,ConstantArray[2,numSys],complementaryCover[numSys]};


cohomologicalGrading[func_][multiState__][deg_]:=Module[{multiData,cover},
multiData=multipartiteData[multiState];
cover=multiData[[3]];
(func@@multiData)[fixDegree[deg,cover]]
];


cohomologicalGradingList[func_][multiState__]:=Module[{multiData,cover},
multiData=multipartiteData[multiState];
cover=multiData[[3]];
Reverse@(func@@multiData)
];


cohomologicalGradingObj[func_][multiState__][deg_]:=Module[{multiData,cover,coSimplex,cochainGeneratorList,signFix},
multiData=multipartiteData[multiState];
cover=multiData[[3]];
coSimplex=Complement[Range[Length@cover],#]&;
cochainGeneratorList=(func@@multiData)[fixDegree[deg,cover]];
signFix[obj_][cosimp_]:=(-1)^(Total@Map[#&,cosimp]+1)*obj[cosimp];
Map[(signFix[#]@*coSimplex)&,cochainGeneratorList]
];


generators[func_][multiState__][deg_]:=Module[{numSys,gens},
numSys=Length@multipartiteData[multiState][[3]];
gens=cohomologicalGradingObj[func][multiState][deg];
Table[Map[#->MatrixForm@gens[[b]][#]&,Subsets[Range[1,numSys],{deg+1}]],{b,1,Length@gens}]
];


comCohomologyRk:=cohomologicalGrading[comHomologyRk];


GNSCohomologyRk:=cohomologicalGrading[GNSHomologyRk];


GCohomologyRk:=GNSCohomologyRk;


(* comCohomologyVect:=cohomologicalGrading[comHomologyVect]; *)


(* GNSCohomologyVect:=cohomologicalGrading[GNSHomologyVect]; *)


comCohomologyObj:=cohomologicalGradingObj[comHomologyObj];


GNSCohomologyObj:=cohomologicalGradingObj[GNSHomologyObj];


GNSCohomologyGens:=generators[GNSHomologyObj];


comCohomologyGens:=generators[comHomologyObj];


comCohomologyRkList:=cohomologicalGradingList[comHomologyRkList];


GNSCohomologyRkList:=cohomologicalGradingList[GNSHomologyRkList];


ranksToPoly[list_][y_]:=y^(Range[0,Length@list-1]).list;


comPoly:=ranksToPoly@*comCohomologyRkList;


GNSPoly:=ranksToPoly@*GNSCohomologyRkList;


(* ::Title:: *)
(*End Matter*)


End[];


EndPackage[]
