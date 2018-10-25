(* ::Package:: *)

AppendTo[$Path,NotebookDirectory[]];


BeginPackage["BasicStable`"];


(* ::Title:: *)
(*Descriptions of Public (Unhidden) Functions*)


(* ::Section:: *)
(*Synonyms and Sugar*)


(* ::Subsection::Closed:: *)
(*Synonyms*)


tensor::usage="Takes tensor product of an arbitrary number of vectors/matrices. 
 Utilizes 'KroneckerProduct'.";


Id::usage="Id[n] is a synonym for IdentityMatrix[n].";


IdSparse::usage="IdSparse[n] generates an n x n identity matrix as a sparse array.";


(* ::Subsection:: *)
(*Sugar*)


IdState::usage="identityState[n] generates the trace-1'identity state'/maximal von Neumann entropy
density state on n-systems: Diag(1/n,...,1/n).";


pureToDensity::usage="pureToDensity[v_] takes in a pure state 'v' and outputs its representation
as a density matrix.";


ketPure::usage="ketPure[list_] takes in an expression 'list' consisting of 0's and 1's and outputs
the standard qubit ket vector |list> as an element of C^(length of list).  Here 0 is associated to the
unit vector {1,0} and 1 is associated to {0,1}.  Adding an optional (list valued) argument 'dims', 
ketpure[list,dims] takes in an expression 'list' such that list[[k]] <= dims[[k]]-1 and outputs the 
ket vector |list>.  Here a list element 'r' is associated to the unit vector with a 1 in the rth position.";


ketPureSp::usage="ketPureSp[list] (ketPureSp[list,dims]) outputs the same result as ketPure[list] (ketPure[list,dims]) in sparse matrix form.";


ketDens::usage="See ketPure.  ketDens[list] (ketDens[list,dims]) converts the result from ketPure[list] (ketPure[list,dims])
into its corresponding density matrix representation.";


ketDensSp::usage="Sparse version of ketDens.";


ghzPure::usage="ghzPure[N] takes in an integer N and outputs the GHZ state on the N-qubit system.  ghzPure[dims]
takes in a list of dimensions 'dims' for an N=Length[dims] qudit system and is the result of
lifting the qubit vector ghzPure[N] in C^{2^N} to a qudit vector in C^{dims[[1]]} otimes C^{dims[[2]]}...
otimes C^{dims[[N]]} using the embedding C^2 ---> C^{dims[[k]]} via [1,0] |--> [1,0,0,...0] and 
[0,1] |--> [0,1,0,0....0].";


ghzDens::usage="See ghzPure.  ghzDens[N] (ghzDens[dims]) converts the result from ghzPure[N] (ghzPure[dims])
into its corresponding density matrix representation.";


ghzDensSp::usage="Sparse version of ghzDens.";


wStatePure::usage="wStatePure[N] takes in an integer N and outputs the 'W'-state on the N-qubit system. 
wStatePure[dims] takes in a list of dimensions 'dims' for an N=Length[dims] qudit system and is the result 
of lifting the qubit vector wStatePure[N] in C^{2^N} to a qudit vector in C^{dims[[1]]} otimes C^{dims[[2]]}...
otimes C^{dims[[N]]} using the embedding C^2 ---> C^{dims[[k]]} via [1,0] |--> [1,0,0,...0] and 
[0,1] |--> [0,1,0,0....0].";


wStateDens::usage="See wStatePure.  wStateDens[N] (wStateDens[dims]) converts the result from wStatePure[N] 
(wStatePure[dims]) into its corresponding density matrix representation.";


(* ::Section:: *)
(*Basic Operations*)


(* ::Subsection::Closed:: *)
(*Permute Subsystems *)


sysPermute::usage="sysPermute[rho,perm,dim] applies a permutation 'perm' to the tensor product
components of rho given the dimension vector 'dim'.";


(* ::Subsection::Closed:: *)
(*Partial Trace*)


partialTrace::usage="partialTrace[rho,systr,dim] takes the partial trace
of rho over the spaces in list 'systr' which are among those forming the tensor
product corresponding to the list dim.  E.g., let dim={2,3,4}, syslist = {1,3}; the 2
and 4 dimensional spaces are traced out to leave a 3x3 matrix.  'systr' must be a numeric list
compatible with the locations of the dimensions in 'dim'.";

reducedDensityMat::usage="reducedDensityMat[rho,subsys,dim] outputs the reduced density matrix
on the system 'subsys' given the list of dimensions 'dim' on the full system.  'subsys' must be a numeric
list compatible with the locations of the dimensions in 'dim'.  This function just applies partialTrace
to the complement of 'subsys'.";


(* ::Subsection::Closed:: *)
(*Extension of Operator by the Identity/Partial Cotrace*)


extendOp::usage ="extendOp[op,sysop,sysfull,dimfull] extends the operator 'op' on the subsystem
'sysop' to the system 'sysfull' with dimension vector 'dimfull', by tensoring by the identity
matrix/operator on the complement of sysop in sysfull.  sysop and sysfull can be lists of symbolic 
or numeric quantities.";


extendOpSparse::usage ="extendOpSparse[op,sysop,sysfull,dimfull] is the sparse version 
of extendOp: instead of tensoring by identity matrices, this function uses sparse identity 
matrices.";


partialCoTraceSparse::usage = "partialCoTraceSparse[op,syscotr,dim] is the sparse version of partialCoTrace:
instaed of tensoring by identity matrices, this function uses sparse identity matrices.  This function
refers directly to extendOpSparse.";


partialCoTrace::usage = "partialTrace[op,syscotr,dim] takes the partial cotrace of 'op'
(tensoring by identity operators in positions specified by 'syscotr') given the dimension vector
'dim' of the target.  E.g., let dim={2,3,4}, syscotr = {1,3}, and op be an
operator matrix on the second subsystem (of dimension 3); then partial cotrace returns
1_{2} otimes op otimes 1_{4}. 'syscotr' must be a numeric list compatible with the locations
of the dimensions in 'dim'.  This function refers directly to extendOp.";


(* ::Section::Closed:: *)
(*Supplementary Operations*)


(* ::Subsection::Closed:: *)
(*Partial Expectation Value*)


partialExpectation::usage="partialExpectation[rho][op,systr,dimfull] computes the partial trace
Tr_systr[rho*op] given the dimension vector 'dimfull'.";


(* ::Section:: *)
(*Entropies and Related*)


(* ::Subsection:: *)
(*von Neumann Entropy and Interaction Information*)


vonNeumann::usage="vonNeumann[rho] calculates the von Neumann entropy of a density matrix rho.";

interactionInfo::usage="interactionInfo[rho,primsys,dimprim] calculates the interaction information of a
density matrix rho defined on the primitive subsystems 'primsys' with dimension vectors 'dimprim'.  Here
interaction information is defined via the alternating sum Sum[(-1)^(N-|lambda|-1) S_(von Neumann)(rho_lambda)]
, where the sum runs over all subsets lambda of primitive subystems, N is the number of primitive subystems, and
 rho_lambda is the reduced density matrix on lambda.";


(* ::Subsection:: *)
(*Euler Characteristics of GNS and Commutant Complexes*)


eulerCharG::usage="eulerCharG[rho_,primsys_,dimprim_] computes the Euler characteristic of the GNS complex of a multipartite density state rho
defined on the primitive subsystems 'primsys' with dimension vectosr 'dimprim'.";


eulerCharGNS::usage="synonym for eulerCharG'.


eulerCharE::usage="eulerCharG[rho_,primsys_,dimprim_] computes the Euler characteristic of the commutant complex of a multipartite density state rho
defined on the primitive subsystems 'primsys' with dimension vectosr 'dimprim'.";


eulerCharCom::usage="synonym for eulerCharE'.


(* ::Subsection:: *)
(*q-deformed quantities*)


(* ::Subsubsection:: *)
(*Tsallis Entropy, the q-deformed (Tsallis) Interaction information, and related quantities*)


tsallis::usage="tsallis[rho,q] calculates the (q-deformed) tsallis entropy of a density matrix rho.";

qinteractionInfo::usage="qinteractionInfo[rho,primsys,dimprim] calculates the q-deformed interaction information
 of a density matrix rho defined on the primitive subsystems 'primsys' with dimension vectors 'dimprim'.  Here
the q-deformed interaction information is defeind via the alternating sum Sum[(-1)^(N-|lambda|-1) S_(tsallis)(rho_lambda,q)]
, where the sum runs over all subsets lambda of primitive subystems, N is the number of primitive subystems, and
 rho_lambda is the reduced density matrix on lambda.";

qPartitionFunc::usage="qPartitionFunc[rho,q] calculates Tr[rho^q]";

qEulerChar::usage="qEulerChar[rho,primsys,dimprim][q] calculates the alternating sum Sum[(-1)^(N-|lambda|-1) Tr[rho^q],
 where the sum runs over all subsets lambda of primitive subystems, N is the number of primitive
 subystems, and rho_lambda is the reduced density matrix on lambda.  
(q-1)*qEulerChar[rho,primsys,dimprim,q] is equal to qinteractionInfo[rho,primsys,dimprim,q].";


qrEulerChar::usage="qrEulerChar[rho,primsys,dimprim][q,r] calculates the two parameter alternating sum Sum_{lambda}(-1)^(N-|lambda|-1) Tr[rho^q]^r,
 where the sum runs over all subsets lambda of primitive subystems, N is the number of primitive
 subystems, and rho_lambda is the reduced density matrix on lambda.  The specialization at r=1 gives qEulerChar.";


stateIndex::usage="qrEulerChar[rho,primsys,dimprim][q,r] calculates the two parameter alternating sum Sum_{lambda}(-1)^(N-|lambda|-1) Tr[rho^q]^r,
 where the sum runs over all subsets lambda of primitive subystems, N is the number of primitive
 subystems, and rho_lambda is the reduced density matrix on lambda.  The specialization at r=1 gives qEulerChar.";


(* ::Subsubsection:: *)
(*Renyi Entropy*)


renyi::usage="renyi[rho,alpha] calculates the renyi alpha-entropy of a density matrix rho.";


(* ::Title:: *)
(*Function Definitions*)


Begin["`Private`"];


(* ::Section:: *)
(*Synonyms and Sugar (Useful Definitions)*)


(* ::Subsection::Closed:: *)
(*Synonyms*)


tensor[args__?MatrixQ]:=KroneckerProduct[args];


tensor[args__?VectorQ]:=Flatten@KroneckerProduct[args];


Id:=IdentityMatrix;


IdSparse[n_]:= SparseArray[{i_,i_}->1,{n,n}];


(* ::Subsection::Closed:: *)
(*Sugar*)


IdState[N_]:=1/N*IdentityMatrix[N];


IdStateSparse[n_]:=SparseArray[{i_,i_}->1/n,{n,n}];


(*pureToDensity[v_?VectorQ]:=KroneckerProduct[{v},ConjugateTranspose[{v}]] *)


pureToDensity[v_?VectorQ]:=Transpose@KroneckerProduct[v,Conjugate@v];


ketPure[list_,dims_]:=If[list==={}&&dims==={},{1},
(*else*)
tensor[UnitVector[First@dims,First@list+1],#]&@ketPure[Rest@list,Rest@dims]
];

ketPure[list_]:=ketPure[list,ConstantArray[2,Length@list]];


ketPureSp[list_,dims_]:=Module[{sparseUnitVect},

spUnitVect[dimension_,k_]:=SparseArray[{k}->1,{dimension}];

If[list==={}&&dims==={},spUnitVect[1,1],
(*else*)
tensor[spUnitVect[First@dims,First@list+1],#]&@ketPureSp[Rest@list,Rest@dims]
]
];

ketPureSp[list_]:=ketPureSp[list,ConstantArray[2,Length@list]];


ketDens:=pureToDensity@*ketPure;


ketDensSp:=pureToDensity@*ketPureSp;


ghzPure[dims_?VectorQ]:=Module[{N=Length@dims},
(1/Sqrt[2])*(ketPure[ConstantArray[0,N],dims]+ketPure[ConstantArray[1,N],dims])
];

ghzPure[N_?IntegerQ]:=ghzPure[ConstantArray[2,N]];


ghzPureSp[dims_?VectorQ]:=Module[{N=Length@dims},
(1/Sqrt[2])*(ketPureSp[ConstantArray[0,N],dims]+ketPureSp[ConstantArray[1,N],dims])
];

ghzPureSp[N_?IntegerQ]:=ghzPureSp[ConstantArray[2,N]];


ghzDens:=pureToDensity@*ghzPure;


ghzDensSp:=pureToDensity@*ghzPureSp;


wStatePure[dims_?VectorQ]:=Module[{N=Length@dims},
(1/Sqrt[N])*Sum[ketPure[UnitVector[N,k],dims],{k,1,N}]
];

wStatePure[N_?IntegerQ]:=wStatePure[ConstantArray[2,N]];


wStateDens:=pureToDensity@*wStatePure;


(* ::Section:: *)
(*Basic Operations*)


(* ::Subsection::Closed:: *)
(*Permute Subsystems*)


sysPermute[rho_,perm_,dim_]:=Module[{permTensor},
permTensor=Join[perm,Length@dim+perm];
ArrayReshape[#,Dimensions@rho]&@Transpose[#,permTensor]&@ArrayReshape[rho,Join[dim,dim]]
];


(* ::Subsection::Closed:: *)
(*Partial Trace and Reduced Density Matrix*)


partialTrace[rho_,systr_,dim_]:=Module[{keep,permTensor,rhoint,dimKeep,dimTrace},
keep=Complement[Range@Length@dim,systr];
dimTrace=Times@@dim[[systr]];
dimKeep=Times@@dim/(dimTrace);
permTensor=InversePermutation@Join[keep,Length@dim+keep,systr,Length@dim+systr];
rhoint=ArrayReshape[#,{dimKeep,dimKeep,dimTrace^2}]&@Transpose[#,permTensor]&@ArrayReshape[rho,Join[dim,dim]];
Sum[rhoint[[All,All,k]],{k,1,dimTrace^2,dimTrace+1}]
];


reducedDensityMat[rho_,subsys_,dim_]:=partialTrace[rho,Complement[Range@Length[dim],subsys],dim];


(* ::Subsection::Closed:: *)
(*Extension of Operator by the Identity/Partial Cotrace*)


extendOp[op_,sysop_,sysfull_,dimfull_]:=
If[sysop==={},Id[Times@@dimfull],
(*else*) Module[{sysExtend,dimCurrent,dimExtend,sysExtendPos,sysOpPos,currentPerm},
sysExtend=Complement[sysfull,sysop];
sysExtendPos=Flatten@Map[Position[sysfull,#]&,sysExtend];
sysOpPos=Complement[Range@Length[sysfull],sysExtendPos];
currentPerm=Join[sysOpPos,sysExtendPos];
dimCurrent=dimfull[[currentPerm]];
dimExtend=Times@@dimfull[[sysExtendPos]];
sysPermute[tensor[op,Id[dimExtend]], currentPerm, dimCurrent] ]
];


extendOpSparse[op_,sysop_,sysfull_,dimfull_]:=
If[sysop==={},Id[Times@@dimfull],
(*else*) Module[{sysExtend,dimCurrent,dimExtend,sysExtendPos,sysOpPos,currentPerm},
sysExtend=Complement[sysfull,sysop];
sysExtendPos=Flatten@Map[Position[sysfull,#]&,sysExtend];
sysOpPos=Complement[Range@Length[sysfull],sysExtendPos];
currentPerm=Join[sysOpPos,sysExtendPos];
dimCurrent=dimfull[[currentPerm]];
dimExtend=Times@@dimfull[[sysExtendPos]];
sysPermute[tensor[op,IdSparse[dimExtend]], currentPerm, dimCurrent] ]
];


partialCoTrace[op_,syscotr_,dim_]:=Module[{sysop,sysfull},
sysfull=Range@Length@dim;
sysop=Complement[sysfull,syscotr];
extendOp[op,sysop,sysfull,dim]
];


partialCoTraceSparse[op_,syscotr_,dim_]:=Module[{sysop,sysfull},
sysfull=Range@Length@dim;
sysop=Complement[sysfull,syscotr];
extendOpSparse[op,sysop,sysfull,dim]
];


(* ::Section:: *)
(*Supplementary Operations*)


(* ::Subsection::Closed:: *)
(*Join (tensor) States/Operators on disjoint/independent subsystems*)


joinAndOrder[A_,B_,sysA_,sysB_,dimA_,dimB_]:=Module[{perm},
perm=InversePermutation@Ordering@Join[sysA,sysB]
sysPermute[tensor[A,B],perm,Join[dimA,dimB]]
];


(* ::Subsection::Closed:: *)
(*Partial Expectation Value*)


partialExpectation[rho_][op_,systr_,dimfull_]:=partialTrace[rho*op,systr,dimfull];


(* ::Subsection::Closed:: *)
(*Extension of state by specified state on independent system(s)/Co-partial expectation value*)


extendStateBy[newState_][rho_,sysrho_,sysfull_,dimfull_]:=
If[sysrho==={},newState,
(*else*) Module[{sysExtend,dimCurrent,dimExtend,sysExtendPos,sysrhoPos,currentPerm},
sysExtend=Complement[sysfull,sysrho];
sysExtendPos=Flatten@Map[Position[sysfull,#]&,sysExtend];
sysrhoPos=Complement[Range@Length[sysfull],sysExtendPos];
currentPerm=Join[sysrhoPos,sysExtendPos];
dimCurrent=dimfull[[currentPerm]];
dimExtend=Times@@dimfull[[sysExtendPos]];
sysPermute[tensor[rho,newState], currentPerm, dimCurrent] ]
];


coPartialExpectation[rho_][chi_,sysextend_,dimfull_]:=Module[{sysfull,rhoExtend,syschi},
sysfull=Range@Length@dimfull;
rhoExtend=reducedDensityMat[rho,sysextend];
syschi=Complement[sysfull,sysextend];
extendStateBy[rho][chi,syschi,sysfull,dimfull]
];


(* ::Subsection::Closed:: *)
(*Extension of Density State by Identity States *)


extendState[rho_,sysstate_,sysfull_,dimfull_]:=
If[sysstate==={},IdState[Times@@dimfull],
(*else*) Module[{sysExtend,dimCurrent,dimExtend,sysExtendPos,sysOpPos,currentPerm},
sysExtend=Complement[sysfull,sysstate];
sysExtendPos=Flatten@Map[Position[sysfull,#]&,sysExtend];
sysOpPos=Complement[Range@Length[sysfull],sysExtendPos];
currentPerm=Join[sysOpPos,sysExtendPos];
dimCurrent=dimfull[[currentPerm]];
dimExtend=Times@@dimfull[[sysExtendPos]];
sysPermute[tensor[rho,IdState[dimExtend]], currentPerm, dimCurrent] ]
];


extendStateSparse[rho_,sysstate_,sysfull_,dimfull_]:=
If[sysstate==={},IdStateSparse[Times@@dimfull],
(*else*) Module[{sysExtend,dimCurrent,dimExtend,sysExtendPos,sysOpPos,currentPerm},
sysExtend=Complement[sysfull,sysstate];
sysExtendPos=Flatten@Map[Position[sysfull,#]&,sysExtend];
sysOpPos=Complement[Range@Length[sysfull],sysExtendPos];
currentPerm=Join[sysOpPos,sysExtendPos];
dimCurrent=dimfull[[currentPerm]];
dimExtend=Times@@dimfull[[sysExtendPos]];
sysPermute[tensor[rho,IdStateSparse[dimExtend]], currentPerm, dimCurrent] ]
];


(* ::Subsection::Closed:: *)
(* Adjoints under State-Operator Pairing*)


stdMatrix[i_,j_,size_]:=SparseArray[{i,j}->1,size];


mapofOpsToMat[map_]:=


genAdjoint[M_,inprodsrc_,inprodtgt_]:= Module[{ipMatSrc, ipMatTgt},
ipMatSrc=inprodToMat[inprodsrc,IdentityMatrix[Dimensions[M][[1]]]];
ipMatTgt=inprodToMat[inprodtgt,IdentityMatrix[Dimensions[M][[2]]]];
PseudoInverse[ipMatSrc].Transpose[M].ipMatTgt
];


stateMapToOperatorMap[statemap_]:=


















(* ::Section:: *)
(*Entropies and Related*)


(* ::Subsection:: *)
(*Inclusion-Exclusion Sum/Alternating Sum of Evaluations on Reduced Density Matrices*)


index[rho_,primsys_,dimprim_,fun_]:=Module[{N=Length@primsys},
Total@Map[(-1)^Length@#*fun@reducedDensityMat[rho,#,dimprim]&,Subsets[primsys]]
];


shiftedIndex[rho_,primsys_,dimprim_,fun_]:=(-1)*index[rho,primsys,dimprim,fun];


(* ::Subsection:: *)
(*von Neumann Entropy and Interaction Information*)


vNKernel[lambda_]:=If[lambda==0,0,-lambda*Log[lambda]];

vNMatrixKernel:=MatrixFunction[vNKernel,#]&;

vonNeumann[rho_?PositiveSemidefiniteMatrixQ]:=Tr[vNMatrixKernel@rho];


interactionInfo[rho_,primsys_,dimprim_]:=shiftedIndex[rho,primsys,dimprim,vonNeumann];


(* ::Subsection:: *)
(*Euler Characteristics of GNS and Commutant Complexes*)


GNSDim[rho_?SquareMatrixQ]:=(Dimensions[rho][[1]])*MatrixRank[rho];


ComDim:=MatrixRank[#]^2&;


eulerCharG[rho_,primsys_,dimprim_]:=shiftedIndex[rho,primsys,dimprim,GNSDim];


eulerCharE[rho_,primsys_,dimprim_]:=shiftedIndex[rho,primsys,dimprim,ComDim];


eulerCharGNS:=eulerCharG;


eulerCharCom:=eulerCharE;


(* ::Subsection:: *)
(*q-deformed quantities*)


powerMod[q_]:=If[#==0,0,Power[#,q]]&;
matrixPowerMod:=MatrixFunction[powerMod[#2],#1]&

qLog[rho_?NumberQ,q_]:=If[q===1,Log[rho],
(*else*)
1/(1-q)*(powerMod[1-q][rho]-1)
];

qLog[rho_?PositiveSemidefiniteMatrixQ,q_]:=If[q===1,Log[rho],
(*else*)
1/(1-q)*(matrixPowerMod[rho,1-q]-IdentityMatrix@Dimensions[rho])
];


(* ::Subsubsection:: *)
(*Tsallis Entropy, the q-deformed Interaction information, and related quantities*)


tsallis[rho_?PositiveSemidefiniteMatrixQ, q_]:=If[q===1,vonNeumann[rho],
(*else*)
-Tr[rho.qLog[rho,q]]
];


qInteractionInfo[rho_,primsys_,dimprim_,q_]:=shiftedIndex[rho,primsys,dimprim,tsallis[#,q]&];


qPartitionFunc[rho_,q_]:=Tr[matrixPowerMod[rho,q]];


finStateDim[rho_][q_,r_,\[Alpha]_]:=Module[{hilbDim},
hilbDim[M_]:=Dimensions[M][[1]];
hilbDim[rho]^\[Alpha]*qPartitionFunc[rho,q]^r
];


qEulerChar[rho_,primsys_,dimprim_][q_]:=shiftedIndex[rho,primsys,dimprim,qPartitionFunc[#,q]&];


qrEulerChar[rho_,primsys_,dimprim_][q_,r_]:=shiftedIndex[rho,primsys,dimprim,qPartitionFunc[#,q]^r&];


stateIndex[rho_,primsys_,dimprim_][q_,r_,\[Alpha]_,w_]:=w^(Length@primsys)*index[rho,primsys,dimprim, finStateDim[#][q,r,\[Alpha]]&];


(* ::Subsubsection:: *)
(*Renyi Entropy*)


renyi[rho_,q_]:=If[q===1,vonNeumann[rho],1/(1-q)*Log@qPartitionFunc[rho,q]];


(* ::Title:: *)
(*End Matter*)


End[];


EndPackage[]
