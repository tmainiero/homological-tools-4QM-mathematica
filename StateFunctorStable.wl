(* ::Package:: *)

AppendTo[$Path,NotebookDirectory[]];


BeginPackage["StateFunctorStable`",{"BasicStable`"}];


(* ::Title:: *)
(*Descriptions of Public (Unhidden) Functions*)


endFunctorObj::usage="endFunctorObj[rho,dimprim][subsys] returns a list of objects of the image of the endomorphism functor associated to the subsystem 'subsys' given the state 'rho'
and list of (primitive) system dimensions dimprim.";


endFunctorMor::usage="endFunctorMor[rho,dimprim][sys1,sys2] returns the morphism endFunctorObj[sys1] -> 
endFunctorObj[sys2] which is given by extension 'by identity operator tensoring' then compressing to the appropriate
subspace of operators.";


GNSFunctorObj::usage="GNSFunctorObj[rho,dimprim][subsys] returns a list of objects of the image of the  GNS functor associated to the subsystem 'subsys' given the state 'rho'
and list of (primitive) system dimensions dimprim.";


GNSFunctorMor::usage="GNSFunctorObj[rho,dimprim][subsys] returns a list of objects of the image of the GNS functor associated to the subsystem 'subsys' given the state 'rho'
and list of (primitive) system dimensions dimprim.";


nullFunctorObj::usage="nullFunctorObj[rho,dimprim][subsys] returns a list of objects of the image of the endomorphism functor associated to the subsystem 'subsys' given the state 'rho'
and list of (primitive) system dimensions dimprim.";


nullFunctorMor::usage="nullFunctorMor[rho,dimprim][sys1,sys2] returns the morphism endFunctorObj[sys1] -> 
nullFunctorObj[sys2] which is given by extension 'by identity operator tensoring' then compressing to the appropriate
subspace of operators.";


localInProd::usage="localInProd[subsys]:";


(* ::Title:: *)
(*Function Definitions*)


Begin["`Private`"];


(* ::Section:: *)
(*Linear Algebraic Operations*)


matrixUnit::usage ="matrixUnit[v_,w_] returns the matrix v otimes w^{vee} where w^{vee} = <w,-> using the standard metric on C^n.";

matrixUnit[v_,w_]:=KroneckerProduct[v,Conjugate@w];


homSpace::usage="homSpace[vspace1,vspace2] takes in vector spaces vspacei (i = 1,2) in the form of a list of basis vectors {v^1_1,v^1_2....} 
and outputs the vector space of homomorphisms = vspace_1 otimes (v_space_2)^{vee} from vspace1 to vspace2 in the form of a list of basis vectors v^i_1 otimes (v^j_2)^{vee}
 where (v^j_2)^{vee} is the dual basis vector of v^j_2";

homSpace[vspace1_,vspace2_]:=Flatten[#,1]&@Outer[matrixUnit,vspace2,vspace1,1];


froebenius::usage="froebenius[A,B] takes in matrices A and B and returns the Froebenius inner product Tr[ConjugateTranspose[A].B].";

froebenius[A_,B_]:=Tr[ConjugateTranspose[A].B];


eigImage::usage = "eigImage[M] returns all eigenvectors of M associated to nonzero eigenvalues.";

eigImage[M_]:= Module[{Ev=Eigensystem[M]},
Pick[Ev[[2]],(#>0)&/@Ev[[1]]]];


stdBasis::usage="stdBasis[dim] takes in an integer dim and returns the standard basis on a vector space of dimension dim in the form {v_1,v_2,...v_dim}
where each of the v_i are 1 x dim arrays.";

stdBasis[dim_]:=If[dim==0, {} , IdentityMatrix[dim]];


(* ::Section:: *)
(*Defining the Functors*)


(* ::Subsection:: *)
(*Support / Null space projections*)


powerMod[q_]:=If[#==0,0,Power[#,q]]&;
matrixPowerMod:=MatrixFunction[powerMod[#2],#1]&


suppProj[rho_]:=matrixPowerMod[rho,0];


(* nullProj[rho_]:=Module[{dimSys},
dimSys:=Dimensions[rho][[1]];
IdentityMatrix[dimSys]-suppProj[rho]
]; *)


(* ::Subsection:: *)
(*Commutant Functor*)


endFunctorObj[rho_,dimprim_]:=homSpace[#,#]&@*(eigImage@reducedDensityMat[rho,#,dimprim]&);


endFunctorMor[rho_,dimprim_][sys1_,sys2_]:=Module[{supp},
supp=suppProj@reducedDensityMat[rho,sys2,dimprim];
supp.extendOp[#,sys1,sys2,dimprim[[sys2]]].supp&
];


(* ::Subsection:: *)
(*GNS Functor*)


GNSFunctorObj[rho_,dimprim_][subsys_]:=Module[{reducedState}, 
reducedState=reducedDensityMat[rho,subsys,dimprim];
dimSys=Dimensions[reducedState][[1]];
homSpace[eigImage@reducedState,stdBasis[dimSys]]
];


GNSFunctorMor[rho_,dimprim_][sys1_,sys2_]:=extendOp[#,sys1,sys2,dimprim[[sys2]]].suppProj[reducedDensityMat[rho,sys2,dimprim]]&;


(* ::Subsection:: *)
(*Null ideal/Left kernel Functor*)


nullFunctorObj[rho_,dimprim_][subsys_]:=Module[{reducedState,hSpan}, 
reducedState=reducedDensityMat[rho,subsys,dimprim];
dimSys=Dimensions[reducedState][[1]];
hSpan=homSpace[NullSpace@reducedState,stdBasis[dimSys]];
If[hSpan==={}, {ConstantArray[0,{dimSys,1}]},hSpan]
];


nullFunctorMor[rho_,dimprim_][sys1_,sys2_]:=extendOp[#,sys1,sys2,dimprim[[sys2]]]&;


(* ::Section:: *)
(*Local inner products on objects in image of functorObj*)


localInProd[subsys_]:=froebenius[#1,#2]&;


(* ::Title:: *)
(*End Matter*)


End[];


EndPackage[]
