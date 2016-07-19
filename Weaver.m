(* ::Package:: *)

BeginPackage["Weaver`"]

Print["Weaver: a package for explicit exact-arithmetic calculations in W-algebras"];


Clear["Weaver`*"];
Clear["Weaver``Private`*"];


(* User Defined Symbols *)
OperatorList::usage="This symbol is a list of names of operators.";
NullDescendents::usage="NullDescendents[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] find the subspace at \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\) that are descendents of null vectors at lower levels.";
FermionQ::usage="FermionQ[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)] returns True if \!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\) is a fermion and False otherwise.";
\[CapitalDelta]::usage="\[CapitalDelta][\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)] returns the conformal weight of \!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)";
Grading::usage="Grading[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)] returns the grading of \!\(\*
StyleBox[\"operator\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"name\",\nFontSlant->\"Italic\"]\), usually 0 for bosons or
R-sector fermions and 1/2 for NS-sector fermions.";
Commutator::usage="A relation defining the (anti)commutation relations of the theory.  For every operator X,Y, Commutator
should be defined using the following format:\n
Commutator[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n_\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"m_\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level_\",\nFontSlant->\"Italic\"]\)]:=(n-m)Operator[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"+\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)]+c/12(n^3-n)\[Delta][\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"+\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)]\n
using Operator[\!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] for Subscript[X, l], and NOperator[X,Y,l,level] for :XYSubscript[:, l].";

MinimalGeneratorsCreation::"usuage"="MinimalGeneratorsCreation should be set to be a list of operators that
generate the negative level subalgebra of the algebra.  Each operators is written in the form {operator_name,level}.";
MinimalGeneratorsAnnihilation::"usuage"="MinimalGeneratorsAnnihilation should be set to be a list of operators that
generate the positive level subalgebra of the algebra.  Each operators is written in the form {operator_name,level}.";

HighestWeightDimension::usage="HighestWeightDimension should be set to the dimension of the highest weight states.";

(* Fundamental algebra manipulations *)
FSign::usage="FSign[\!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\)] returns -1 if both \!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\) are fermions and 1 otherwise";

AddOperator::usage="AddOperator[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"fermionness\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"grading\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"conformal_weight\",\nFontSlant->\"Italic\"]\)] initializes a new 
operator named \!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\), which follows fermionic statistics is \!\(\*
StyleBox[\"fermionness\",\nFontSlant->\"Italic\"]\) is set to True, and with 
levels that are in (\!\(\*
StyleBox[\"grading\",\nFontSlant->\"Italic\"]\) + Z).  The paramater \!\(\*
StyleBox[\"conformal_weight\",\nFontSlant->\"Italic\"]\) is optional, and will set the conformal
weight of the operator. This is neccessary is normal-ordered proucts exist in the commutation relations.
AddOperator is equivalent to adding operator_name to OperatorList, defining FermionQ, \[CapitalDelta], and Grading on 
operator_name, and defining the Hermetian conjugate: HC[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\),n,\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)]";
Operator::usage="Operator[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"operator_level\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a matrix representing the action
of \!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"operator_level\",\nFontSlant->\"Italic\"]\)\)]\) acting on states in the Verma module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
NOperator::usage="NOperator[\!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"op_level\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] return a matrix representing the action of the normal-ordered product,
:\!\(\*
StyleBox[\"XY\",\nFontSlant->\"Italic\"]\)\!\(\*SubscriptBox[\(:\), \(\!\(\*
StyleBox[\"op_level\",\nFontSlant->\"Italic\"]\)\)]\)acting on the states in the Verma module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).
NOperator[\!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"Y\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"op_level\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the normal ordered product :\!\(\*SuperscriptBox[\(\[PartialD]\), \(\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\)]\)\!\(\*
StyleBox[\"X\",\nFontSlant->\"Italic\"]\)\!\(\*SuperscriptBox[\(\[PartialD]\), \(\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\)]\)Y\!\(\*SubscriptBox[\(:\), \(\!\(\*
StyleBox[\"op_level\",\nFontSlant->\"Italic\"]\)\)]\).";
EvaluatedOperators::usage="EvaluatedOperators is a list of all the operators that have previously been evaluated.";
\[Delta]::usage="\[Delta][\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the 0 operator is \!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)!=0 and the identity operator if \!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)=0 acting on states in the Verma
module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
OperatorNames::usage="OperatorNames[\!\(\*
StyleBox[\"-\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a list of all of the operators at a given negative level that are in the universal covering algebra
of the negative level operators.  Each operator is given as a product of
operators in OperatorList, listing those operators by {operator_name,operator_level}.";
PositiveOperatorNames::usage="PositiveOperatorNames[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a list of all of the operators at a given positive level that are in the universal covering algebra
of the positive level operators.  Each operator is given as a product of
operators in OperatorList, listing those operators by {operator_name,operator_level}.";
StateNames::usage="StateNames[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a list of states in the Verma module at a given level.  The states are products of operators
times one of the highest-weight vectors.";
LLen::usage="LLen[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the dimension of the subspace of the Verma module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
StateVector::usage="StateVector[\!\(\*
StyleBox[\"name\",\nFontSlant->\"Italic\"]\)] returns the vector corresponding to a state with a given name in the format
given by StateNames.";
HC::usage="HC[\!\(\*
StyleBox[\"operator_name\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a matrix representing the Hermetian conjugate acting on the \!\(\*
StyleBox[\"right\",\nFontWeight->\"Bold\"]\)
on the Verma module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)";
OperatorChain::usage="OperatorChain[\!\(\*
StyleBox[\"operator_list\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the matrix corresponding to the the product of
the operators given in operator_list.  The format of operator_list is as the output of OperatorNames.";
HCOperatorChain::usage="HCOperatorChain[\!\(\*
StyleBox[\"operator_list\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] gives the Hermitian conjugate of the product of
operators in operator_list, acting on the right on the subspace of the Verma module at level level.";
DisplayName::usage="DisplayName[\!\(\*
StyleBox[\"state_name\",\nFontSlant->\"Italic\"]\)] returns a more human-readable form of an operator or state name."

(* Functions for interpretation of results *)
GramMatrix::usage::"Gram Matrix[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the Gram Matrix (the inner product for a unitary representation)
 at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).  GramMatrix assumes that the highest weight subspace is already decomposed into an orthonormal basis.";
KacDeterminant::usage="KacDeterminant[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the Kac determinant at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
GramNullspace::usage="GramNullspace[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a basis for the null vectors at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
GramRank::usage="GramRank[\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns the rank of the Gram matrix at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).";
MinimalRowReducedForm::usage="MinimalRowReducedForm[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] returns the row reduced form of \!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\) with the zero rows removed.";
Descendents::usage="Descendents[\!\(\*
StyleBox[\"source\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"source_level\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"target_level\",\nFontSlant->\"Italic\"]\)] finds a basis for the descendents of the vector, \!\(\*
StyleBox[\"source\",\nFontSlant->\"Italic\"]\), at
level \!\(\*
StyleBox[\"target_level\",\nFontSlant->\"Italic\"]\) via the operators in MinimalGeneratorsCreation.";
DescendentsOfSpace::usage="Descendents[\!\(\*
StyleBox[\"source_basis\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"source_level\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"target_level\",\nFontSlant->\"Italic\"]\)] finds a basis for the descendents of the subspace,
 \!\(\*
StyleBox[\"source_basis\",\nFontSlant->\"Italic\"]\), at level \!\(\*
StyleBox[\"target_level\",\nFontSlant->\"Italic\"]\) via the operators in MinimalGeneratorsCreation.";
(*FindEigenspace::usage="FindEigenspace[X,val,basis,level] finds the subspace of the space generated by basis at level level with \!\(\*SubscriptBox[\(X\), \(0\)]\) eigenvalue of val.
If basis is omitted, FindEigenspace assumes the space to be the entire subspace of the Verma module at that level";*)
FindEigenbasis::usage="FindEigenbasis[\!\(\*
StyleBox[\"Xs\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a list of simultanious eigenvectors and eigenvalues of the operator \!\(\*FormBox[SubscriptBox[\(X\), \(\(0\)\(\\\ \)\)],
TraditionalForm]\)for X in \!\(\*
StyleBox[\"Xs\",\nFontSlant->\"Italic\"]\) acting on the subspace of the Verma module at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\) generated by \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\).  If \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\) is omitted, FindEigenbasis uses the entire space at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).  FindEigenvalue may return unexpected results if \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\) is not the union of eigenspaces of \!\(\*FormBox[SubscriptBox[\(X\), \(0\)],
TraditionalForm]\) or if the operators in \!\(\*
StyleBox[\"Xs\",\nFontSlant->\"Italic\"]\) do not commute.";
IntersectVectorSpaces::usage="IntersectVectorSpaces[\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"1\",\nFontSlant->\"Italic\"]\)\)]\),\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"2\",\nFontSlant->\"Italic\"]\)\)]\),...,\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\)]\)] returns a basis for the intersection of the spaces
generated by the \!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)\)]\)s.";
UnionVectorSpaces::usage="IntersectVectorSpaces[\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"1\",\nFontSlant->\"Italic\"]\)\)]\),\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"2\",\nFontSlant->\"Italic\"]\)\)]\),...,\!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\)]\)] returns a basis for the union of the spaces
generated by the \!\(\*SubscriptBox[\(\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)\), \(\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)\)]\)s.";
EMatrixRank::usage="EMatrixRank[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] returns MatrixRank[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] if \!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\) has at least one row and returns 0 if \!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)={}.";
FindSingularSubspace::usage="FindSingularSubspace[\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)] returns a basis for the set of singular vectors in the span of \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\)
at level \!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\).  If \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\) is omitted, FindSingularSubspace acts over the nullspace at that level.";
InVectorSpaceQ::usage="InVectorSpaceQ[\!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"vector\",\nFontSlant->\"Italic\"]\)] returns True is \!\(\*
StyleBox[\"vector\",\nFontSlant->\"Italic\"]\) is in the span of \!\(\*
StyleBox[\"basis\",\nFontSlant->\"Italic\"]\) and False otherwise.";
DescendentQ::usage="DescendentQ[\!\(\*
StyleBox[\"source\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"sourcelevel\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"target\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"targetlevel\",\nFontSlant->\"Italic\"]\)] returns True if \!\(\*
StyleBox[\"target\",\nFontSlant->\"Italic\"]\) is a descendent of \!\(\*
StyleBox[\"source\",\nFontSlant->\"Italic\"]\).";
ConstructLadderDiagram::usage="ConstructLadderDiagram[\!\(\*
StyleBox[\"cartan\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"max_level\",\nFontSlant->\"Italic\"]\)]={singulars,parentage} returns a list of singular vectors
(including the highest weight vector) up to level \!\(\*
StyleBox[\"max_level\",\nFontSlant->\"Italic\"]\) labeled under eigenvalues of the operators listing in \!\(\*
StyleBox[\"cartan\",\nFontSlant->\"Italic\"]\).
parentage is the adjacency matrix of the corresponding embedding diagram: parentage[[i,j]]=1 if j is a descendent of i.";
LadderGraph::usage="LadderGraph[\!\(\*
StyleBox[\"singulars\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"parentage\",\nFontSlant->\"Italic\"]\)] displays the embedding diagram calculated by ConstructLadderDiagram.  Call it as
LadderGraph@@ConstructLadderDiagram[\!\(\*
StyleBox[\"cartan\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"level\",\nFontSlant->\"Italic\"]\)].";
FindAnnihilatingOperators::usage="FindAnnihilatingOperators[\!\(\*
StyleBox[\"operator_level\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"state\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"state_level\",\nFontSlant->\"Italic\"]\)] finds a basis for the space of operators
at level \!\(\*
StyleBox[\"operator_level\",\nFontSlant->\"Italic\"]\) that annihilate \!\(\*
StyleBox[\"state\",\nFontSlant->\"Italic\"]\).";


Begin["`Private`"]


If[!ListQ[OperatorList],OperatorList={}];
EvaluatedOperators={};

OperatorOrder[op_]:=(For[i=1,i<=Length[OperatorList],i++,OperatorOrder[OperatorList[[i]]]=i;];OperatorOrder[op])
FSign[x_,y_]:=(-1)^Boole[FermionQ[x]&&FermionQ[y]]
\[Delta][n_,level_]:=Which[n==0,IdentityMatrix[LLen[level]],n!=0,Table[0,{i,1,LLen[level-n]},{j,1,LLen[level]}]]


DisplayName[name_]:=Module[{},
If[name=={},Return["Id"];];
If[name[[1]]=={},Return["|"<>ToString[name[[2]],InputForm]<>">"];];

If[ListQ[name[[1,1]]],
	ops = name[[1]];
	appendstring = "|"<>ToString[name[[2]]]<>">";
,
	ops = name;
	appendstring = "";
];
out = "";
Do[
out = out<>"\!\("<>ToString[op[[1]]]<>"\_\("<>ToString[op[[2]],InputForm]<>"\)\)"
,{op,ops}];
Return[out<>appendstring];
]

AddOperator[name_,fermion:True|False,grading_,conformalweight_]:=Module[{},
	If[!ListQ[OperatorList],OperatorList={}];
	AppendTo[OperatorList,name];
	OperatorList=DeleteDuplicates[OperatorList];
	FermionQ[name]=fermion;
	Grading[name]=grading;
	\[CapitalDelta][name]=conformalweight;
	HC[name,n_,level_]:= Operator[name,-n,level-n]
]
Operator[X_,n_,level_]:=Operator[X,n,level]=
	(AppendTo[EvaluatedOperators,{X,n,level}];Which[
	n==0&&level==0&&MatrixQ[X[0]],X[0],
	n==0&&level==0,{{X[0]}},
	n>level,Table[0,{i,1,0},{j,1,LLen[level]}],
	True, Transpose@Table[ApplyOperator[X,n,level,state],{state,StateNames[level]}]
])

ApplyOperator[X_,n_,level_,state_]:= (* Applies the operator, calculating recursively *)
Which[state[[1]]=={},StateVector[{{{X,n}},state[[2]]},n],True,
Module[{A,m,Bs,out1,out2},
(*Print["Evaluating ApplyOperator[",X,", ",n,", ",level,", ",state,"]"];*)
{A,m} = state[[1,1]];
Bs=state[[1,2;;]];
If[n<m,Return@StateVector[{Prepend[state[[1]],{X,n}],state[[2]]},level+n];];
If[n==m&&X==A&&FermionQ[X]==False,Return@StateVector[{Prepend[state[[1]],{X,n}],state[[2]]},level+n];];
If[n==m&&X==A&&FermionQ[X]==True,Return@(1/2 Commutator[X,n,X,n,level+m].StateVector[{Bs,state[[2]]}]);];
If[n==m&&OperatorOrder[X]<OperatorOrder[A],Return@StateVector[{Prepend[state[[1]],{X,n}],state[[2]]},level+n];];
If[n>level,Return@{};];
If[n>level+m,Return[Commutator[X,n,A,m,level+m].StateVector[{Bs,state[[2]]}]];];
(*Print[X,n,A,m,level+m];*)
out1 = Commutator[X,n,A,m,level+m].StateVector[{Bs,state[[2]]}];
out2 = FSign[X,A]Operator[A,m,level+m-n].(Operator[X,n,level+m].StateVector[{Bs,state[[2]]}]);
Return@(out1+out2);
]]

nopcoefficient[X_,i_,n_]:=Product[x,{x,-n-\[CapitalDelta][X]-i+1,-n-\[CapitalDelta][X]}]
NOperator[X_,Y_,n_,level_]:=(Sum[(*Print["NOp ",X,p,Y,n-p];*)Reorder[X,p,Y,n-p,level],{p,prevorequal[-\[CapitalDelta][X],X],n-level,-1}]+FSign[X,Y]Sum[(*Print["NOp ",Y,n-p,X,p];*)Reorder[Y,n-p,X,p,level],{p,next[-\[CapitalDelta][X],X],level}])
NOperator[X_,x_,Y_,y_,n_,level_]:=
(Sum[nopcoefficient[X,x,p]nopcoefficient[Y,y,n-p]Reorder[X,p,Y,n-p,level],{p,prevorequal[-\[CapitalDelta][X],X],n-level,-1}]
+FSign[X,Y]Sum[nopcoefficient[X,x,p]nopcoefficient[Y,y,n-p]Reorder[Y,n-p,X,p,level],{p,next[-\[CapitalDelta][X],X],level}])

next[del_,X_]:=Module[{grading,y},grading = Grading[X];y = Mod[del,1,grading];Return[del+1+grading-y];]
prevorequal[del_,X_]:=Module[{grading,y},grading = Grading[X];y = Mod[del,1,grading];Return[del+grading-y];]

Reorder[X_,n_?NumericQ,Y_,m_?NumericQ,level_]:=Which[
(* Keep the ordering *)n<m||(n==m &&OperatorOrder[X]<OperatorOrder[Y]),Operator[X,n,level-m].Operator[Y,m,level],
(* Reverse and preserve both terms *)n<=level,FSign[X,Y]Operator[Y,m,level-n].Operator[X,n,level]+Commutator[X,n,Y,m,level],
(* Reverse, but remove a term *)True,Commutator[X,n,Y,m,level]
]

HC[name_]:=HC[name,0]
HC[{},level_]:=IdentityMatrix[LLen[level]]
OperatorChain[list_,level_]:=OperatorChain[list,level]=Which[Length[list]==0,IdentityMatrix[LLen[level]],True,
Operator[list[[1,1]],list[[1,2]],level-Total[list[[2;;,2]]]].OperatorChain[list[[2;;]],level]
]

(* Acts on the right of level *)
HCOperatorChain[list_,level_]:=HCOperatorChain[list,level]=Which[Length[list]==0,IdentityMatrix[LLen[level]],True,
HCOperatorChain[list[[2;;]],level].HC[list[[1,1]],list[[1,2]],level-Total[list[[2;;,2]]]]
]


StateNames[level_] := Flatten[Table[{i,j},{i,OperatorNames[level]},{j,1,HighestWeightDimension}],1]
negativeify[op_]:=({#[[1]],-#[[2]]})&/@op
PositiveOperatorNames[level_]:=negativeify/@OperatorNames[level]
OperatorNames[level_]:=OperatorNames[level,level,1]
OperatorNames[level_,maxlevel_,index_]:=Module[{gcd,newmaxlevel,newindex},
If[level<0,Return[{}]];
	If[level==0,Return@{{}}];
If[maxlevel<=0,Return[{}]];

	gcd = GCD@@(Grading/@OperatorList);If[gcd==0,gcd=1;];
	newmaxlevel = maxlevel;
	newindex = index; If[newindex>Length[OperatorList],newmaxlevel = newmaxlevel-gcd; newindex=1;];
	If[level<newmaxlevel,newmaxlevel=level;newindex=1;];

	While[Mod[Grading[OperatorList[[newindex]]]-newmaxlevel,1]!=0,
		If[newindex==Length[OperatorList],newindex=1;newmaxlevel-=gcd;,newindex++;];
		If[newmaxlevel<=0,Return[{}]];
	];
	(* Apply the first operator we find *)
	If[FermionQ[OperatorList[[newindex]]],
		(* Fermion *)
		Return@Join[
			(Prepend[#,{OperatorList[[newindex]],-newmaxlevel}]&)/@OperatorNames[level-newmaxlevel, newmaxlevel, newindex+1],
			OperatorNames[level, newmaxlevel, newindex+1]
		];
	,
		(* Boson *)	
		Return@Join[
			(Prepend[#,{OperatorList[[newindex]],-newmaxlevel}]&)/@OperatorNames[level-newmaxlevel, newmaxlevel, newindex],
			OperatorNames[level, newmaxlevel, newindex+1]
		];
	];
]
LLen[level_]:=Length[OperatorNames[level]]*HighestWeightDimension


StateVector[name_]:=StateVector[name,-Total[name[[1,All,2]]]]
StateVector[name_,level_]:=Module[{i,allnames},
allnames = StateNames[level];
For[i=1,i<=Length[allnames],i++,
StateVector[allnames[[i]],level]=StateVector[allnames[[i]]]=UnitVector[Length[allnames],i];
];Return@StateVector[name];]

	

GramMatrix[level_]:=Table[HCOperatorChain[name[[1]],0][[name[[2]],All]],{name,StateNames[level]}]//Simplify//Expand
GramNullspace[level_]:=NullSpace@GramMatrix@level
KacDeterminant[level_]:=Det[GramMatrix[level]]
GramRank[level_]:=MatrixRank[GramMatrix[level]]

Descendents[source_,sourcelevel_,targetlevel_]:=Descendents[source,sourcelevel,targetlevel]=
Which[targetlevel==sourcelevel,{source},True,MinimalRowReducedForm[Join@@Table[Which[
targetlevel-sourcelevel<-op[[2]],{},
True,
Descendents[source,sourcelevel,targetlevel+op[[2]]].Transpose[Operator[op[[1]],op[[2]],targetlevel+op[[2]]]]
],{op,MinimalGeneratorsCreation}]]]
DescendentsOfSpace[vecs_,startlevel_,endlevel_]:=RemoveZeroRows@RowReduce@(Join@@Table[Descendents[vec,startlevel,endlevel],{vec,vecs}])

FindEigenbasis[x_,level_]:=FindEigenbasis[x,IdentityMatrix[LLen[level]],level]


FindEigenbasis[xss_,ops_,level_]:=Module[{indices,b,binv,matrix,vals,vecs,xs,i,j,out,bases,newbases,valslist,scalarout},
If[!ListQ[xss],xs={xss};scalarout=True;,xs=xss;scalarout=False];
If[Length[xs]==1,
b =MinimalRowReducedForm[ops];
binv=PseudoInverse[b];
matrix = Transpose[binv].Operator[xs[[1]],0,level].Transpose[b];
{vals,vecs} = Eigensystem[matrix];
(*Print[MatrixForm@matrix,vals,vecs];*)

out = {};
valslist = DeleteDuplicates[vals];
Do[
indices = Select[Range@Length@vals,vals[[#]]==val&];
AppendTo[out,{Which[scalarout,val,True,{val}],vecs[[indices]].b}];
,{val,valslist}];

Return[out];
];
(* There is more than one operator. *)

(* Check for mutual commutivity. *)
For[i=1,i<=Length[xs],i++,
For[j=i+1,j<=Length[xs],j++,
If[!IsZeroMatrix[(Operator[xs[[i]]0,level].Operator[xs[[j]],0,level]-Operator[xs[[j]]0,level].Operator[xs[[i]],0,level]).Transpose[ops]],
Print["Operators ",xs[[i]]," and "xs[[j]]" do not commute at level ",level " ."];
Return[{}];
];
];
];
(* Everything should commute.  Now recurse. *)
bases = FindEigenbasis[xs[[2;;]],ops,level];
out = {};
Do[
newbases = FindEigenbasis[{xs[[1]]},base[[2]],level];
Do[AppendTo[out,{Join[k[[1]],base[[1]]],k[[2]]}],{k,newbases}]
,{base,bases}];
Return[out];
];



IntersectVectorSpaces[b1_,b2_,b3__]:=IntersectVectorSpaces[IntersectVectorSpaces[b1,b2],b3]
IntersectVectorSpaces[b1_,b2_]:=Module[{null,nullleft,out},
(* We assume that b1 and b2 are both bases.  We will use the minimality condition here. *)
null = NullSpace[Join[b1,b2]//Transpose];
nullleft = null[[All,1;;Length[b1]]];
out = MinimalRowReducedForm[-nullleft.b1];
Return@out;
]
UnionVectorSpaces[b___]:=MinimalRowReducedForm@Join[b]

IsZeroVector[vector_]:=And@@((#==0)&/@vector)
IsZeroMatrix[matrix_]:=And@@(IsZeroVector/@vector)
RemoveZeroRows[matrix_]:=Select[matrix,!IsZeroVector[#]&]
MinimalRowReducedForm[matrix_]:=Which[matrix=={},{},True,RemoveZeroRows@RowReduce@matrix]
EMatrixRank[matrix_]:=Which[matrix=={},0,True,MatrixRank@matrix];

FindSingularSubspace[level_]:=FindSingularSubspace[level]=FindSingularSubspace[NullSpace@GramMatrix[level],level]

FindSingularSubspace[basis_,level_]:=Module[{i,currentbasis,null},
currentbasis = basis;
If[Length[currentbasis]==0,Return@{};];
Do[
If[op[[2]]<=level,
null = NullSpace[Operator[op[[1]],op[[2]],level].Transpose[currentbasis]];
If[Length[null]==0,Break[];];
currentbasis=null.currentbasis;
];
,{op,MinimalGeneratorsAnnihilation}
];
If[Length[null]==0,Return@{};];
Return@currentbasis;
]

InVectorSpaceQ[space_,vector_]:=Module[{rowreduced},
rowreduced = RowReduce[space];
Return[MatrixRank[rowreduced]==MatrixRank[Append[rowreduced,vector]]];
]

ConstructLadderDiagram[cartanf_,l_]:=Module[{gcd,eigenbasis,cartan,singulars,m,i,k,pplevel,ppxs,level,xs,ppvecs,j,states,parentmatrix,vecs},
If[!ListQ[cartanf],cartan={cartanf};,cartan=cartanf];
gcd = GCD@@(Grading/@OperatorList);If[gcd==0,gcd=1;];
states={};
(* Find all the null singular vectors, and include all the weight 0 vectors. *)
For[i=0,i<=l,i=i+gcd,
If[i==0,singulars = IdentityMatrix@LLen@0;
,
Print[ "Searching for singular vectors at level ",i];
singulars = FindSingularSubspace[i];
];

If[Length[singulars]>0,
Print["Singular vectors have dimension ",singulars//Length,". "];
eigenbasis = FindEigenbasis[cartan,singulars,i];
If[Or@@(Length[#[[2]]>1]&/@eigenbasis),Print["Duplicate eigenvalues found for singular vector.  ConstructLadderDiagram will merge these eigenvalues."];];
(* Add the singular vectors to the list of states (level, x) *)
Do[
Print["Found singular vector at level ",i,", with eigenvalues ",eig[[1]],"."];
AppendTo[states,{i,eig[[1]],eig[[2]]}];
,{eig,eigenbasis}];
];
]; (*end dealing with the singulars*)

(* Now determine which states are descendents of each other. *)
parentmatrix = Table[Which[temp1<temp2,-1,True,0],{temp1,1,Length[states]},{temp2,1,Length[states]}]; (* parent, child *)

For[k=1,k<=Length[states],k++,
{level,xs,vecs}=states[[k]];

For[j=k-1,j>=1,j--,
	If[parentmatrix[[j,k]]==-1,(* We need to determine if j is a parent of k *)
		{pplevel,ppxs,ppvecs}=states[[j]];
		Print["Checking parentage of ",ppxs," to ",xs];
		If[IntersectQ[DescendentsOfSpace[ppvecs,pplevel,level],vecs],
			parentmatrix[[j,k]]=1;
			(* Make parents of parents parents *)
			For[m=1,m<=k,m++,
				If[parentmatrix[[m,j]]>0,parentmatrix[[m,k]]=2];
			];
			,
			parentmatrix[[j,k]]=0;
		];
	];
];
];
Return[{states,parentmatrix}];
];

LadderGraph[states_,parentmatrix_]:=Module[{edges,vertexlabels},vertexlabels=Table[i->Placed[ToString[states[[i,2]],InputForm],Top],{i,Length[states]}];
edges = {};
Do[If[parentmatrix[[i,j]]==1,AppendTo[edges,i->j]],{i,Length[states]},{j,Length[states]}];
Graph[edges,EdgeShapeFunction->GraphElementData["FilledArrow","ArrowSize"->0.05],VertexLabels->vertexlabels]
]

NullDescendents[level_]:=NullDescendents[level]=(Join@@Table[Which[level+ops[[2]]<0,{},Length[GramNullspace[level+ops[[2]]]]==0,{},True,GramNullspace[level+ops[[2]]].Transpose[Operator[ops[[1]],ops[[2]],level+ops[[2]]]]],{ops,MinimalGeneratorsCreation}]);


FindNonDescendents[level_]:=Module[{rank,matrix},
matrix=(Join@@Table[Which[level+ops[[2]]<0,{},Length[GramNullspace[level+ops[[2]]]]==0,{},True,GramNullspace[level+ops[[2]]].Transpose[Operator[ops[[1]],ops[[2]],level+ops[[2]]]]],{ops,MinimalGeneratorsCreation}]);
If[matrix=={},Return@MatrixRank[GramNullspace[level]];];
rank = MatrixRank@matrix;
Print[MatrixRank[GramNullspace[level]],", ",rank];
Return[MatrixRank[GramNullspace[level]]-rank];
]

FindAnnihilatingOperators[oplevel_,state_,level_]:=Module[{matrix},
If[oplevel<0,
matrix = Table[OperatorChain[list,level].state,{list,OperatorNames[-oplevel]}];
,
matrix = Table[OperatorChain[list,level].state,{list,PositiveOperatorNames[oplevel]}];
]
If[matrix=={},Return[{}];];
Return@NullSpace@Transpose@matrix;
]

InVectorSpaceQ[space_,vector_]:=Module[{rowreduced},
rowreduced = RowReduce[space];
Return[MatrixRank[rowreduced]==MatrixRank[Append[rowreduced,vector]]];
]

IntersectQ[space1_,space2_]:=Module[{s1,s2,rowreduced},
s1 = MinimalRowReducedForm[space1];
s2 = MinimalRowReducedForm[space2];
rowreduced = RowReduce[Join[s1,s2]];
Return[IsZeroVector[rowreduced[[-1]]]];
]

End[]
EndPackage[]
