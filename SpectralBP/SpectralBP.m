(* ::Package:: *)

(* ::Title:: *)
(*SpectralBP*)


(* ::Input::Initialization:: *)
BeginPackage["SpectralBP`"];


(* ::Chapter:: *)
(*Preamble*)


(* ::Subsection:: *)
(*Function usage*)


(* ::Subsubsection:: *)
(*Processing methods*)


(* ::Input::Initialization:: *)
DeconstructEquation::usage = "DeconstructEquation[eqn,Nmax]. Breaks down a set of ODEs into a list of matrices using a Bernstein polynomial basis of order Nmax"
MToLinear::usage = "MToLinear[M]. Converts the polynomial eigenproblem to a linear eigenproblem"
MListToLinear::usage = "MListToLinear[MList]. Converts a eigenproblem of several dependent variables to a linear eigenproblem."
GetDerMatrices::usage = "GetDerMatrices[Nmax,dermax,lbpower,ubpower,method,R,prec]. Computes the generic form of the function and its derivatives as a spectral matrix."
FuncToVector::usage = "FuncToVector[f,var,Nmax,ulb,uub,lbpower,ubpower,method,prec]. Computes the corresponding vector form of the function."


(* ::Subsubsection:: *)
(*Eigenvalue methods*)


(* ::Input::Initialization:: *)
GetModes::usage = "GetModes[eqn,n]. Generates the spectrum of eqn using a Bernstein polynomial basis of order n."
GetAccurateModes::usage = "GetAccurateModes[eqn,n1,n2]. Compares the spectrum of two different grid, and keeps the eigenvalues that agree up to a default value of 2 digits."
CompareModes::usage = "CompareModes[spectrum1,spectrum2]. Compares two spectra, and keeps the eigenvalues that agree up to a default value of 2 digits."
MagPrecision::usage = "MagPrecision[\[Omega]test,\[Omega]ref]. Returns the maximum number of digits the magnitude of \[Omega]test agrees with its closest counterpart in the list \[Omega]ref."
IndivPrecision::usage = "IndivPrecision[\[Omega]test,\[Omega]ref]. Returns the maximum number of digits that the real and imaginary parts of 
\[Omega]test agrees with its closest counterpart in the list \[Omega]ref."


(* ::Subsubsection:: *)
(*Eigenfunction methods*)


(* ::Input::Initialization:: *)
GetEigenvectors::usage = "GetEigenvectors[eqn,spectrum,Nmax]. Assumes that the list spectrum contains the eigenvalues of eqn 
using a spectral basis of order N, and returns the corresponding eigenvectors."
GetEigenvector::usage = "GetEigenvector[MList,\[Omega],normalization,prec]. Calculates the eigenvector of the eigenproblem (MList, \[Omega]), to machine precision prec."
GetEigenfunctions::usage = "GetEigenfunctions[eqn,spectrum,Nmax,var]. Calculates the eigenfunctions as a sum of Bernstein polynomials, from solving the eigenvectors of a list of eigenvalues."
CompareEigenfunctions::usage="CompareEigenfunctions[eqn,spectra,Nmaxes]. Compares the L2-norm of the difference between eigenfunctions as a test for convergence."


(* ::Input::Initialization:: *)
PrintAll::usage = "PrintAll[eqn,{eigenhigh,eigenlow},N]. Prints a frequency plot in the complex plane, a table of frequencies common to both eigenhigh and eigenlow,
 and plot of corresponding eigenfunctions."
PrintFrequencies::usage = "PrintFrequencies[spectrum]. Prints a frequency plot in the complex plane."
PrintTable::usage = "PrintTable[{eigenhigh,eigenlow}]. Prints a table of frequencies common to both eigenhigh and eigenlow, keeping the digits that they agree in.
Separates damped and overdamped modes."
PrintEigenfunctions::usage = "PrintEigenfunctions[eqn,spectrum,N]. Assumes that the list spectrum contains the eigenvalues of eqn 
using a spectral basis of order N. Returns a plot of the real and imaginary parts of the corresponding eigenfunctions."
CleanDigits::usage = "CleanDigits[\[Omega]test,\[Omega]ref]. Returns the eigenvalue \[Omega]test to the precision that it shares the most digits of some eigenvalue in \[Omega]ref."


(* ::Subsubsection:: *)
(*Error handling*)


(* ::Input::Initialization:: *)
ThrowAbort::usage = "ThrowAbort[message,args]."
CatchAbort::usage = "CatchAbort[exp]."


(* ::Subsection:: *)
(*Global variables*)


(* ::Input::Initialization:: *)
Begin["`Private`"];


(* ::Chapter:: *)
(*Calculations*)


(* ::Section:: *)
(*Processing methods*)


(* ::Text:: *)
(*In this section, we break down the equation so that it may converted into a matrix eigenproblem. Current things to do: (1) make it accept a coupled set of ODE's (2) make it accept non-polynomial type coefficient functions (3) make it accept numerical backgrounds.*)


(* ::Input::Initialization:: *)
(*02/03/2019*)
DeconstructEquation[equation_List,Nmax_,opts : OptionsPattern[{BoundaryOptions,NumericalOptions}]]:=Block[
{eqn,funclist,f,fVector,var,\[Omega],\[Omega]max,dermax,funcmax,\[Omega]coefficients,dercoefficients,MList,DerMatrices,ulb,uub,testlist,prec = If[OptionValue[FunctionPrecision] == 0, Max[2*Nmax,$MinPrecision], OptionValue[FunctionPrecision]],
method=OptionValue[CollocationGrid],lbpower=OptionValue[LBPower],ubpower=OptionValue[UBPower]},

ulb = SetPrecision[OptionValue[LowerBound],prec];uub=SetPrecision[OptionValue[UpperBound],prec];
eqn = Expand@equation;

(*find independent variable*)
testlist = Union[Cases[eqn,\[Phi]_Symbol[u_Symbol]-> u, \[Infinity]],Cases[eqn,Derivative[n_][\[Phi]_Symbol][u_Symbol]-> u,\[Infinity]]];
Which[Length[testlist]<1,ThrowAbort[DeconstructEquation::novariables];,
Length[testlist]>1,ThrowAbort[DeconstructEquation::manyvariables,testlist];,
True,var = Last@testlist;];

(*find eigenvariable*)
testlist = Union[Cases[eqn, _Symbol, \[Infinity]]] //DeleteCases[#,var]&;
Which[Length[testlist]<1,ThrowAbort[DeconstructEquation::noeigenvalues];,
Length[testlist]>1,ThrowAbort[DeconstructEquation::manyeigenvalues,testlist];,
True, \[Omega] = Last@testlist;];

(*find dependent functions*)
testlist =  Union[Cases[eqn,\[Phi]_Symbol[var]/;Context[\[Phi]]!= "System`"->  \[Phi],\[Infinity]],Cases[eqn,Derivative[n_][\[Phi]_Symbol][var]-> \[Phi],\[Infinity]]];

Which[Length[testlist]<Length[eqn],ThrowAbort[DeconstructEquation::underdetermined,testlist,Length[eqn]];,
Length[testlist]>Length[eqn],ThrowAbort[DeconstructEquation::overdetermined,testlist,Length[eqn]];,
True, funclist = testlist;];

funcmax = Length[funclist];
\[Omega]max = Max@Cases[eqn,\[Omega]^(n_:1)->n,\[Infinity]];
dermax = (Cases[eqn,Derivative[n_][#][var]->n,\[Infinity]]&/@funclist//Max)/.-Infinity->0;

\[Omega]coefficients = CoefficientList[eqn,\[Omega],\[Omega]max+1];

dercoefficients = CatchAbort[
Table[FuncToVector[Coefficient[\[Omega]coefficients[[eqnindex,\[Omega]index]],Derivative[derindex][funclist[[funcindex]]][var]],var,Nmax,ulb,uub,lbpower,ubpower,method,prec],{eqnindex,funcmax},{funcindex,funcmax},{\[Omega]index,\[Omega]max+1},{derindex,0,dermax}]
];

DerMatrices =CatchAbort[
GetDerMatrices[Nmax, dermax,lbpower,ubpower,method,(uub-ulb),prec]
];

MList= Table[Sum[
Switch[dercoefficients[[eqnindex,funcindex,\[Omega]index,derindex]],
0,0,
_,Table[dercoefficients[[eqnindex,funcindex,\[Omega]index,derindex,j+1]]*DerMatrices[[derindex,j+1,k+1]],{j,0,Nmax},{k,0,Nmax}]
],
{derindex,dermax+1}],
{eqnindex,funcmax},{funcindex,funcmax},{\[Omega]index,\[Omega]max+1}];
MList
]
DeconstructEquation[equation_,Nmax_,opts : OptionsPattern[{BoundaryOptions,NumericalOptions}]]:=DeconstructEquation[{equation},Nmax,opts]


(* ::Input::Initialization:: *)
DeconstructEquation::novariables= "No independent variables detected."
DeconstructEquation::manyvariables = "Too many independent variables detected. `1` detected, expecting 1."
DeconstructEquation::noeigenvalues = "No free parameter detected."
DeconstructEquation::manyeigenvalues = "Too many free parameters detected. `1` detected, expecting 1."
DeconstructEquation::underdetermined = "Problem is underdetermined. `1` dependent variables detected, expecting `2`."
DeconstructEquation::overdetermined = "Problem is overdetermined. `1` dependent variables detected, expecting `2`."


(* ::Text:: *)
(*This function converts the polynomial eigenproblem to a linear eigenproblem.*)


(* ::Input::Initialization:: *)
(*02/03/2019*)
MToLinear[M_,IsLinear_:False]:= Block[{A,B,Nmax,\[Omega]max,Mfinalinv,Mfirstinv},
\[Omega]max = Length[M]-1;Nmax = Length[M[[1]]];$MinPrecision=Precision[M[[1]]];
Which[Det[M[[\[Omega]max+1]]]!= 0&&IsLinear,
(Mfinalinv = Inverse[-M[[\[Omega]max+1]]];
A = Table[If[i==\[Omega]max,Mfinalinv.M[[j]],KroneckerDelta[i+1,j]*IdentityMatrix[Nmax]],{i,\[Omega]max},{j,\[Omega]max}]//ArrayFlatten;
B = 1;),
Det[M[[1]]]!= 0&&IsLinear,
(Mfirstinv = Inverse[M[[1]]];
A = 1;
B = Table[If[i==1,-Mfirstinv.M[[j+1]],KroneckerDelta[i,j+1]*IdentityMatrix[Nmax]],{i,\[Omega]max},{j,\[Omega]max}]//ArrayFlatten;),
True,
(A = Table[If[i==1,M[[j]],KroneckerDelta[i,j]*IdentityMatrix[Nmax]],{i,\[Omega]max},{j,\[Omega]max}] //ArrayFlatten;
B = Table[If[i==1 && j == \[Omega]max , -M[[\[Omega]max+1]],KroneckerDelta[i,j+1]*IdentityMatrix[Nmax]],{i,\[Omega]max},{j,\[Omega]max}]//ArrayFlatten;)];

$MinPrecision=0;
{A,B}
]


(* ::Input::Initialization:: *)
MListToLinear[MList_]:=Block[{MlinearList,A,B,funcmax,IsLinear},
funcmax = Length[MList];
If[funcmax == 1,IsLinear=True;,IsLinear=False;];
MlinearList = Table[MToLinear[MList[[eqnindex,funcindex]],IsLinear],{eqnindex,funcmax},{funcindex,funcmax}];
A = Table[MlinearList[[eqnindex,funcindex,1]],{eqnindex,funcmax},{funcindex,funcmax}]//ArrayFlatten;
B = Table[MlinearList[[eqnindex,funcindex,2]],{eqnindex,funcmax},{funcindex,funcmax}]//ArrayFlatten;
{A,B}
]


(* ::Text:: *)
(*Calculates how \[Phi](u), d\[Phi]/du, (d^2 \[Phi])/du^2,etc. is converted into a spectral matrix.*)


(* ::Input::Initialization:: *)
(*12/03/2019*)
GetDerMatrices[Nmax_, dermax_,lbpower_,ubpower_,method_,R_,prec_] := Block[{CMatrix,DerMatrices,Ntotal,klowerbound,kupperbound},
Ntotal= Nmax + lbpower + ubpower; klowerbound = Max[-dermax,-lbpower]; kupperbound = Nmax + Min[ubpower,dermax];
CMatrix = Switch[method,
"Chebyschev", Table[Switch[k,-lbpower,N[(1+Cos[(j+lbpower)/Ntotal \[Pi]])^Ntotal,prec],Nmax + ubpower,N[(1-Cos[(j+lbpower)/Ntotal \[Pi]])^Ntotal,prec],_,N[(1-Cos[(j+lbpower)/Ntotal \[Pi]])^(k+lbpower) (1+Cos[(j+lbpower)/Ntotal \[Pi]])^(Nmax+ubpower - k),prec]],{j,0,Nmax},{k,klowerbound,kupperbound}],
"Equal", Table[Switch[k,-lbpower,N[(1-(j+lbpower)/Ntotal)^Ntotal,prec],Nmax + ubpower, N[((j+lbpower)/Ntotal)^Ntotal,prec],_,N[((j+lbpower)/Ntotal)^(k+lbpower) (1-(j+lbpower)/Ntotal)^(Nmax+ubpower-k),prec]],{j,0,Nmax},{k,klowerbound,kupperbound}],
_,ThrowAbort[GetDerMatrices::method]
];
DerMatrices = Table[Sum[If[klowerbound <= k+l+m-n<= kupperbound,(1/R)^n (-1)^l Binomial[n,l]Binomial[n,m]Binomial[Ntotal,k+l+lbpower]Pochhammer[k+l+lbpower+1-n,n]CMatrix[[j+1,k+l+m-n+1-klowerbound]],0],{l,0,n},{m,0,n}],{n,0,dermax},{j,0,Nmax},{k,0,Nmax}];
If[MemberQ[DerMatrices,Infinity]||MemberQ[DerMatrices,Indeterminate],ThrowAbort[GetDerMatrices::infinity]];
DerMatrices
]


(* ::Input::Initialization:: *)
GetDerMatrices::method = "Unknown method."
GetDerMatrices::infinity = "Infinity or indeterminate quantity encountered in generating function spectral matrices."


(* ::Text:: *)
(*Converts the coefficient function f(u) into the corresponding vector for the spectral method.*)


(* ::Input::Initialization:: *)
(*13/03/2019*)
FuncToVector[f_,var_,Nmax_,ulb_,uub_,lbpower_,ubpower_,method_,prec_]:=Block[{funcVector,Ntotal = Nmax+lbpower+ubpower},
If[f[var]==0,Return[0]];
funcVector = Switch[method,
"Chebyschev",Table[N[f /. var -> N[(uub-ulb)/2 (1-Cos[(j + lbpower)/Ntotal \[Pi]])+ulb,prec],prec],{j,0,Nmax}],
"Equal",Table[N[f /. var -> N[((uub-ulb)(j+lbpower))/Ntotal+ulb,prec],prec],{j,0,Nmax}],
_,ThrowAbort[FuncToVector::method]];
If[MemberQ[funcVector,Infinity]||MemberQ[funcVector,Indeterminate],ThrowAbort[FuncToVector::infinity]];
funcVector
]


(* ::Input::Initialization:: *)
FuncToVector::method = "Unknown method."
FuncToVector::infinity= "Infinity or indeterminate quantity encountered in evaluating coefficient function vectors."


(* ::Section:: *)
(*Calculation of modes*)


(* ::Text:: *)
(*Calculates the spectrum of eqn from a Bernstein basis of order n, with machine precision prec. *)


(* ::Input::Initialization:: *)
(*04/23/2019*)
GetModes[eqn_,{n_,prec_},opts:OptionsPattern[{BoundaryOptions,NumericalOptions}]]:= Block[{MList,A,B,EVs},
MList = DeconstructEquation[eqn,n,opts];
{A,B} = MListToLinear[SetPrecision[MList,prec]];
A=SetPrecision[A,prec];B=SetPrecision[B,prec];
Which[B == IdentityMatrix[First@Dimensions[B]],A=SetPrecision[A,prec];EVs = Eigenvalues[A];,
A == IdentityMatrix[First@Dimensions[A]], If[MemberQ[EVs,0],Message[GetModes::infinity];EVs=DeleteCases[EVs,0];];EVs = 1/Eigenvalues[B];, 
True,EVs = Eigenvalues[{A,B}];
];
If[MemberQ[EVs,Infinity]||MemberQ[EVs,Indeterminate],Message[GetModes::infinity]; EVs=DeleteCases[DeleteCases[EVs,Infinity],Indeterminate];];
Sort[EVs]
]
GetModes[eqn_,n_Integer,opts:OptionsPattern[{BoundaryOptions,NumericalOptions}]]:= GetModes[eqn,{n,n/2},opts]


(* ::Input::Initialization:: *)
GetModes::infinity = "Infinite or indeterminate eigenvalues encountered. Calculation continues, but we have removed these eigenvalues.";


(* ::Text:: *)
(*Calculates the spectrum of eqn using two Bernstein basis of order n1 and n2, calculated with machine precision prec1 and prec2 respectively, and keeps modes shared by both spectra that agree to a Cutoff number of digits. One may set whether the agreement is in the magnitude of the modes, or in the digits shared by the real and imaginary parts separately, toggled by TestMagnitude.*)


(* ::Input::Initialization:: *)
(*04/23/2019*)
GetAccurateModes[eqn_,grids_List,opts : OptionsPattern[{ComparisonOptions,BoundaryOptions,NumericalOptions}]] :=Block[{sortedgrid,spectra},
sortedgrid = SortBy[grids,First];
spectra = CatchAbort[Table[GetModes[eqn,sortedgrid[[k]],FilterRules[{opts},{LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid}]],{k,Length[sortedgrid]}]];
CompareModes[spectra,FilterRules[{opts},{Cutoff,TestMagnitude}]]
]

GetAccurateModes[eqn_,{n1_,prec1_},{n2_,prec2_},opts : OptionsPattern[{ComparisonOptions,BoundaryOptions,NumericalOptions}]] := GetAccurateModes[eqn,{{n1,prec1},{n2,prec2}},opts]
GetAccurateModes[eqn_,n1_Integer,n2_Integer, opts : OptionsPattern[{ComparisonOptions,BoundaryOptions,NumericalOptions}]] := GetAccurateModes[eqn,{{n1,n1/2},{n2,n2/2}},opts]


(* ::Text:: *)
(*Compares two spectra, and keeps modes shared by both spectra that agree to a Cutoff number of digits. One may set whether the agreement is in the magnitude of the modes, or in the digits ahred by the real and imaginary parts separately, toggled by TestMagnitude.*)


(* ::Input::Initialization:: *)
(*05/03/2019*)
CompareModes[spectra_List,opts : OptionsPattern[ComparisonOptions]] :=Block[{sortedspectra,filteredspectra,
cutoff=OptionValue[Cutoff],testmagnitude = OptionValue[TestMagnitude]},
sortedspectra= Reverse[SortBy[spectra,Length[#]&]];
filteredspectra = Table[If[testmagnitude,Select[sortedspectra[[k]],MagPrecision[#,sortedspectra[[1;;k-1]]~Join~sortedspectra[[k+1;;Length[sortedspectra]]]]>= cutoff &],Select[sortedspectra[[k]],IndivPrecision[#,sortedspectra[[1;;k-1]]~Join~sortedspectra[[k+1;;Length[sortedspectra]]]]>= cutoff &]],{k,Length[sortedspectra]}]
]
CompareModes[spectrum1_List,spectrum2_List,opts : OptionsPattern[ComparisonOptions]]  := (CompareModes[{spectrum1,spectrum2},opts])


(* ::Text:: *)
(*Compares a frequency \[Omega]test to a list of frequencies \[Omega]ref by taking the closest frequency from the list, returning the number of agreed digits.*)


(* ::Input::Initialization:: *)
(*02/03/2019*)
MagPrecision[\[Omega]test_,\[Omega]ref_] := Block[{avemag,agreedDigits,\[Omega]samelist,\[Omega]same},
\[Omega]samelist = Table[First@MinimalBy[\[Omega]ref[[k]],Abs[#-\[Omega]test]&],{k,Length[\[Omega]ref]}];
\[Omega]same = First@MaximalBy[\[Omega]samelist,Abs[#-\[Omega]test]&];
avemag = Floor[Log10[Abs[(\[Omega]test+\[Omega]same)/2]]/. -Infinity->0];
agreedDigits = (-Ceiling[Log10[Abs[\[Omega]test-\[Omega]same]]/. Indeterminate -> -Min[Precision/@{\[Omega]test,\[Omega]same}]])+avemag;
agreedDigits
]


(* ::Text:: *)
(*Compares a frequency \[Omega]test to a list of frequencies \[Omega]ref by taking the closest frequency from the list, returning the smaller of the agreed digits of the real and imaginary parts.*)


(* ::Input::Initialization:: *)
(*05/03/2019*)
IndivPrecision[\[Omega]test_,\[Omega]ref_] := Block[{reave\[Omega],imave\[Omega],re\[Omega]test,im\[Omega]test,re\[Omega]same,im\[Omega]same,reagreedDigits,imagreedDigits,\[Omega]samelist,\[Omega]same},
\[Omega]samelist = Table[First@MinimalBy[\[Omega]ref[[k]],Abs[#-\[Omega]test]&],{k,Length[\[Omega]ref]}];
\[Omega]same = First@MaximalBy[\[Omega]samelist,Abs[#-\[Omega]test]&];
re\[Omega]same = Re[\[Omega]same];
im\[Omega]same = Im[\[Omega]same];
re\[Omega]test = Re[\[Omega]test];
im\[Omega]test = Im[\[Omega]test];
reave\[Omega] = Floor[Log10[Abs[(re\[Omega]test+re\[Omega]same)/2]]/. -Infinity-> 0];

imave\[Omega] = Floor[Log10[Abs[(im\[Omega]test+im\[Omega]same )/2]]/. -Infinity -> 0];
reagreedDigits = (-Floor[Ceiling[Log10[Abs[re\[Omega]test-re\[Omega]same]] /. Indeterminate -> -Min[Precision/@{\[Omega]test,\[Omega]same}],0.5]] )+ reave\[Omega];
imagreedDigits = (-Floor[Ceiling[Log10[Abs[im\[Omega]test-im\[Omega]same]] /. Indeterminate -> -Min[Precision/@{\[Omega]test,\[Omega]same}],0.5]]) + imave\[Omega];
Min[reagreedDigits,imagreedDigits]
]


(* ::Section:: *)
(*Calculation of eigenfunctions*)


(* ::Input::Initialization:: *)
GetL2Norm[eigenvector_,R_,normalization_,lbpower_,ubpower_,prec_]:=Block[{A,Nmax,N1,N2,Atilde,n,m},
Switch[normalization,
_List,If[normalization[[1]]!= "L2Norm",ThrowAbort[GetL2Norm::method]];
Atilde=normalization[[2,1]];n=normalization[[2,2]];m=normalization[[2,3]];,
_,Atilde=1;n=0;m=0;
];

If[n+2*lbpower<0||m+2*ubpower<0,ThrowAbort[GetL2Norm::notnormalizable]];
Nmax = Length[eigenvector];N1 = Nmax + lbpower + ubpower - 1; N2 = 2N1+n+m;
Abs@N[Sqrt[Sum[eigenvector[[j]]Conjugate[eigenvector[[k]]] (Binomial[N1,j-1+lbpower]Binomial[N1,k-1+lbpower])/Binomial[N2, j+ k + 2*lbpower + n-2] (Atilde *R^(n+m+1))/(N2 + 1) ,{j,Nmax},{k,Nmax}]],prec]
]


(* ::Input::Initialization:: *)
GetL2Norm::method="Unknown L2 normalization method";
GetL2Norm::notnormalizable="Weight function makes function not normalizable";


(* ::Input::Initialization:: *)
(*13/03/2019*)
GetEigenvectors[eqn_,spectrum_List,Nmax_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:=Block[{MList,eigenvectors,prec},
prec = If[OptionValue[FunctionPrecision] == 0, Max[2*Nmax,$MinPrecision], OptionValue[FunctionPrecision]];
MList = DeconstructEquation[eqn,Nmax,FilterRules[{opts},{LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid}]];
eigenvectors = CatchAbort[
Table[GetEigenvector[MList,spectrum[[j]],OptionValue[Normalization],prec,FilterRules[{opts},{LowerBound,UpperBound,LBPower,UBPower}]],{j,Length[spectrum]}]
];
eigenvectors
]
GetEigenvectors[eqn_,spectrum_,Nmax_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions,NumericalOptions}]] := 
GetEigenvectors[eqn,{spectrum},Nmax,opts]


(* ::Input::Initialization:: *)
(*13/03/2019*)
GetEigenvector[MList_,\[Omega]_,normalization_,prec_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions}]]:=Block[{M,Minv,eigenvector,eigenvectortest,counter=2,subvector,residue,matrix,Nmax,funcmax,ulb = SetPrecision[OptionValue[LowerBound],2*prec],uub=SetPrecision[OptionValue[UpperBound],2*prec],lbpower=OptionValue[LBPower],ubpower=OptionValue[UBPower],R,jacobian,N1,N2,A},
funcmax = Length[MList];
M = Table[Sum[N[\[Omega]^(j-1) MList[[eqnindex,funcindex,j]],prec],{j,Length[MList[[eqnindex,funcindex]]]}],{funcindex,funcmax},{eqnindex,funcmax}]//ArrayFlatten;

Nmax = Dimensions[M]//First;
Minv= Quiet[Inverse[N[M,prec]]];
eigenvector = Minv.UnitVector[Nmax,1]; eigenvector = SetPrecision[eigenvector/Max@Abs@eigenvector,prec];
eigenvectortest = Minv.eigenvector; eigenvectortest = SetPrecision[eigenvectortest/Max@Abs@eigenvectortest,prec];

While[Log10[Max[Abs[eigenvector-eigenvectortest]]]>-OptionValue[FunctionCutoff]&&counter<OptionValue[MaxIteration],(eigenvector=eigenvectortest; eigenvectortest = Minv.eigenvector; eigenvectortest = SetPrecision[eigenvectortest/Max@Abs@eigenvectortest,prec];);counter++];
eigenvector=eigenvectortest;

Nmax = Length[eigenvector]/funcmax;
eigenvector = Table[(subvector = eigenvector[[(k-1)Nmax+1;;k*Nmax]];
Switch[normalization,
"UB", If[subvector[[-1]]==0,Message[GetEigenvector::boundaryconditionsub];SetPrecision[subvector,prec],SetPrecision[subvector/subvector[[-1]],prec]],
"LB",If[subvector[[1]]==0,Message[GetEigenvector::boundaryconditionslb];SetPrecision[subvector,prec],SetPrecision[subvector/subvector[[1]],prec]],
"L2Norm",(
subvector = subvector/subvector[[1]];
A = GetL2Norm[subvector,uub-ulb,normalization,lbpower,ubpower,prec];
If[A==0,Message[GetEigenvector::boundaryconditionsL2];SetPrecision[subvector,prec],
SetPrecision[subvector/A,prec]]
),
_List,(
subvector = subvector/subvector[[1]];
A = GetL2Norm[subvector,uub-ulb,normalization,lbpower,ubpower,prec];
If[A==0,Message[GetEigenvector::boundaryconditionsL2];SetPrecision[subvector,prec],
SetPrecision[subvector/A,prec]]
),
_,ThrowAbort[GetEigenvector::method]
]),{k,funcmax}];
eigenvector
]


(* ::Input::Initialization:: *)
GetEigenvector::method = "Unknown normalization method. Either 'LB' or 'UB' or a normalization pair";
GetEigenvector::boundaryconditionsub= "Unexpeted vanishing leading polynomial at the upper bound";
GetEigenvector::boundaryconditionslb= "Unexpected vanishing leading polynomial at the lower bound";
GetEigenvector::boundaryconditionsL2norm= "Unexpected vanishing of the L2-norm";


(* ::Input::Initialization:: *)
GetEigenfunctions[eqn_,spectrum_,Nmax_,var_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:=Block[{eigenfunctions,eigenvectors,BPs,R,eigenmax,funcmax,finallbpower,finalubpower,prec,Cmatrix,
ulb = OptionValue[LowerBound],uub=OptionValue[UpperBound],lbpower = OptionValue[LBPower],ubpower=OptionValue[UBPower]},

prec = If[OptionValue[FunctionPrecision] == 0, Max[2*Nmax,$MinPrecision], OptionValue[FunctionPrecision]];
ulb = SetPrecision[OptionValue[LowerBound],prec];uub = SetPrecision[OptionValue[UpperBound],prec];R = uub-ulb;
{finallbpower,finalubpower}=If[OptionValue[FinalAsymptotics]==="Default",{lbpower,ubpower},OptionValue[FinalAsymptotics]];
Print[OptionValue[FinalAsymptotics]];
Print[finallbpower];

BPs = Table[Binomial[Nmax+lbpower+ubpower,l+lbpower]/R^(Nmax+lbpower+ubpower) Switch[l,
-finallbpower,(uub-var)^(Nmax+finallbpower+finalubpower),
Nmax+finalubpower,(var-ulb)^(Nmax+finallbpower+finalubpower),
_,(var-ulb)^(l+finallbpower) (uub-var)^(Nmax+finalubpower-l)],{l,0,Nmax}];

eigenvectors = GetEigenvectors[eqn,spectrum,Nmax,FilterRules[{opts},{Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]];
eigenmax = Dimensions[eigenvectors][[1]];
funcmax = Dimensions[eigenvectors][[2]];

eigenfunctions = Table[Sum[eigenvectors[[eigenindex,funcindex,k]]*BPs[[k]],{k,Nmax+1}],{eigenindex,eigenmax},{funcindex,funcmax}];
eigenfunctions
]


(* ::Input::Initialization:: *)
CompareEigenfunctions[eqn_List,spectra_List,Nmaxes_List,opts : OptionsPattern[{EigenCompOptions,ComparisonOptions,EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:=Block[{Nmax=Max[Nmaxes],R,filteredspectra={},prec,eigenvectors,testeigenvector,subvector,eigenvalues,var,meanerrors,Neigenfunc=Length[eqn],sortedspectra,outputspectra,\[Omega]test,lbpower=OptionValue[LBPower],ubpower=OptionValue[UBPower]},
prec = If[OptionValue[FunctionPrecision] == 0, Max[2*Nmax,$MinPrecision], OptionValue[FunctionPrecision]];
R = SetPrecision[OptionValue[UpperBound] - OptionValue[LowerBound],prec];
sortedspectra = FilterEigenvalues[spectra[[1]],True];
For[i=1,i<=First@Dimensions[sortedspectra],i++,(\[Omega]test=SetPrecision[sortedspectra[[i,3]],prec];
eigenvectors= Table[testeigenvector=GetEigenvectors[eqn,{\[Omega]test},Nmaxes[[k]],FilterRules[{opts},{Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]];

(*Expand to common highest basis order*)
subvector=testeigenvector[[1,funcindex]];
Table[
Sum[If[0<= m-l<= Nmax-Nmaxes[[k]],
N[(Binomial[Nmaxes[[k]]+lbpower+ubpower,l+lbpower]Binomial[Nmax-Nmaxes[[k]],m-l])/Binomial[Nmax+lbpower+ubpower,m+lbpower],prec] subvector[[l+1]],
0] ,
{l,Max[m-Nmax+Nmaxes[[k]],0],Min[m,Nmaxes[[k]]]}],
{m,0,Nmax}]
,{k,Length[Nmaxes]},{funcindex,Neigenfunc}];

(*calculate L2Norm of the difference of eigenvectors*)
meanerrors = Table[subvector=eigenvectors[[1,funcindex]]-eigenvectors[[k,funcindex]];
GetL2Norm[subvector,R,OptionValue[Normalization],lbpower,ubpower,prec]
,{k,2,Length[Nmaxes]},{funcindex,Neigenfunc}];

(*filter out eigenvalues whose eigenvectors did not converge*)
If[Log10[Max[meanerrors]]<-OptionValue[L2Cutoff],filteredspectra=filteredspectra~Join~Which[sortedspectra[[i,1]]=="\[PlusMinus] ",{First@MinimalBy[spectra[[1]],Abs[#-(\[Omega]test)]&],First@MinimalBy[spectra[[1]],Abs[#-(-Conjugate[\[Omega]test])]&]},
sortedspectra[[i,2]]=="\[PlusMinus] ",{First@MinimalBy[spectra[[1]],Abs[#-(\[Omega]test)]&],First@MinimalBy[spectra[[1]],Abs[#-(Conjugate[\[Omega]test])]&]},
True,{First@MinimalBy[spectra[[1]],Abs[#-(\[Omega]test)]&]}];];)];
If[Length[filteredspectra]==0,{},CompareModes[{filteredspectra}~Join~spectra[[2;;-1]],FilterRules[{opts},{Cutoff,TestMagnitude}]]]
]
CompareEigenfunctions[eqn_,spectra_List,Nmaxes_List,opts : OptionsPattern[{EigenCompOptions,ComparisonOptions,EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:=CompareEigenfunctions[{eqn},spectra,Nmaxes,opts]


(* ::Section:: *)
(*Outputting*)


(* ::Input::Initialization:: *)
(*05/03/2019*)
PrintAll[eqn_,{eigenhigh_,eigenlow_},Nmax_,opts:OptionsPattern[{PrintOptions,EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:=Block[{nspectrum=OptionValue[NSpectrum],neigenfunc=OptionValue[NEigenFunc],cutoff=OptionValue[ZeroCutoff],precList,frequencyplot,eigenUndamped,eigenDamped,eigenOverdamped,nUndamped,nDamped,nOverdamped,eigenfuncundamped,eigenfuncdamped,eigenfuncoverdamped,table,results},

frequencyplot = PrintFrequencies[Reverse[eigenhigh],Sequence[FilterRules[{opts},{FreqName,NSpectrum}]]];
table = PrintTable[{eigenhigh,eigenlow},Sequence[FilterRules[{opts},{ZeroCutoff,FreqName}]]];

results = {{"Eigenvalues"},{frequencyplot},{"Table of Eigenvalues"},{table}};

eigenUndamped = Sort[Select[eigenhigh,(Log10[Abs[Im[#]]]/.Indeterminate -> -Infinity)<-cutoff&],Re[#1]<Re[#2]&];
nUndamped = Length[eigenUndamped];
If[nUndamped>0,
results = results~Join~{{"Undamped Eigenfunctions"}}~Join~{{PrintEigenfunctions[eqn,eigenUndamped[[1;;Min[neigenfunc,nUndamped]]],Nmax,FilterRules[{opts},{FinalAsymptotics,Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]]}}
];

eigenOverdamped = Sort[Select[eigenhigh,(Log10[Abs[Re[#]]]/.Indeterminate -> -Infinity)<-cutoff&&(Log10[Abs[Im[#]]]/.Indeterminate -> -Infinity)>-cutoff &],Im[#1]<Im[#2]&];
nOverdamped = Length[eigenOverdamped];
If[nOverdamped > 0,results=results~Join~{{"Overdamped Eigenfunctions"}}~Join~{{PrintEigenfunctions[eqn,eigenOverdamped[[1;;Min[neigenfunc,nOverdamped]]],Nmax,FilterRules[{opts},{FinalAsymptotics,Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]]}}
];

eigenDamped =FilterEigenvalues[Select[eigenhigh,(Log10[Abs[Re[#]]]/.Indeterminate -> -Infinity)>-cutoff&&(Log10[Abs[Im[#]]]/.Indeterminate -> -Infinity)>-cutoff &],False];
nDamped = Length[eigenDamped];
If[nDamped>0,results = results~Join~{{"Damped Eigenfunctions"}}~Join~{{PrintEigenfunctions[eqn,eigenDamped[[1;;Min[neigenfunc,nDamped]]],Nmax,FilterRules[{opts},{FinalAsymptotics,Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]]}}
];
results = Grid[results];
results
]


(* ::Input::Initialization:: *)
(*14/03/2019*)
PrintFrequencies[spectrum_, opts: OptionsPattern[PrintOptions]] := Block[{nspectrum = OptionValue[NSpectrum],plot},
plot = ListPlot[(Tooltip[{Re[#1],Im[#1]}]&)/@(Sort[spectrum][[1;;nspectrum]]),PlotRange-> All,Axes->True,ImageSize-> 800,Frame->True, BaseStyle->{FontSize->24},LabelStyle->{FontSize->24},FrameLabel->{"Re "<>OptionValue[FreqName],"Im "<>OptionValue[FreqName]}];
plot
]


(* ::Input::Initialization:: *)
PrintTable[{eigenhigh_,eigenlow_},opts:OptionsPattern[PrintOptions]]:=Block[{table,tableUndamped,tableDamped,tableOverdamped,n=Length[eigenhigh],nReal,nImaginary,nComplex,eigensame,eigensameReal,eigensameImaginary,eigensameComplex,cutoff = OptionValue[ZeroCutoff]},
eigensame = Table[CleanDigits[eigenhigh[[i]],eigenlow],{i,n}];

eigensameReal = Sort[Select[eigensame,(Log10[Abs[Im[#]]]/.Indeterminate -> -Infinity)<-cutoff&],Re[#1]<Re[#2]&];
eigensameImaginary = Sort[Select[eigensame,(Log10[Abs[Re[#]]] /.Indeterminate -> -Infinity)< -cutoff &&( Log10[Abs[Im[#]]]/.Indeterminate -> -Infinity)>-cutoff&],Im[#1]>Im[#2]&];
eigensameComplex = FilterEigenvalues[Select[eigensame,(Log10[Abs[Im[#]]] /.Indeterminate -> -Infinity)> -cutoff && (Log10[Abs[Re[#]]]/.Indeterminate -> -Infinity)>-cutoff&],True];

nReal = Length[eigensameReal];
nImaginary = Length[eigensameImaginary];
nComplex = Length[eigensameComplex];

table = {};
If[nReal >0, table = table~Join~{{"Real Eigenvalues",SpanFromLeft,SpanFromLeft},{"n", "Re " <>OptionValue[FreqName],"Im "<>OptionValue[FreqName]}}~Join~Table[{i,ToString[Re[eigensameReal[[i]]],TraditionalForm],ScientificForm[Im[eigensameReal[[i]]],1]},{i,nReal}]];
If[nImaginary > 0, table = table~Join~{{"Imaginary Eigenvalues",SpanFromLeft,SpanFromLeft},{"n", "Re " <>OptionValue[FreqName],"Im "<>OptionValue[FreqName]}}~Join~Table[{i,ScientificForm[Re[eigensameImaginary[[i]]],1],Im[eigensameImaginary[[i]]]},{i,nImaginary}]];
If[nComplex >0, table = table~Join~{{"Complex Eigenvalues",SpanFromLeft,SpanFromLeft},{"n", "Re " <>OptionValue[FreqName],"Im "<>OptionValue[FreqName]}}~Join~Table[{i,eigensameComplex[[i,1]]<>ToString[Re[eigensameComplex[[i,3]]],TraditionalForm],eigensameComplex[[i,2]]<>ToString[Im[eigensameComplex[[i,3]]],TraditionalForm]},{i,nComplex}]];

table = Grid[table,Frame-> All];
table
]


(* ::Input::Initialization:: *)
(*14/03/2019*)
PrintEigenfunctions[eqn_List,spectrum_,Nmax_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions,NumericalOptions}]]:= Block[{plot,equation,testlist,var,varstring,funclist,funcmax,eigenvectors,eigenfunctions,reeigenfunclist,imeigenfunclist,ulb = OptionValue[LowerBound],uub=OptionValue[UpperBound],prec = If[OptionValue[FunctionPrecision] == 0, Max[2*Nmax,$MinPrecision], OptionValue[FunctionPrecision]]},

equation = eqn//Expand;

(*find independent variable*)
testlist = Union[Cases[eqn,\[Phi]_Symbol[u_Symbol]-> u, \[Infinity]],Cases[eqn,Derivative[n_][\[Phi]_Symbol][u_Symbol]-> u,\[Infinity]]];
Which[Length[testlist]<1,ThrowAbort[GetEigenfunctions::novariables];,
Length[testlist]>1,ThrowAbort[GetEigenfunctions::manyvariables,testlist];,
True,var = testlist//Last;];

(*find dependent functions*)
testlist =  Union[Cases[eqn,\[Phi]_Symbol[var]->\[Phi],\[Infinity]],Cases[eqn,Derivative[n_][\[Phi]_Symbol][var]-> \[Phi],\[Infinity]]];
Which[Length[testlist]<Length[eqn],ThrowAbort[GetEigenfunctions::underdetermined,testlist,Length[eqn]];,
Length[testlist]>Length[eqn],ThrowAbort[GetEigenfunctions::overdetermined,testlist,Length[eqn]];,
True, funclist = testlist;];

funcmax = Length[funclist];
varstring=ToString[var];
Clear[var];

eigenfunctions = GetEigenfunctions[eqn,spectrum,Nmax,var,FilterRules[{opts},{Normalization,LowerBound,UpperBound,LBPower,UBPower,FunctionPrecision,CollocationGrid,FunctionCutoff,MaxIteration}]];

reeigenfunclist = Table[Tooltip[Re[eigenfunctions[[eigenindex,funcindex]]],ToString[SetPrecision[spectrum[[eigenindex]],5]]],{funcindex,funcmax},{eigenindex,Length[eigenfunctions]}];
imeigenfunclist = Table[Tooltip[Im[eigenfunctions[[eigenindex,funcindex]]],ToString[SetPrecision[spectrum[[eigenindex]],5]]],{funcindex,funcmax},{eigenindex,Length[eigenfunctions]}];

plot = Quiet[
Table[Plot[Evaluate@(reeigenfunclist[[funcindex]]~Join~imeigenfunclist[[funcindex]]),{var,ulb,uub},PlotRange->All,Axes->True,ImageSize->800,Frame->True,BaseStyle->{FontSize->24},LabelStyle->{FontSize->24},FrameLabel->{varstring,ToString[funclist[[funcindex]]]},PlotStyle->Table[{Blue},{Length[eigenfunctions]}]~Join~Table[{Red},{Length[eigenfunctions]}],WorkingPrecision->prec],{funcindex,funcmax}]
];

plot
]
PrintEigenfunctions[eqn_,spectrum_,Nmax_,opts:OptionsPattern[{EigenfunctionOptions,BoundaryOptions,NumericalOptions}]] := PrintEigenfunctions[{eqn},spectrum,Nmax,opts]


(* ::Input::Initialization:: *)
GetEigenfunctions::novariables= "No independent variables detected."
GetEigenfunctions::manyvariables = "Too many independent variables detected. `1` detected, expecting 1."
GetEigenfunctions::underdetermined = "Problem is underdetermined. `1` dependent variables detected, expecting `2`."
GetEigenfunctions::overdetermined = "Problem is overdetermined. `1` dependent variables detected, expecting `2`."


(* ::Section:: *)
(*Options, Error Handling, and Other Miscellaneous Functions*)


(* ::Input::Initialization:: *)
Attributes[CatchAbort]=Attributes[ThrowAbort]={HoldAll};
ThrowAbort[message_,args___]:= (Message[message,args];Throw[Abort[],$Failed])
CatchAbort[exp_] := Catch[exp,$Failed]


(* ::Input::Initialization:: *)
Options[BoundaryOptions]={LowerBound-> 0,UpperBound->1,LBPower-> 0,UBPower->0};
Options[NumericalOptions] = {FunctionPrecision-> 0,CollocationGrid->"Chebyschev"};
Options[ComparisonOptions] = {Cutoff-> 3,TestMagnitude->True};
Options[EigenCompOptions] = {L2Cutoff-> 3};
Options[EigenfunctionOptions]={Normalization->"UB",FinalAsymptotics->"Default",FunctionCutoff->12,MaxIteration->50};
Options[PrintOptions] = {FreqName -> "\[Omega]",NSpectrum-> -1,NEigenFunc->-1,ZeroCutoff->20};


(* ::Input::Initialization:: *)
FilterEigenvalues[spectrum_,tags_]:=Block[{unique=spectrum,conjugates1={},conjugates2={},testspec1,testspec2,\[Omega]test,prec,agreedDigits,filteredspec},

(* \[Omega] = -\[Omega]^* *)
testspec1 = Select[unique,Re[#]>0&];
testspec2 = Select[unique,Re[#]<= 0&];
For[i=1,i<= Length[testspec1]&&Length[testspec2]>0,i++,(\[Omega]test = First@MinimalBy[testspec2,Abs[#-(-Conjugate[testspec1[[i]]])]&];prec=Min[Precision/@{\[Omega]test,testspec1[[i]]}];
agreedDigits = (-Ceiling[Log10[Abs[\[Omega]test-(-Conjugate[testspec1[[i]]])]]/. Indeterminate -> -prec]+Floor[Log10[Abs[(\[Omega]test+(-Conjugate[testspec1[[i]]]))/2]]/. -Infinity->0]);
If[agreedDigits>=prec/2,conjugates1=Append[conjugates1,testspec1[[i]]];testspec2 = DeleteCases[testspec2,\[Omega]test];];
)];
unique = Complement[testspec1,conjugates1]~Join~testspec2;

(* \[Omega] = \[Omega]^* *)
testspec1 = Select[unique,Im[#]>0&];
testspec2 = Select[unique,Im[#]<= 0&];
For[i=1,i<= Length[testspec1]&&Length[testspec2]>0,i++,(\[Omega]test = First@MinimalBy[testspec2,Abs[#-(Conjugate[testspec1[[i]]])]&];prec=Min[Precision/@{\[Omega]test,testspec1[[i]]}];
agreedDigits = (-Ceiling[Log10[Abs[\[Omega]test-(Conjugate[testspec1[[i]]])]]/. Indeterminate -> -prec] + Floor[Log10[Abs[(\[Omega]test+(Conjugate[testspec1[[i]]]))/2]]/. -Infinity->0]);
If[agreedDigits>=prec/2,conjugates2=Append[conjugates2,testspec1[[i]]];testspec2 = DeleteCases[testspec2,\[Omega]test];];
)];
unique=Complement[testspec1,conjugates2]~Join~testspec2;

filteredspec = If[tags,Sort[Table[{"\[PlusMinus] ","",conjugates1[[i]]},{i,Length[conjugates1]}]~Join~Table[{"","\[PlusMinus] ",conjugates2[[i]]},{i,Length[conjugates2]}]~Join~Table[{"","",unique[[i]]},{i,Length[unique]}],Re[#1[[3]]]>Re[#2[[3]]]&],Sort[conjugates1~Join~conjugates2~Join~unique,Re[#1]>Re[#2]&]];
filteredspec
]


(* ::Input::Initialization:: *)
(*02/03/2019*)
CleanDigits[\[Omega]test_,\[Omega]ref_] := Block[{\[Omega]same,re\[Omega]test,im\[Omega]test,re\[Omega]same,im\[Omega]same,reave\[Omega] ,imave\[Omega] ,result,prec},
\[Omega]same = MinimalBy[\[Omega]ref,Abs[#-\[Omega]test]&]//First;
prec = Min[Precision/@{\[Omega]test,\[Omega]same}];
re\[Omega]test = Re[\[Omega]test];
im\[Omega]test = Im[\[Omega]test];
re\[Omega]same = Re[\[Omega]same];
im\[Omega]same = Im[\[Omega]same];

reave\[Omega] = Floor[Log10[Abs[(re\[Omega]test+re\[Omega]same)/2]]/. -Infinity-> 0];

imave\[Omega] = Floor[Log10[Abs[(im\[Omega]same+im\[Omega]test)/2]]/. -Infinity -> 0];
result= SetPrecision[re\[Omega]test,((-Floor[Ceiling[Log10[Abs[re\[Omega]test-re\[Omega]same]],0.2]] /. Indeterminate -> prec) +reave\[Omega])/._?Negative->0] + I SetPrecision[im\[Omega]test,((-Floor[Ceiling[Log10[Abs[im\[Omega]test-im\[Omega]same]],0.2]]/. Indeterminate -> prec)+imave\[Omega] )/._?Negative->0];
result
]


(* ::Input::Initialization:: *)
End[];
EndPackage[];
