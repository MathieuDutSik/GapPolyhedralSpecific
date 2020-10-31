CheckLLL_cppImplementation:=function(GramMat)
  local FileIn, FileOut, output, TheCommand, recLLL, res, TheRemainder, TheTrans;
  FileIn:=Filename(POLYHEDRAL_tmpdir, "LAT_lll.in");
  FileOut:=Filename(POLYHEDRAL_tmpdir, "LAT_lll.out");
  output:=OutputTextFile(FileIn, true);
  CPP_WriteMatrix(output, GramMat);
  CloseStream(output);
  TheCommand:=Concatenation(FileLattLll, " ", FileIn, " ", FileOut);
  Exec(TheCommand);
  recLLL:=ReadAsFunction(FileOut)();
  RemoveFileIfExist(FileIn);
  RemoveFileIfExist(FileOut);
  #
  res:=LLLReducedGramMat(GramMat);
  TheRemainder:=res.remainder;
  TheTrans:=res.transformation;
  if TheRemainder<>recLLL.GramMat or TheTrans<>recLLL.Pmat then
    Error("Find inconsistency for LLL between GAP and CPP versions (supposed to be the same)");
  fi;
end;


SAGE_NullspaceMat:=function(eMat)
  local FileIn, FileOut, output, IsFirst, eLine, eVal, nbRow, nbCol, TheRes, TheCommand;
  FileIn:=Filename(POLYHEDRAL_tmpdir, "Nullspace.sage");
  FileOut:=Filename(POLYHEDRAL_tmpdir, "Nullspace.output");
  output:=OutputTextFile(FileIn, true);
  nbRow:=Length(eMat);
  nbCol:=Length(eMat[1]);
  AppendTo(output, "A = MatrixSpace(RationalField(),", nbRow, ",", nbCol, ")");
  AppendTo(output, "([");
  IsFirst:=true;
  for eLine in eMat
  do
    for eVal in eLine
    do
      if IsFirst=false then
        AppendTo(output, ",");
      fi;
      IsFirst:=false;
      AppendTo(output, eVal);
    od;
  od;
  AppendTo(output, "])\n");
  AppendTo(output, "eKer = A.kernel()\n");
  AppendTo(output, "eKerB = eKer.basis()\n");
  AppendTo(output, "dim=len(eKerB)\n");
  AppendTo(output, "strO=\"return [\"\n");
  AppendTo(output, "for i in range(dim):\n");
  AppendTo(output, "    if (i>0):\n");
  AppendTo(output, "        strO += \",\"\n");
  AppendTo(output, "    strO += \"[\"\n");
  AppendTo(output, "    for j in range(", nbRow, "):\n");
  AppendTo(output, "        if (j>0):\n");
  AppendTo(output, "            strO += \",\"\n");
  AppendTo(output, "        strO += str(eKerB[i][j])\n");
  AppendTo(output, "    strO += \"]\"\n");
  AppendTo(output, "strO += \"];\"\n");
  AppendTo(output, "print strO\n");
  CloseStream(output);
  #
  TheCommand:=Concatenation("sage ", FileIn, " > ", FileOut);
  Exec(TheCommand);
  TheRes:=ReadAsFunction(FileOut)();
  RemoveFileIfExist(FileOut);
  RemoveFileIfExist(FileIn);
  return List(TheRes, RemoveFraction);
end;


NullspaceMatMultiProg:=function(eMat)
  local nbRowMin, nbColMin;
  nbRowMin:=700;
  nbColMin:=700;
  if Length(eMat) > nbRowMin and Length(eMat[1]) > nbColMin then
    return SAGE_NullspaceMat(eMat);
  fi;
  return NullspaceMat(eMat);
end;


Pari_WriteMatrix:=function(output, TheMat)
  local nbLine, iLine, nbCol, iCol;
  AppendTo(output, "[");
  nbLine:=Length(TheMat);
  for iLine in [1..nbLine]
  do
    nbCol:=Length(TheMat[iLine]);
    for iCol in [1..nbCol]
    do
      if iCol>1 then
        AppendTo(output, ",");
      fi;
      AppendTo(output, TheMat[iLine][iCol]);
    od;
    if iLine<nbLine then
      AppendTo(output, ";\n");
    fi;
  od;
  AppendTo(output, "]");
end;



PariB_WriteMatrix:=function(output, TheMat)
  local nbLine, iLine, nbCol, iCol;
  AppendTo(output, "[");
  nbLine:=Length(TheMat);
  for iLine in [1..nbLine]
  do
    nbCol:=Length(TheMat[iLine]);
    AppendTo(output, "[");
    for iCol in [1..nbCol]
    do
      if iCol>1 then
        AppendTo(output, ",");
      fi;
      AppendTo(output, TheMat[iLine][iCol]);
    od;
    AppendTo(output, "]");
    if iLine<nbLine then
      AppendTo(output, ",\n");
    fi;
  od;
  AppendTo(output, "]");
end;



PariC_WriteListMatrix:=function(output, NameListMat, ListMat)
  local nbMat, iMat, eMat, nbLine, iLine, eLine, nbCol, iCol;
  nbMat:=Length(ListMat);
  for iMat in [1..nbMat]
  do
    eMat:=ListMat[iMat];
    nbLine:=Length(eMat);
    AppendTo(output, "U", iMat, "= [");
    for iLine in [1..nbLine]
    do
      AppendTo(output, "[");
      eLine:=eMat[iLine];
      nbCol:=Length(eLine);
      for iCol in [1..nbCol]
      do
        AppendTo(output, eLine[iCol]);
        if iCol<nbCol then
          AppendTo(output, ",");
        fi;
      od;
      AppendTo(output, "]");
      if iLine<nbLine then
        AppendTo(output, ",");
      fi;
    od;
    AppendTo(output, "]\n");
  od;
  AppendTo(output, NameListMat, " = [");
  for iMat in [1..nbMat]
  do
    AppendTo(output, "U", iMat);
    if iMat<nbMat then
      AppendTo(output, ",");
    fi;
  od;
  AppendTo(output, "]");
end;


Pari_WriteVector:=function(output, eVect)
  local iCol, nbCol;
  AppendTo(output, "[");
  nbCol:=Length(eVect);
  for iCol in [1..nbCol]
  do
    if iCol>1 then
      AppendTo(output, ",");
    fi;
    AppendTo(output, eVect[iCol]);
  od;
  AppendTo(output, "]");
end;


BuildSetMultipleB:=function(ListSiz)
  local ListS, eSiz;
  ListS:=[];
  for eSiz in ListSiz
  do
    Add(ListS, [1..eSiz]);
  od;
  return BuildSetMultiple(ListS);
end;


Matrix_GetNNZ_L1norm:=function(eMat)
  local eNorm_L1, eNNZ, eLine, eVal;
  eNorm_L1:=0;
  eNNZ:=0;
  for eLine in eMat
  do
    for eVal in eLine
    do
      eNorm_L1:=eNorm_L1 + AbsInt(eVal);
      if eVal<>0 then
	eNNZ:=eNNZ+1;
      fi;
    od;
  od;
  Print([eNorm_L1, eNNZ]);
  return [eNorm_L1, eNNZ];
end;


RandomScrambling:=function()
  local eNB, eSum, i;
  eNB:=GetDate() mod 10000;
  eSum:=0;
  for i in [1..eNB]
  do
    eSum:=eSum + Random([1..10]);
  od;
  Print("RandomScrambling eSum=", eSum, "\n");
end;

CanonicalizeNullspace:=function(eMat)
  local k, len, eSub, eSubMat;
  k:=Length(eMat);
  len:=Length(eMat[1]);
  eSub:=RandomSubset([1..len],k);
  eSubMat:=List(eMat, x->x{eSub});
  return List(Inverse(eSubMat) * eMat, RemoveFraction);
end;


WriteVectorComma:=function(outputarg, eLine)
  local eVal, len, i;
  len:=Length(eLine);
  WriteAll(outputarg, String(eLine[1]));
  for i in [2..len]
  do
    WriteAll(outputarg, Concatenation(", ", String(eLine[i])));
  od;
  WriteAll(outputarg, "\n");
end;


CddOutput:=function(Filename, Vect)
  local eVal, eV, output;
  output:=OutputTextFile(Filename, true);;
  for eV in Vect
  do
    for eVal in eV
    do
      AppendTo(output, " ", eVal);
    od;
    AppendTo(output, "\n");
  od;
  CloseStream(output);
end;


SaturationDeterminant:=function(ListVect)
  local n, eBasis, NSP, TotalBasis, ListSol, eRedMat;
  n:=Length(ListVect[1]);
  eBasis:=GetZbasis(ListVect);
  if RankMat(eBasis)=n then
    return DeterminantMat(eBasis);
  fi;
  NSP:=NullspaceIntMat(TransposedMat(ListVect));
  TotalBasis:=NullspaceIntMat(TransposedMat(NSP));
  ListSol:=List(ListVect, x->SolutionMat(TotalBasis, x));
  eRedMat:=BaseIntMat(ListSol);
  return DeterminantMat(eRedMat);
end;


ProjectingCone:=function(OBJ, SET)
  return List(OBJ, x->x{SET});
end;


GraphBySubgraph:=function(Graph, SubGraph)
  local nba, A, i, j, eVert;
  nba:=Graph.order+1;
  A:=NullGraph(Group(()), nba);
  for i in [1..nba-1]
  do
    for j in [1..nba-1]
    do
      if IsEdge(Graph,[i,j])=true then
	AddEdgeOrbit(A, [i,j]);
      fi;
    od;
  od;
  for eVert in SubGraph
  do
    AddEdgeOrbit(A, [eVert, nba]);
    AddEdgeOrbit(A, [nba, eVert]);
  od;
  return A;
end;


TestIsomorphySubGraph:=function(Graph, SubGraph1, SubGraph2)
  local nba, A, B, i, j, eVert;
  nba:=Graph.order+1;
  A:=NullGraph(Group(()), nba);
  B:=NullGraph(Group(()), nba);
  for i in [1..nba-1]
  do
    for j in [1..nba-1]
    do
      if IsEdge(Graph,[i,j])=true then
        AddEdgeOrbit(A, [i,j]);
        AddEdgeOrbit(B, [i,j]);
      fi;
    od;
  od;
  for eVert in SubGraph1
  do
    AddEdgeOrbit(A, [eVert, nba]);
    AddEdgeOrbit(A, [nba, eVert]);
  od;
  for eVert in SubGraph2
  do
    AddEdgeOrbit(B, [eVert, nba]);
    AddEdgeOrbit(B, [nba, eVert]);
  od;
  return IsIsomorphicGraph(A,B);
end;


GetIndexRealizintTheSort:=function(eList)
  local len, eListExt, SetListExt, TheValue, ListIdx, i;
  len:=Length(eList);
  eListExt:=List([1..len], x->[eList[x], x]);
  SetListExt:=Set(eListExt);
  #
  TheValue:="unset";
  if Position(eList, TheValue)<>fail then
    Error("An unclever bug");
  fi;
  ListIdx:=[];
  for i in [1..len]
  do
    if TheValue<>SetListExt[i][1] then
      TheValue:=SetListExt[i][1];
      Add(ListIdx, SetListExt[i][2]);
    fi;
  od;
  return ListIdx;
end;


GetAllPosition:=function(eVect, eVal)
  local len, TheRet, i;
  len:=Length(eVect);
  TheRet:=[];
  for i in [1..len]
  do
    if eVect[i]=eVal then
      Add(TheRet, i);
    fi;
  od;
  return TheRet;
end;


GetBlockPermutationGroup:=function(GRP, ListPart)
  local ListPermGens, eGen, eList, ePerm;
  ListPermGens:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=List(ListPart, x->Position(ListPart, OnSets(x, eGen)));
    ePerm:=PermList(eList);
    Add(ListPermGens, ePerm);
  od;
  return Group(ListPermGens);
end;


TransferToHexadecimal:=function(eVal)
  local eValWork, eRetStr, ListLetter, res;
  eValWork:=eVal;
  eRetStr:="";
  ListLetter:=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"];
  while(true)
  do
    if eValWork=0 then
      break;
    fi;
    res:=eValWork mod 16;
    eRetStr:=Concatenation(ListLetter[res+1], eRetStr);
    eValWork:=(eValWork - res) / 16;
  od;
  return eRetStr;
end;


TransferFromHexadecimal:=function(estr)
  local ListLetter, len, ThePow, TheSum, i, eLetter, pos;
  ListLetter:=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"];
  len:=Length(estr);
  ThePow:=1;
  TheSum:=0;
  for i in Reversed([1..len])
  do
    eLetter:=estr{[i]};
    pos:=Position(ListLetter, eLetter);
    TheSum:=TheSum + ThePow*(pos-1);
    ThePow:=ThePow*16;
  od;
  return TheSum;
end;


GetFunctionalityStandardLattice:=function(n,d)
  local ListTrans, GetPosition;
  ListTrans:=Set(BuildSet(n, [0..d-1]));
  GetPosition:=function(eTrans)
    local TheRes;
    TheRes:=List(eTrans, x->x mod d);
    return Position(ListTrans, TheRes);
  end;
  return rec(ListTrans:=ListTrans, GetPosition:=GetPosition);
end;

# A n-simplex is a n-dimensional simplex defined by n+1 points
# v(1), v(2), ...., v(n+1).
# We define the differences w(1)=v(2) - v(1), ...., w(n)=v(n+1) - v(1).
# We make the summation over all coordinate subsets in order to get the volume.
# It is
# W = sum_S   (det(w)_S)^2
# The volume is expressed as
#    V = sqrt(W) / Factorial(n!)
# W is the entry returned by this function
FuncSquareMvolume:=function(SetVector)
  local m, Dim, LS, u, SQR, TrLS, Mat, eSet, l;
  m:=Length(SetVector)-1;
  Dim:=Length(SetVector[1]);
  LS:=[];
  for u in [2..Length(SetVector)]
  do
    AddSet(LS, SetVector[u]-SetVector[1]);
  od;
  SQR:=0;
  TrLS:=TransposedMat(LS);
  for eSet in Combinations([1..Dim], m)
  do
    Mat:=[];
    for u in eSet
    do
      Add(Mat, TrLS[u]);
    od;
    l:=DeterminantMat(Mat);
    SQR:=SQR+l*l;
  od;
  return SQR;
end;


SYMPOL_PrintMatrixCppCode:=function(output, TitleMat, eMat)
  local nbRow, nbCol, eRec, iRow, iCol;
  nbRow:=Length(eMat);
  nbCol:=Length(eMat[1]);
  eRec:=RemoveFractionMatrixPlusCoef(eMat);
  AppendTo(output, "template<typename T>\n");
  AppendTo(output, "MyMatrix<T> ", TitleMat, "()\n");
  AppendTo(output, "{\n");
  AppendTo(output, "  T eDen=", eRec.TheMult, ";\n");
  AppendTo(output, "  std::vector<std::vector<T>> ListListNum={");
  for iRow in [1..nbRow]
  do
    if iRow>1 then
      AppendTo(output, ",");
    fi;
    AppendTo(output, "{");
    for iCol in [1..nbCol]
    do
      if iCol>1 then
        AppendTo(output, ",");
      fi;
      AppendTo(output, eRec.TheMat[iRow][iCol]);
    od;
    AppendTo(output, "}");
  od;
  AppendTo(output, "};\n");
  AppendTo(output, "  int nbRow=", nbRow, ";\n");
  AppendTo(output, "  int nbCol=", nbCol, ";\n");
  AppendTo(output, "  MyMatrix<T> RetMat(nbRow, nbCol);\n");
  AppendTo(output, "  for (int iRow=0; iRow<nbRow; iRow++)\n");
  AppendTo(output, "    for (int iCol=0; iCol<nbCol; iCol++) {\n");
  AppendTo(output, "      T eVal=ListListNum[iRow][iCol];\n");
  AppendTo(output, "      RetMat(iRow, iCol)=eVal/eDen;\n");
  AppendTo(output, "    }\n");
  AppendTo(output, "  return RetMat;\n");
  AppendTo(output, "}\n");
end;


RandomPolytopePoint:=function(EXT)
  local siz, TheDim, eSumPt, eSum, eEXT, h, eInsidePt, ePrism;
  siz:=10;
  TheDim:=Length(EXT[1]);
  eSumPt:=ListWithIdenticalEntries(TheDim,0);
  eSum:=0;
  for eEXT in EXT
  do
    h:=Random([0..siz]);
    eSum:=eSum + h;
    eSumPt:=eSumPt + h*eEXT;
  od;
  eInsidePt:=eSumPt/eSum;
  return eInsidePt;
end;


OrderJordanHolderDecomposition:=function(Grp)
  local LN, ListQuotientSequenceSubgroup, GRP;
  ListQuotientSequenceSubgroup:=[];
  GRP:=Grp;
  while(Order(GRP)>1)
  do
    LN:=MaximalNormalSubgroups(GRP);
    Add(ListQuotientSequenceSubgroup, Order(GRP)/Order(LN[1]));
    GRP:=LN[1];
  od;
  return ListQuotientSequenceSubgroup;
end;


QuotientPermutationGroup:=function(Grp, NormalSubGroup)
  local O, Gens, ListNewGen, PList, eElt, eGen, eO, jLin, Quotient, CanonicSurjection;
  if IsNormal(Grp, NormalSubGroup)=false then
    return false;
  fi;
  O:=Orbits(NormalSubGroup, Grp, OnRight);
  Gens:=GeneratorsOfGroup(Grp);
  ListNewGen:=[];
  for eGen in Gens
  do
    PList:=[];
    for eO in O
    do
      eElt:=eO[1]*eGen;
      for jLin in [1..Length(O)]
      do
        if eElt in O[jLin] then
          Add(PList, jLin);
        fi;
      od;
    od;
    Add(ListNewGen, PermList(PList));
  od;
  Quotient:=Group(ListNewGen);
  CanonicSurjection:=GroupHomomorphismByImages(Grp, Quotient, Gens, ListNewGen);
  if Kernel(CanonicSurjection)<>NormalSubGroup then
    Error("We have inconsistency in the quotient computation");
  fi;
  return rec(Classes:=ShallowCopy(O), Quotient:=Quotient, CanonicSurjection:=CanonicSurjection);
end;


SecondReduceGroupActionPlusHom:=function(TheGroup, SmallSet)
  local SecondGRP, phi, ListGens, eGen, NewListGens;
  NewListGens:=[];
  ListGens:=GeneratorsOfGroup(TheGroup);
  for eGen in ListGens
  do
    Add(NewListGens, PermList(List(SmallSet, x->Position(SmallSet, OnPoints(x, eGen)))));
  od;
  SecondGRP:=PersoGroupPerm(NewListGens);
  phi:=GroupHomomorphismByImagesNC(TheGroup, SecondGRP, ListGens, NewListGens);
  return rec(SecondGRP:=SecondGRP, phi:=phi);
end;


HilbertMatrix:=function(n)
  local TheMat, i, j;
  TheMat:=NullMat(n,n);
  for i in [1..n]
  do
    for j in [1..n]
    do
      TheMat[i][j]:=1/(i+j-1);
    od;
  od;
  return TheMat;
end;
