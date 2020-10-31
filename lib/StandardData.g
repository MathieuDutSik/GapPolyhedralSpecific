# set of points in Zn satisfying
# sum x_i = 0 (mod k)
LatticeLnk:=function(n,k)
  local TheBasis, i, eVect;
  TheBasis:=[];
  for i in [2..n]
  do
    eVect:=ListWithIdenticalEntries(n,0);
    eVect[i]:=1;
    eVect[1]:=-1;
    Add(TheBasis, eVect);
  od;
  eVect:=ListWithIdenticalEntries(n,0);
  eVect[1]:=k;
  Add(TheBasis, eVect);
  return rec(TheBasis:=TheBasis, TheGram:=RemoveFractionMatrix(TheBasis*TransposedMat(TheBasis)));
end;



ReadGramMatFile:=function(TheFile)
  local PMat, iLin, iCol, GramMat;
  PMat:=ReadVectorFile(TheFile);
  if Length(PMat[1])=1 then
    GramMat:=NullMat(Length(PMat), Length(PMat));
    for iLin in [1..Length(PMat)]
    do
      for iCol in [1..Length(PMat)]
      do
        if iCol<iLin then
          GramMat[iCol][iLin]:=PMat[iLin][iCol];
        else
          GramMat[iCol][iLin]:=PMat[iCol][iLin];
        fi;
      od;
    od;
    return GramMat;
  else
    return PMat;
  fi;
end;



FuncFormAnr:=function(n, r)
  local q, GramMat, i;
  if r=1 then
    return false;
  fi;
  q:=(n+1)/r;
  if IsInt(q)=false then
    return false;
  fi;
  GramMat:=NullMat(n,n);
  for i in [1..n-1]
  do
    GramMat[i][i]:=1;
  od;
  for i in [1..n-2]
  do
    GramMat[i+1][i]:=-1/2;
    GramMat[i][i+1]:=-1/2;
  od;
  GramMat[q][n]:=-1/2;
  GramMat[n][q]:=-1/2;
  GramMat[n][n]:=(1/2)*q*(1-(1/r));
  return RemoveFractionMatrix(GramMat);
end;




MobiusLadderGraph:=function(N)
  local GRA, i;
  GRA:=NullGraph(Group(()), 2*N);
  for i in [1..N]
  do
    AddEdgeOrbit(GRA, [2*i-1, 2*i]);
    AddEdgeOrbit(GRA, [2*i, 2*i-1]);
  od;
  for i in [2..N]
  do
    AddEdgeOrbit(GRA, [2*i, 2*(i-1)]);
    AddEdgeOrbit(GRA, [2*(i-1), 2*i]);
    AddEdgeOrbit(GRA, [2*i-1, 2*(i-1)-1]);
    AddEdgeOrbit(GRA, [2*(i-1)-1, 2*i-1]);
  od;
  AddEdgeOrbit(GRA, [2*N, 2 - 1]);
  AddEdgeOrbit(GRA, [2 - 1, 2*N]);
  AddEdgeOrbit(GRA, [2, 2*N-1]);
  AddEdgeOrbit(GRA, [2*N-1, 2]);
  return GRA;
end;

