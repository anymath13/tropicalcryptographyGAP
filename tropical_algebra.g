
# Generates a random matrix A (interval matrices)
GenerateRandomMatrix := function(n, mm)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return 0;
else
  return  Random(mm, 2*mm );
  fi;
  end));
end;

GenerateRandomMatrixW := function(n, mm,MM)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return Random(mm, MM );
else
  return  Random(mm, MM );
  fi;
  end));
end;


GenerateRandomMatrixLindePuente := function(n, mm,k1)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return (-1)*k1;
else
  return  Random(mm, mm );
  fi;
  end));
end;

GenerateRandDiagonal := function(n, mm,mM)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return Random(mm, mM );
else
  return infinity  ;
  fi;
  end));
end;

#Generate List Rational Number
rational1:=function(D) 
return (List([0..D],i->i/D));
end;
#Generate List
listrational:=function(D)
local N;
N:= Union(List([1..D],rational1));
return List([1..Length(N)], j->[j,N[j]]) ;
end;

# Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM].
GenerateRandomPolynomial := function(D, pm, pM)
  return List([1..Random(1,Length(listrational(D)))], i -> [i, Random(pm, pM)]);
end;

#LR=listrational;
PowerOfMatrixMinPlus:=function(A,t,D)
local LR;
LR:=listrational(D);
if t=0 then 
return IdentityMatrixMinPlus(Length(A));
elif t=1 then 
return A;
else
return PowersMinPlus(A,LR[t][2]);
fi;
end;


ElementaryMatrixMinPlus :=function(n,k,l)
return List([1..n], i -> List([1..n], function(j) 
      if i = k and j=l then 
        return 0; 
      else 
        return infinity; 
      fi;
    end));
end;


#Return generator
GenMinPlus:=function(n)
local k,i,j,a,b,result;
#find i,j such that k=n*i+j
i:=[];
j:=[];
for k in [1..n^2] do
  #i:=result[k];
  #j:=result[k];
  i[k]:=QuoInt(k,n);
  j[k]:=RemInt(k,n);
  #Print('i',i);
  #Print('j',j);
od;
#return [i,j];
#result:=[];
for a in [i] do
for b in [j] do
if a +1=b then
  return IdentityMatrixMinPlus(n);
else
  return ElementaryMatrixMinPlus(n,a,b);
fi;
od;
od;
#return result;
#od;
end;
#GenMinPlusb:=function(n)
#local k,j;
#find i,j such that k=n*i+j
#j:=[];
#j:=[];
#for k in [1..n^2] do
  #i:=result[k];
  #j:=result[k];
  #i[k]:=QuoInt(k,n);
  #j[k]:=RemInt(k,n);
  #Print('i',i);
  #Print('j',j);
#od;
#return j;
#end;
#return [j];
#result:=function(n)
#local result,i,j,a,b;
#i:=GenMinPlus(n);
#j:=GenMinPlusb(n);
#result:=[];
#for a in [i] do
#for b in [j] do
#if a[i] +1=b[i] then
 # return IdentityMatrixMinPlus(n);
#else
 # return ElementaryMatrixMinPlus(n,a,b);
#fi;
#od;
#od;
#return result;
#od;
#end;


nrparts:= function(n)
local np;
np:= function(n, m)
local i, res;
if n = 0 then
return 1;
fi;
res:= 0;
for i in [1..Minimum(n,m)] do
 res:= res + np(n-i, i);
od;
return res;
end;
return np(n,n);
end;

# Returns a \otimes b.
ProductOfTwoScalarMinPlus := function(a, b)
  if IsInfinity(a) or IsInfinity(b) then
    return infinity;
  else 
    return a + b;
  fi;
end;


# Returns A \otimes b.
ProductOfTwoMatricesMinPlus := function(A, B)
  local i, j, n, result; 
  n := Length(A);
  result := [];
  for i in [1..n] do
    Add(result, []);
    for j in [1..n] do
      Add(result[i], Minimum(List([1..n], k -> ProductOfTwoScalarMinPlus(A[i][k], B[k][j]))));
    od;
  od;
  return result;
end;


# Returns zero matrix of size nxn over min-plus algebra.
ZeroMatrixMinPlus := function(n)
  return List([1..n], i -> List([1..n], j -> infinity));
end;


# Returns ident matrix of size nxn over min-plus algebra.
IdentityMatrixMinPlus := function(n)
  return List([1..n], i -> List([1..n], function(j) 
      if i = j then 
        return 0; 
      else 
        return infinity; 
      fi;
    end));
end;




# Returns A + B.
SumOfTwoMatricesMaxPlus := function(A, B)
  return ListN(A, B, function(a, b) return ListN(a, b, function(x, y) return Minimum(x, y); end); end);
end;


# Returns s \otimes A.
ProductOfScalarAndMatrixMinPlus := function(s, A)
  return List(A, a -> List(a, x -> ProductOfTwoScalarMinPlus(s, x)));
end;


# Retruns p(A).
ApplyPolynomialMinPlus := function(p, A,D)
  if Length(p) = 0 then
    return ZeroMatrixMinPlus(Length(A));
  fi;
  return Iterated(List(p, c -> ProductOfScalarAndMatrixMinPlus(c[2], PowerOfMatrixMinPlus(A, c[1],D))), SumOfTwoMatricesMaxPlus);
end;


# Returns A - B. Some elements of the matrices can be infinity.
MinusMatrixFromMatrix := function(A, B)
  return ListN(A, B, function(a, b) return ListN(a, b, function (x, y)
      if y = infinity then
        return fail;
      elif x = infinity then
        return infinity;
      else  
        return x - y;
      fi;
  end); end);
end;
