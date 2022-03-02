# Implementation of an attack on key exchange protocol using matrices over tropical algebra from 
# Grigoriev, Shpilrain, Tropical cryptography, http://arxiv.org/pdf/1301.1195.pdf.
# 
# Matvei Kotov <mkotov@stevens.edu>, Alexander Ushakov <aushakov@stevens.edu>

Read("tropical_algebraquasi.g");

# Runs our attack.
#n=size matrix, a=minimum entries,b=maximum entries, D=degree, pm=minimal coefficient poly, pM=maximal coefficient poly
TestAttack := function(n,D,pm, pM,ApplyAttack)
  local A,B, PP,DD, DB, DA, p1, p2, q1, q2, u, v, KA, KB, KC, AA,BB,attack_result;
  PP := GenerateRandomMatrix(n, Random(1,10));
  DD := GenerateRandomMatrix(n, Random(1,10));
  DA:=GenerateRandDiagonal(n,-Random(1,10),Random(1,10));
  DB:=GenerateRandDiagonal(n,-Random(1,10),Random(1,10));
  A:= ProductOfTwoMatricesMinPlus(DA,PP);
  B:=ProductOfTwoMatricesMinPlus(DB,DD);
  p1 := GenerateRandomPolynomial(D, pm, pM);
  p2 := GenerateRandomPolynomial(D, pm, pM);
  u := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(p1, A,D), ApplyPolynomialMinPlus(p2, B,D));
  q1 := GenerateRandomPolynomial(D, pm, pM);
  q2 := GenerateRandomPolynomial(D, pm, pM);
  v := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(q1, A,D), ApplyPolynomialMinPlus(q2, B,D));
  KA := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(p1, A,D), ProductOfTwoMatricesMinPlus(v, ApplyPolynomialMinPlus(p2, B,D)));
  KB := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(q1, A,D), ProductOfTwoMatricesMinPlus(u, ApplyPolynomialMinPlus(q2, B,D)));
  if KA <> KB then
    return false;
  fi;

  attack_result := ApplyAttack(A, B, u, D, pm);
  if attack_result = fail then
    return false;
  fi;
  KC := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(attack_result[1], A,D), 
      ProductOfTwoMatricesMinPlus(v, ApplyPolynomialMinPlus(attack_result[2], B,D)));
  if KA <> KC then
    return false;
  fi;
  return true;
end;


# Runs a set of tests.
TestSuite := function(n,D, pm, pM, ApplyAttack, numberOfTests)
  local st, et, i, ok, fl, maxCovers;
  st := Runtime();
  ok := 0;
  fl := 0;
  for i in [1..numberOfTests] do
    if TestAttack(n,D, pm, pM, ApplyAttack) then
      Print("OK\n");
      ok := ok + 1;
    else
      Print("FAIL\n");
      fl := fl + 1;
    fi;
  od;
  et := Runtime();
  Print(et - st, "\n");
  Print("OK: ", ok, "\n");
  Print("FAIL: ", fl, "\n");
end;
