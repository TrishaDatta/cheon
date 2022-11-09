q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
F = GF(q)
K.<x> = PolynomialRing(F)
nrF = -1
F2.<u> = F.extension(x ^ 2 - nrF)
K2.<y> = PolynomialRing(F2)
nrF2 = u+1
F12.<v> = F2.extension(y ^ 6 - nrF2)
print(pow(nrF2, (1 - 1) / 6, q),
      "##", pow(nrF2, (q^3 - 1) / 6, q), 
      "##",pow(nrF2, (q^6 - 1) / 6, q), 
      "##",pow(nrF2, (q^10 - 1) / 6, q))
