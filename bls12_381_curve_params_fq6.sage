# Verified coeffs/params of Fq6 used in BLS12_381 impl in arkworks
q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
F = GF(q)
K.<x> = PolynomialRing(F)
nrF = -1
F2.<u> = F.extension(x ^ 2 - nrF)
K2.<y> = PolynomialRing(F2)
nrF2 = u+F(1)
F6.<v> = F2.extension(y ^ 3 - nrF2)
print("0", pow(nrF2, (1 - 1) / 3, q))
#pow(nrF2, (q - 1) / 3, q)
print("5", pow(nrF2, (q^5 - 1) / 3, q))

print("0", pow(nrF2, 2 * (1 - 1) / 3, q))
print("4", pow(nrF2, 2 * (q^4 - 1) / 3, q))
print("5", pow(nrF2, 2 * (q^5 - 1) / 3, q))
