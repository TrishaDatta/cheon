from sage.rings.finite_rings.integer_mod import square_root_mod_prime
# BLS12_381 prime
# q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
# BLS12_Cheon prime
q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
F = GF(q)
K.<x> = PolynomialRing(F)
# nrF for BLS12_381 extending Fq to Fq2 (Towering part1)
# nrF = -1
# nrF for BLS12_Cheon
nrF = 29
F2.<u> = F.extension(x ^ 2 - nrF)
# Checking that a. nrF is a non quadratic residue, b. F2 is a field
nrF_sqrt = square_root_mod_prime(Mod(nrF,q))
is_valid_sqrt = ((nrF_sqrt**2 % q) == nrF)
if is_valid_sqrt:
    print(nrF_sqrt, ": valid sqrt(",nrF,") : BAD, THIS Shouldn't exist!!" )
else:
    print("Sqrt(",nrF,") doesn't exist, GOOD. This is a valid nqr!!")
print(F2, " is F2 a field? ", F2.is_field(), " (BAD if not a field)")

K2.<y> = PolynomialRing(F2)
# nrF2 for BLS12_381 for extending Fq2 to Fq2 (Towering part2)
# nrF2 = u+1
# nrF2 for BLS12_Cheon for extending Fq2 to Fq12
nrF2 = 7*u
F12.<v> = F2.extension(y ^ 6 - nrF2)

# Checks for a. nrF2 should be a nonresidue (no 6th root), b. F12 should be a field
try:
    # This times out if a root exists, and is quick if no root exists.
    print(nrF2.nth_root(6), " 6thRoot(",nrF2,") exists: THIS IS BAD!")
except:
    print("NO 6th root(",nrF2,") : this is a GOOD nrF2!!")
print(F12, " is F12 a field? ", F12.is_field(), " (BAD if not field)")

print("FROBENIUS_COEFF_FP12_C1[0]: ", pow(nrF2, (1 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[1]: ", pow(nrF2, (q^1 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[2]: ", pow(nrF2, (q^2 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[3]: ", pow(nrF2, (q^3 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[4]: ", pow(nrF2, (q^4 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[5]: ", pow(nrF2, (q^5 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[6]: ", pow(nrF2, (q^6 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[7]: ", pow(nrF2, (q^7 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[8]: ", pow(nrF2, (q^8 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[9]: ", pow(nrF2, (q^9 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[10]: ", pow(nrF2, (q^10 - 1) / 6, q))
print("FROBENIUS_COEFF_FP12_C1[11]: ", pow(nrF2, (q^11 - 1) / 6, q))
