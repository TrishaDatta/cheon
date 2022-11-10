from sage.rings.finite_rings.integer_mod import square_root_mod_prime
# Verified coeffs/params of Fq6 used in BLS12_381 impl in arkworks
# q for BLS12_381:
# q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
# q for BLS12_cheon:
q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
F = GF(q)
K.<x> = PolynomialRing(F)
# nrF for BLS12_381 for extending F to Fq2
# nrF = -1
# nrF for BLS12_Cheon for extending F to Fq2
# Using 29, same as the one used in Fq2 script:
nrF = 29
nrF_sqrt = square_root_mod_prime(Mod(nrF,q))
is_valid_sqrt = ((nrF_sqrt**2 % q) == nrF)
if is_valid_sqrt:
    print(nrF_sqrt, ": valid sqrt(",nrF,") : BAD, THIS Shouldn't exist!!" )
else:
    print("Sqrt(",nrF,") doesn't exist, GOOD. This is a valid nqr!!")
F2.<u> = F.extension(x ^ 2 - nrF)
K2.<y> = PolynomialRing(F2)
print(F2, " is F2 a field? ", F2.is_field(), " (BAD if not a field)")

# nrF2 for BLS12_381 for extending F2 to F6, as per https://github.com/arkworks-rs/curves/blob/99831650f8021cb6a16481bac674420bc6c1a5a1/bls12_381/src/fields/fq6.rs#L14
# nrF2 = u+F(1)
# nrF2 for BLS12_Cheon for extending F2 to F6, (this is an elem in Fq2), u+F(1) didn't work
# 7u, 7u+F(1) works! Using 7u
nrF2 = 7*u
try:
    # This times out if a root exists, and is quick if no root exists.
    print(nrF2.nth_root(3), " CUBERoot(",nrF2,") exists: THIS IS BAD!")
except:
    print("NO cuberoot(",nrF2,") : this is a GOOD nrF2!!")

F6.<v> = F2.extension(y ^ 3 - nrF2)
print(F6, " is F6 a field? ", F6.is_field(), " (BAD if not field)")

# Printing the coefficients, each coefficient is an elem in Fq2, i.e. of form a+bu.
# these coeffs are initialized as Fq2::new(a,b) in arkworks.
print("FROBENIUS_COEFF_FP6_C1[0]: ", pow(nrF2, (1 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C1[1]: ", pow(nrF2, (q^1 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C1[2]: ", pow(nrF2, (q^2 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C1[3]: ", pow(nrF2, (q^3 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C1[4]: ", pow(nrF2, (q^4 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C1[5]: ", pow(nrF2, (q^5 - 1) / 3, q))

print("FROBENIUS_COEFF_FP6_C2[0]: ", pow(nrF2, 2 * (1 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C2[1]: ", pow(nrF2, 2 * (q^1 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C2[2]: ", pow(nrF2, 2 * (q^2 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C2[3]: ", pow(nrF2, 2 * (q^3 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C2[4]: ", pow(nrF2, 2 * (q^4 - 1) / 3, q))
print("FROBENIUS_COEFF_FP6_C2[5]: ", pow(nrF2, 2 * (q^5 - 1) / 3, q))
