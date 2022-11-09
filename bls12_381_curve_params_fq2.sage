# q for BLS12_381:
#q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
# our prime for BLS12_Cheon:
q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
F = GF(q)
K.<x> = PolynomialRing(F)
# nr for BLS12_381: (as per https://github.com/arkworks-rs/curves/blob/99831650f8021cb6a16481bac674420bc6c1a5a1/bls12_381/src/fields/fq2.rs#L13)
# nr = F(-1)
# nr for BLS12_Cheon:
nr = F(-1)
F2.<u> = F.extension(x ^ 2 - nr)
print(F2, " this is F2")
print((7 * u).nth_root(3))
print(pow(nr, (1 - 1)/2, q), " this is FROBENIUS_COEFF_FP2_C1[0]")
print(pow(nr, (q - 1)/2, q), " this is FROBENIUS_COEFF_FP2_C1[1]")
