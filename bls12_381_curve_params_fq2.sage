from sage.rings.finite_rings.integer_mod import square_root_mod_prime
# q for BLS12_381:
# q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
# our prime for BLS12_Cheon:
q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
F = GF(q)
K.<x> = PolynomialRing(F)
# nr for BLS12_381: (as per https://github.com/arkworks-rs/curves/blob/99831650f8021cb6a16481bac674420bc6c1a5a1/bls12_381/src/fields/fq2.rs#L13)
# there must not exist a sqrt of nr in F, i.e. this should be a non-quadratic residue.
# nr = F(-1)
# nr for BLS12_Cheon: (Tried -1,7,17,19,26,23,27,53 but Didnt work). Using 29.
nr = F(29)
try:
    # this is an efficient algo (Tonelli&Shanks) for computing square root mod p.
    # if no sqrt exists, then this returns an (incorrect) answer without checking
    nr_sqrt = square_root_mod_prime(Mod(nr,q))
    is_valid_Sqrt = (nr_sqrt**2 % q) == nr
    if is_valid_Sqrt:
        print(nr_sqrt, ": valid sqrt(",nr,") : BAD, THIS Shouldn't exist!!" )
    else:
        print("Sqrt(nr) doesn't exist, GOOD. This is a valid nqr!!")
except:
    print("Sqrt(nr) doesn't exist, GOOD. This is a valid nqr!!")
F2.<u> = F.extension(x ^ 2 - nr)
print(F2, " this is F2, is this a field? ", F2.is_field(), " (BAD IF NOT A FIELD!!)")
print(pow(nr, (1 - 1)/2, q), " this is FROBENIUS_COEFF_FP2_C1[0]")
print(pow(nr, (q - 1)/2, q), " this is FROBENIUS_COEFF_FP2_C1[1]", " q-coeff[1]: ", q - pow(nr, (q - 1)/2, q))
