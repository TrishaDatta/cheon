%%python3
def p(z):
    return (z - 1)^2 * (z^4 - z^2 + 1)/3 + z
def r(z):
    return z^4 - z^2 + 1

#x = 35403282275642611849
#rr = 1112136703329951464128597
# print(r(x), factor(r(x))) : 1570992498071942335289871609737960907960655536146175638441805600265637996638801 337 * 60169 * 1113997 * 1292953 * 5723582620213 * 8450424523622277446317 * 1112136703329951464128597
# print()
#a = 0
#b = 1
#F = GF(p(x))
#E = EllipticCurve(F, [a,b])
#print(r(x), E.cardinality(), E.cardinality()%(r(x)), int(r(x)).bit_length())
#print("bit-length of trisha r: " , rr.bit_length())

# find a z s.t. q=p(z) is prime, r=r(z) is composite & has an 80bit factor
def generate_prime_order(zbits):
    while True:
        z = randint(2^(zbits - 1), 2^zbits - 1)
        pz = int(p(z))
        # this condition is too stringent.
	#if (z% 72) != 7 and (z%72) != 16 and (z%72) != 31 and (z%72) != 64:
	    #print("TRYING z=", z, " Sad, z%72: ", z%72)
	    #continue
	if not is_prime(pz):
            continue
        rz = int(r(z))
        # check if there's a 80bit prime factor
        print("TRYING z = ", z, " got prime q: ", pz, " r: ", rz)
        r_facs = list(factor(rz))
        found = false
        for fac in r_facs:
            if fac[0].bit_length() == 60:
                found = true
                print("Yayy found a good factor! ", rz, r_facs, fac)
                break
        if not found:
            continue
    	print("FOUND ", z, pz, rz, r_facs)
    	K = GF(pz)
    	barr = [-1, 1, 2, 3, 4, 5, 6, 7]
    	found = false
    	for b in barr:
            curve = EllipticCurve(K, [0, b])
            card = curve.cardinality()
            if card % rz == 0:
                break
	        found = true
        if found:
    	    print("FOUND curve ", b)
    	    return curve
    	else:
	    continue

# will generate z's with bitlength 64bits
generate_prime_order(60)
