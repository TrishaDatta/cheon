q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
r = 1112136703329951464128597 
a = 0
b = 1
F = GF(q)
K.<x> = PolynomialRing(F)
nrF = 29
F2.<u> = F.extension(x ^ 2 - nrF)
K2.<y> = PolynomialRing(F2)
nrF2 = 7*u

newb = b*nrF2
# define a curve over F2: b' = b*nrF2
print("b' : ", newb)
Eq2 = EllipticCurve(F2, [0, 0, 0, 0, newb])

# G2 Cofactor:
cofactor = Eq2.cardinality()//r
print("card of Eq2 (E'):", Eq2.cardinality(), " is card multiple of r? [BAD if mod is NOT 0]", Eq2.cardinality()%r, "G2 cofactor: ", cofactor)
cofactor_inv = inverse_mod(Integer(cofactor), r) # pow(cofactor, -1, r)
print("G2 cofactor inverse: ", cofactor_inv, " Verifying: ", cofactor_inv*cofactor % r)
from textwrap import wrap
print(str(',\n').join(list(map(lambda x: '0x' + x[::-1], wrap(hex(cofactor)[2:][::-1], 16)))))

# Getting G2 generator: (both X and Y coordinates are elements of Fq2)
# Sampling a point on Eq2 (E')
# Randomly trying x, Need one s.t. res is a square. u+1 works!
x = u + 2
res = pow(x, 3)  + newb
print("for x:", x, " is res square: ", res.is_square())
y = res.sqrt()
Pq2 = Eq2(x, y)

G2_gen = Pq2*cofactor

print(G2_gen, G2_gen*r)
