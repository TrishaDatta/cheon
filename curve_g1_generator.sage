deg = 12
# BLS12_cheon prime
# p = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
# BLS12_cheon-Curve2
p = 5739212180072054691886037078074598258234828766015465682035977006377650233269727206694786492328072119776769214299497
# BLS12_Cheon r:
#r = 1112136703329951464128597
# BLS12_Cheon r (order of G1)
r = 1114157594638178892192613
K = GF(p)
# Initializing EC: a=0,b=1
a = 0
b = 4 # b = 1 for BLS12_cheon, but b=4 for BLS12_Cheon-Curve2
Ek = EllipticCurve(K, [a, b])
#random_point = Ek.random_point()
print(" ## card: ", Ek.cardinality(), " ## card%r (cardinality must be divisible by r): ", Ek.cardinality()%r , " ## order: " , Ek.order(), " ## order - card (order should be = card.): ", Ek.order() - Ek.cardinality() )
print("p - card: ", p - Ek.cardinality())
# Getting a point on the curve, for x=1:
# BLS12_Cheon-Curve2: x=2
x = 2
n = x**3 + b
from sage.rings.finite_rings.integer_mod import square_root_mod_prime
n_sqrt = square_root_mod_prime(Mod(n,p))
print(" Sqrt(", n,"):", n_sqrt, " is this correct? :", ((n_sqrt**2)%p == n) )
point = Ek(x,n_sqrt)
cofactor = Ek.cardinality()/r
print("G1 Cofactor: ", cofactor, xz)
cofactor_inv = inverse_mod(Integer(cofactor), r)
print("G1 Cofactor inverse: ", inverse_mod(Integer(cofactor), r), " Is inverse correct?", (cofactor_inv*cofactor)%r )

from textwrap import wrap
print(str(',\n').join(list(map(lambda x: '0x' + x[::-1], wrap(hex(cofactor)[2:][::-1], 16)))))
print("--- Cofactor inv:")
print(str(',\n').join(list(map(lambda x: '0x' + x[::-1], wrap(hex(cofactor_inv)[2:][::-1], 16)))))

# Multiplying any point on the curve by cofactor gives a point in G1
# this would be a generator if its not 0 (i.e. if its order != 1)
g1_gen = point*(cofactor)
print( "Point on EC: " , point, " Point in G1:" , g1_gen, " order: ", g1_gen.order(), " its a G1 generator if order != 1 ## is order=r? ", (g1_gen.order() == r) )
