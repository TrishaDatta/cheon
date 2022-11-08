deg = 12

#p = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
#K.<z12> = GF(p^deg)
#Ek = EllipticCurve(K, [0, 1])
#random_point = Ek.random_point()
#print('order: ', random_point.order())
#print("2 * rp", random_point)


# toy example
p = 727
q = 241
deg = 12
t = 5

#alpha is tau; here we are calculating the params given to the puzzlers
alpha = mod(36, q)
d = 16
alpha_d = alpha ^ d
d1 = 9
d2 = 7
alpha_d1 = alpha ^ d1
alpha_d2 = alpha ^ d2

#G1 will be defined over E1
E1 = EllipticCurve(GF(p), [0, 7])
#card1 = E1.cardinality()
#P1 = (card1 / q) * E1.random_point() 
#print(P1)
# P1 is a generator of G1
P1 = E1(245,431)

# G2 will be defined over Ek
K.<z12> = GF((p,deg))
Ek = EllipticCurve(K, [0, 7])
#cardk = Ek.cardinality()
#Q = (cardk / q^2) * Ek.random_point()
#print('q: ', Q)
#print(Q.order())


# This code was used to generate the generator of G2
#tr = Q
#x = Q[0]
#y = Q[1]
#for i in range(1, deg):
#    tr += Ek(x^(p^i), y^(p^i))
#print(tr)
#P2 = Q - tr
#print(P2.order())

# P2 is a generator of G2
P2 = Ek(622*z12^11 + 500*z12^10 + 610*z12^9 + 558*z12^8 + 396*z12^7 + 57*z12^6 + 301*z12^5 + 431*z12^3 + 699*z12^2 + 97*z12 + 695, 22*z12^11 + 294*z12^10 + 359*z12^9 + 565*z12^8 + 572*z12^7 + 302*z12^6 + 200*z12^5 + 540*z12^4 + 418*z12^3 + 80*z12^2 + 38*z12 + 127)

# We convert P1 to a point on Ek so we can do pairings
P1 = Ek(P1)

alpha_P1 = Integer(alpha) * P1
alpha_d1_P1 = Integer(alpha_d1) * P1
alpha_d1_P2 = Integer(alpha_d2) * P2

# Use pairings to calculate the params we will need for Cheon's attack
g = P1.tate_pairing(P2, q, deg)
g_1 = (alpha_P1).tate_pairing(P2, q, deg)
g_d = alpha_d1_P1.tate_pairing(alpha_d1_P2, q, deg)

# Just checking that pairings are working as expected
print(g ^ (alpha^d) == g_d)
print((g ^ alpha) == g_1)

# Run built-in discrete log algorithm; should be slower for much larger p
print(discrete_log(g_1, g))

# The rest of the program implements Cheon's attack; the variable names are 
# generally the same as the ones in the original paper

p_1 = q - 1
# primite_root gives us a generator of F_q
zeta = mod(primitive_root(q), q)
zeta_hat = zeta ^ d
m = ceil(sqrt(1.0 * p_1 / d))
m_hat = floor(p_1 / (m * d))
print('m: ', m)
print('m_hat', m_hat)


lookup_1 = {}

zeta_hat_inv = zeta_hat ^ -1
mult = mod(1, q)

for u in range(0, m):   
    point = g_d ^ Integer(mult)
    lookup_1[point] = u
    mult = mult * zeta_hat_inv
    
k_0 = None


zeta_hat_m = zeta_hat ^ m
mult = 1
for v in range(0, m_hat):
    test = g ^ Integer(mult)
    if test in lookup_1:
        u = lookup_1[test]
        k_0 = u + m * v
        break
    mult = mult * zeta_hat_m


m_prime = ceil(sqrt(1.0 * d))
m_prime_hat = floor(d / m_prime)

zeta_upside_down_hat = zeta ^ (p_1 / d)

lookup_2 = {}
for u_prime in range(0, m_prime):
    point = g_1 ^ (Integer((zeta ^ -k_0) * (zeta_upside_down_hat ^ -u_prime)))
    lookup_2[point] = u_prime

k_1 = None

for v_prime in range(0, m_prime_hat):
    test = g ^ (Integer(zeta_upside_down_hat ^ (m_prime * v_prime)))
    if test in lookup_2:
        u_prime = lookup_2[test]
        k_1 = u_prime + m_prime * v_prime
        break

if k_0 == None or k_1 == None:
    print("SAD")
else:
    print('alpha is: ', zeta ^ (k_0 + k_1 * (p_1/d)))


