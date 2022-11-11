# is_primitive_root() is too slow in sage.
# More efficient method for finding generators:
def find_gen(p,phi,phi_factors):
	for i in range(2,60):
		gen = Integer(i)
		good_gen = True
		mod_dict = {}
		for phif in phi_factors:
			mod_dict[phif] = pow(gen, phi/phif, p )
			if mod_dict[phif] == 1:
				good_gen = False
		if good_gen:
			print("FOR ", p, phi, phi_factors, "FOUND Generator: ", i, gen, mod_dict) 
			return i,gen
# Generator for Fq:
# BLS12_cheon prime
q = 656356683693923479863117532358453498706487471262312031713236062035153349754850600248902117707113718656412873932748617
# BLS12_381 prime: Much faster
# p = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
Fq = GF(q)
phi_q = q - 1
# factorization done online, since sage was too slow
phi_q_factors = [2,3,7,123803,70244607689767087,3972859184419,107657349576566315882257393053833,350121160965644007076583724011565423320646682657]
gq_i, gq_gen = find_gen(q,phi_q, phi_q_factors)
print(Fq(gq_i).order())

# Generator for Fr:
r = 1112136703329951464128597
Fr = GF(r)
phi_r = r-1
print(factor(phi_r)) # 2^2 * 3^2 * 17 * 31 * 3347 * 21630391 * 809701559
phi_r_factors = [2,3,17,31,3347,21630391,809701559]
gr_i, gr_gen = find_gen(r,phi_r,phi_r_factors)
print(Fr(gr_i).order())
