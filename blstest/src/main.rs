use ark_ec::{Group, pairing::Pairing, CurveGroup};
use ark_ff::{Field, One, MontFp};
use ark_ff::biginteger::BigInteger128;
use ark_bls12_cheon::Bls12Cheon;
use  ark_bls12_cheon::{
     G1Projective as G1, G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12 as TargetField
};
// use ark_bls12_381::Bls12_381 as Bls12Cheon;
// use ark_bls12_381::{G1Projective as G1,  G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12 as TargetField};

use ark_std::{Zero, UniformRand, ops::Mul, start_timer, end_timer};

use std::time::{Instant, Duration};
use std::ops::Div;

use std::collections::HashMap;
use num_traits::pow;

// bitlen: bit length of exp, cant be more than 126 
fn pow_sp<S: Field>(p: S, exp: u128, bitlen: u32) -> S
{
    let mut res = S::one();
    // iterate over all bits in exp, in Big endian order.
    for b in (128 - bitlen - 2)..128
    {
        res = res*(res);
        let bi = exp & (1 << (127-b) ) > 0;
        // println!("index: {}, next bit: {}, exp: {}", b, bi, exp);
        if bi
        {
            res = res * p;
        }
    }
    return res;
}

// returns p^exp, p^(exp^2), ...., p^(exp^(n))
// assumes exp is 80bits:
fn pow_sp2<S: Field>(p: S, exp: u128, n: u64, k: u32) -> Vec<S>
{
    let mut vec = Vec::with_capacity(n as usize);

    // populate lookup table:
    // let k = ((n as f64).log2()) as u64;
    // 2^k sized array for 80/k buckets
    let m : u32 = 81;
    let num_bkt : u32 = (m/k);
    let p2k = pow(2,k as usize);

    println!("#####IN pow_sp2: kf: {} k: {}, m: {}, num_bkt: {}, p2k: {}", (n as f64).log2(), k, m, num_bkt, p2k);

    let pow_sp2_start = Instant::now();

    // let mut lookup_table: [ [S; 2^k]; num_bkt as usize];
    let mut lookup_table: Vec<Vec<S>> = Vec::new(); //with_capacity(num_bkt as usize);
    // All buckets should take ~equal time to fill.
    for bi in 0..num_bkt
    {
        let mut bkt_i: Vec<S> = Vec::new(); //with_capacity( p2k as usize );
        bkt_i.push(S::one()); // p^(0*2^(k*bi))
        if bi == 0
        {
            bkt_i.push(p);
        }
        else
        {
            bkt_i.push( pow_sp( lookup_table[(bi-1) as usize][1], p2k, k+1 as u32 ) ); // lookup_table[bi-1][1] ^ (2^k)
        }
        for j in 2..p2k
        {
            bkt_i.push( bkt_i[(j-1) as usize]*bkt_i[1] );
            
            if j % 7000 == 7
            {
                println!("Filling LookupTable, bkt#: {}, index: {}, Time since start: {:?}", bi, j, pow_sp2_start.elapsed());
            }
        }
        println!("FILLED bkt# {} (len: {}) in LookupTable!! Time: {:?}", bi, bkt_i.len(),  pow_sp2_start.elapsed());
        lookup_table.push(bkt_i);
    }
    println!("FILLED ALL BUCKETS in LookupTable!!! (len: {}) Time: {:?}", lookup_table.len(), pow_sp2_start.elapsed());

    // for verification, uisng pow_sp to get powers:
    let  mut pow_sp_ans:  Vec<S> = Vec::with_capacity( n as usize );
    let pow_start = p2k - 5; // 2*p2k
    let pow_end = p2k - 5 + n as u128; // 2*p2k + n
    for pi in pow_start..pow_end
    {
        pow_sp_ans.push(pow_sp(p, pi as u128, 81));
    }
    println!("FINIshed running pow_sp() for powers from {} to {}, time since start: {:?}", pow_start, pow_end,  pow_sp2_start.elapsed());

    for pi in pow_start..pow_end
    {
        let mut gpi = S::one();
        for bi in 0..num_bkt
        {
            // ind = pi & ( (p2k-1) << (27*bi) )
            let ind = pi as u128 & ( (p2k - 1) << (27*bi) );
            // get lookup[bi][ind]
            gpi = gpi * lookup_table[bi as usize][ind as usize];
            println!("Raising g to power {}, in bkt {}, ind: {}", pi, bi, ind);
       
        }

        if pi % 3 == 1
        {
            println!("DONE with {} exponentiations!!! (g^0...g^({}) ), Time: {:?}", pi, pi, pow_sp2_start.elapsed());
        }
        vec.push(gpi);
    }
    println!("FINISHED all exponents!!! from {} to {}, time since statr: {:?}", pow_start, pow_end,  pow_sp2_start.elapsed());
    for pi in 0..n
    {
        if vec[pi as usize] != pow_sp_ans[pi as usize]
        {
            println!("DAYUM, Answers for {} DONT MATCH b/w pow_sp and pow_sp2!!", pi);
        }
    }
    return vec;
}

fn average(numbers: &[Duration]) -> f64 {
    numbers.iter().sum::<Duration>().as_micros() as f64 / numbers.len() as f64
}

// p is r (~80bit), d is a 30bit factor of r-1
// g_pair is a QuadExtField elem = e(g1_gen, g2_gen)
// g_1_pair is a QuadExtField elem = e(a*g1_gen, g2_gen) = e(g1_gen, g2_gen)^a
// g_d_pair is a QuadExtField/Fp12 elem = e((a^d1)*g1_gen, (a^d2)*g2_gen) = e(g1_gen, g2_gen)^(a^(d1+d2))
fn cheon_attack(p: u128, d: u32, g_pair: TargetField, g_1_pair: TargetField, g_d_pair: TargetField) -> i64
{
    let attack_start = Instant::now();
    println!("IN cheon_attack, inputs: {},{},{},{},{}", p, d, g_pair, g_1_pair, g_d_pair);
    let p_1 = p - 1;
    let fr_one = Fr::one();
    // generator of Fr for bls12_cheon:
    let zeta = Fr::from(2); // Fr::generator(); No generator() function, but we know 2 is a generator of Fr as per src/fr.rs
    // zeta^d can be done by pow() since d is ~30bits
    let zeta_hat = pow_sp(zeta, d as u128, 32);
    // Since d is a factor of p_1, p_1/d is actually a 51-bit whole number.
    // As per sage: m is a 26bit number
    let m = ( f64::sqrt( (p_1/(d as u128)) as f64 )).ceil() as u32;
    let m_hat =( ((p_1/(d as u128)) as f64)/(m as f64) ) .floor() as u32;

    // lookup dictionary: 
    // keys are points in target group, a quadExtField element [can stringify it if needed]
    let mut lookup1 = HashMap::new();
    let zeta_hat_inv = fr_one.div(Fr::from(zeta_hat));
    println!("In cheon_attack: p_1: {}, zeta: {}, zeta_hat: {}, m: {}, m_hat: {}, zeta_hat_inv: {}, time: {:?}", p_1, zeta, zeta_hat, m, m_hat, zeta_hat_inv, attack_start.elapsed());
    
    // Nov18: All correct till this point.

    let mut mult = Fr::one();
    let mut timing_pow_arr: [Duration; 100] = [Duration::from_secs(0); 100]; 
    let mut timing_ins_arr: [Duration; 100] = [Duration::from_secs(0); 100]; 
    let mut timing_mul_arr: [Duration; 100] = [Duration::from_secs(0); 100];

    // choosing k=27, this is closest number to log2(m) which is also divisible by 81bit = len of
    // elems in Fr
    // testing with powers 0...9 :
    let gd_powers = pow_sp2(g_d_pair, bigInt_to_u128(zeta_hat_inv.into()), 100 as u64, 27 as u32 );

    for u in 0..m
    {
        let start = Instant::now();
        let pointu = pow_sp(g_d_pair, bigInt_to_u128(mult.into()), 81 ); // pow(g_d_pair, mult as usize); // QuadExtField: has add & mul
        timing_pow_arr[(u%100) as usize] = start.elapsed();
        
        lookup1.insert(pointu, u);
        timing_ins_arr[(u%100) as usize] = (start.elapsed()) - timing_pow_arr[(u%100) as usize];
        
        mult = mult.mul(zeta_hat_inv);
        timing_mul_arr[(u%100) as usize] = (start.elapsed()) - timing_ins_arr[(u%100) as usize] - timing_pow_arr[(u%100) as usize];
        if u % 100 == 0
        {
            println!("IN U Loop, u: {}, mult: {}, pointu: {}, time: {:?}", u, mult, pointu, attack_start.elapsed() );
           println!("IN U Loop, u: {}, timing arr's avg: Power: {}, Insert: {}, Mult: {}", u, average(&timing_pow_arr), average(&timing_ins_arr), average(&timing_mul_arr)); 
        }
    }

    let mut k_0 : u128 = p - 1;
    let mut u_found : u32 = m+1;
    let zeta_hat_m = pow_sp(zeta_hat, m as u128, 32);
    let mut mult2 = Fr::one();
    println!("Just before v-loop: k0: {}, ufound: {}, zeta_hat_m: {}", k_0, u_found, zeta_hat_m);
    for v in 0..m_hat
    {
        let test = pow_sp(g_pair, bigInt_to_u128(mult2.into()), 81 );
        match lookup1.get(&test)
        {
            Some(&val) => {
                // u_found is < m
                u_found = val;
                // k_0 < (m + p_1/d) < p_1
                k_0 = u_found as u128 + (m*v) as u128;
                println!("YAYY!!! Key {} found in lookup1!! val: {}, k_0: {}, time: {:?}", test, val, k_0, attack_start.elapsed());
                break;
            },
            _ => println!("Key {} not found in test :(", test),
        }
        mult2 = mult2.mul(Fr::from(zeta_hat_m));
    }

    if k_0 == (p-1)
    {
        println!("DAYUM, No k_0 FOUND!!!");
        return -1;
    }

    let m_prime = ( f64::sqrt(d as f64).ceil() ) as u64;
    let m_prime_hat = ( (d as f64) / (m_prime as f64)).floor() as u64;
    // can use pow() here since p1/d is a 51bit number.
    let zeta_upside_down_hat = pow_sp(zeta, (p_1/(d as u128)), 64 );
    println!("JUST before uprime loop: mprime {}, mprime_hat {}, zeta_upside_down_hat {}, time: {:?}", m_prime, m_prime_hat, zeta_upside_down_hat, attack_start.elapsed());

    let mut lookup2 = HashMap::new();
    for u_prime in 0..m_prime
    {
        let exp = pow_sp(Fr::from( fr_one.div(zeta)),k_0, 81) * pow_sp( fr_one.div( zeta_upside_down_hat), u_prime as u128, 64 ) ;
        let point = pow_sp(g_1_pair, bigInt_to_u128(exp.into()), 81 );
        lookup2.insert(point, u_prime);
        if (u_prime % 100 == 0)
        {
            println!("IN u_prime loop, uprime: {}, exp: {}, point: {}, time: {:?}", u_prime, exp, point, attack_start.elapsed() );
        }
    }

    // k_1 is < mprime + mprime*mprime_hat
    // i.e. < mprime + d < d+d.
    let mut k_1 : u64 = 2*d as u64;
    // u_prime_found < m_prime, which is ~sqrt(d), so < d
    let mut u_prime_found : u64 = d as u64;
    for v_prime in 0..m_prime_hat
    {
        let test = pow_sp(g_pair, bigInt_to_u128(pow_sp( zeta_upside_down_hat, (m_prime as u128)*(v_prime as u128) , 128 ).into()), 81 );
        match lookup2.get(&test)
        {
            Some(&val) => {
                u_prime_found = val;
                k_1 = u_prime_found + m_prime * v_prime;
                println!("YAYYY Key {} found in lookup2, val: {}, k_1: {}, time: {:?}", test, val, k_1, attack_start.elapsed());
                break;
            },
            _ => println!("Key {} not found in lookup2 :(", test),
        }
    }

    if (k_0 == (p-1)) || (k_1 == 2*d as u64)
    {
        println!("DAYUM, No k_0/k_1 FOUND!!!, k0: {}, k1: {}, total time: {:?}", k_0, k_1, attack_start.elapsed());
        return -1;
    }
    else
    {
        // k_0: u128, k_1: u64, p_1: u128, d: u32
        let alpha_log = k_0 + (k_1 as u128)*(p_1/(d as u128));
        let alpha = pow_sp(zeta, alpha_log, 126 );
        println!("FOUND alpha!! {} log: {}, k0: {}, k1: {}, total time: {:?}", alpha, alpha_log, k_0, k_1, attack_start.elapsed());
        return bigInt_to_u128(alpha.into()) as i64;
    }
}

fn bigInt_to_u128(bi: BigInteger128) -> u128
{
    let ans: u128 = (bi.0[0] as u128 + pow(2 as u128,64)*(bi.0[1] as u128) ) as u128;
    if ans % 100 == 7
    {
        println!("input: {}, u128: {}", bi, ans);
    }
    return ans;
}

fn main() {
    println!("Hello, world!");
   
    // TESTING Fr:
    let alpha  = Fr::from(73); // this is the secret
    println!("alpha: {}, alpha^1: {}, alpha^4: {}, alpja^19: {}, alpha^1 (default pow): {}", alpha, pow_sp(alpha,1 as u128, 16), pow_sp(alpha, 4 as u128, 16), pow_sp(alpha, 19 as u128, 32), alpha.pow(&[1 as u64]) );
    let fr_11 = Fr::from(11);
    let fr_one = Fr::one();
    // let fq_11_mul_inv = Fq::from(4173972494597857957735299693145162369625330011647611405117074186456472896923437968505299267147688814383104883126907);
    let fr_11_m_inv = fr_one.div(fr_11); // .div(Fq::one());
    let fr_11_a_inv = Fr::from(-11);
    let fr_zero = Fr::zero();
    println!("Mult Identity in Fr: {}, Element in Fr: {}, Its inverse : {}, minv 11: {} , Fr(-11): {} , Fr(11)+Fr(-11): {}, alpha: {}, alpha^1: {}, alpha*alpha: {}", fr_one, fr_11, fr_11_m_inv, fr_11_m_inv*fr_11, fr_11_a_inv, (fr_11 + fr_11_a_inv) == fr_zero, alpha, pow_sp(alpha,1, 2), alpha*alpha );

    // Testing fq square_in_place, pow funcs:
    let mut alpha1 = Fq::from(73);
    println!("a1:  sq_in_pl: {}", alpha1.square_in_place());
    // println!("TESTING Fq: alpha1: {}, alpha*alpha: {}, alpha1^1: {}, alpha1^10: {}, alpha sq in place: {}", alpha1, alpha1*alpha1, alpha1.pow(&[1 as u64]), alpha1.pow(&[10 as u64]), alpha1.square_in_place());

    let mut rng = ark_std::test_rng();
    let g1_rand0 = G1::rand(&mut rng);
    // this, when run twice, returned the same point
    let g1_rand = G1::rand(&mut rng); // prime_subgroup_generator() not supported
    // println!("point on G1: {}", g1_rand);

    let g1_gen = G1::generator();
    let g2_gen = G2::generator();

    let g2_rand = G2::rand(&mut rng);
    println!("point in G2: {}", g2_rand);

    let g1_rand2 = G1::rand(&mut rng);
    // println!("rand2 point in G2: {}", g1_rand2);

    // Testing group ops on G1:
    let c = g1_rand + g1_rand2;
    let d = g1_rand - g1_rand2;
    assert_eq!(c+d, g1_rand.double());
    let e = -g1_rand;
    assert_eq!(e+g1_rand,G1::zero());
    let scalar = Fr::rand(&mut rng);
    let e = c.mul(scalar);
    let f = e.mul(scalar.inverse().unwrap());
    assert_eq!(f, c);
    // affine repr:
    let g1_rand_aff = g1_rand.into_affine();
    assert_eq!(g1_rand, g1_rand_aff);
    // Testing done: All worked.
    
    // Testing group ops in G2:
    let g2_rand2 = G2::rand(&mut rng);
    let c = g2_rand2 + g2_rand;
    let d = g2_rand2 - g2_rand;
    assert_eq!(c+d, g2_rand2.double());
    assert_eq!(c-d, g2_rand.double());

    // Testing field ops in Fq:
    let fq_11 = Fq::from(11);
    let fq_one = Fq::one();
    // let fq_11_mul_inv = Fq::from(4173972494597857957735299693145162369625330011647611405117074186456472896923437968505299267147688814383104883126907);
    let fq_11_m_inv = fq_one.div(fq_11); // .div(Fq::one());
    let fq_11_a_inv = Fq::from(-11);
    let fq_zero = Fq::zero();
    println!("Mult Identity in Fq: {}, Element in Fq: {}, Its inverse : {}, Fq(-11): {} , Fq(11)+Fq(-11): {}", fq_one, fq_11, fq_11_m_inv, fq_11_a_inv, (fq_11 + fq_11_a_inv) == fq_zero);

    // Testing Fq2:
    // u:
    let fq2_zero = Fq2::new(Fq::ZERO, Fq::ZERO);
    let fq2_one = Fq2::new(Fq::ONE, Fq::ZERO);
    let fq2_a = Fq2::new(Fq::ZERO, Fq::ONE.double());
    // [bls12_cheon] additive inverse of a in Fq2 from Sage
    let fq2_a_ainv = Fq2::new(Fq::ZERO, MontFp!("5739212180072054691886037078074598258234828766015465682035977006377650233269727206694786492328072119776769214299495") );
    let fq2_a_minv = Fq2::new(Fq::ZERO, MontFp!("4366791876141780743826332559404585631265630582837854323288243374417777351400879396398207113727880960699715706532226"));
    println!("a: {}, a_ainv: {}, a+ainv: {}, a_minv: {}, a*minv: {}", fq2_a, fq2_a_ainv, fq2_a+fq2_a_ainv == fq2_zero, fq2_a_minv, fq2_a*fq2_a_minv);

    let fq2_b = Fq2::new(Fq::ONE, Fq::ONE.double());

    //  Testing pairing:
    let e1 = Bls12Cheon::pairing(g1_rand, g2_rand);
    let ml_result = Bls12Cheon::miller_loop(g1_rand, g2_rand);
    let e2 = Bls12Cheon::final_exponentiation(ml_result).unwrap();
    // println!("pairing: {} | {}", e1, e2);
    assert_eq!(e1, e2);
    // Testing done

    // Making an array of 10^6 random G1 elements
    const arr_len: usize = 1000;
    let mut g1_list: [ark_ec::short_weierstrass::Projective<ark_bls12_cheon::g1::Parameters>; arr_len] = [g1_rand; arr_len];
    // Making an equivalent list with AffineRepr:
    let mut g1_list_aff: [ark_ec::short_weierstrass::Affine<ark_bls12_cheon::g1::Parameters>; arr_len] = [g1_rand.into_affine() ; arr_len];

    let tstart = Instant::now();
    let list_time = start_timer!(|| "Generating List of random G1 points");
    for i in 0..arr_len
    {
        g1_list[i] = g1_list[i].mul(Fr::rand(&mut rng));
        g1_list_aff[i] = g1_list[i].into_affine();
    }
    let tduration = tstart.elapsed();
    end_timer!(list_time);

    println!("Time elapsed in sampling {} random pts via [random_scalar_mult] is: {:?}", arr_len, tduration);

    let mut g1_sum = G1::zero();

    let start =  Instant::now();
    let addn_time = start_timer!(|| format!("About to do {} CurveGrp additions", arr_len) );
    for i in 0..arr_len
    {
        g1_sum = g1_sum + g1_list_aff[i];
    }
    let duration = start.elapsed();
    end_timer!(addn_time);
    println!("Time elapsed in {} additions via [CurveGrp] is {:?}", arr_len, duration);

    // for bls12_cheon:
    let r : u128 = 1114157594638178892192613;
    // for bls12_cheon70:
    // let r : u128 = 765612544343321548849;

    // Testing degeneracy:
    let epgen = Bls12Cheon::pairing(g1_gen, g2_gen);
    println!("Testing degeneracy: e(g1,g2):  {}, Equal to One? {}, epgen^r == One? {}", epgen.0, epgen.0 == TargetField::one(), pow_sp(epgen.0, r, 81) == TargetField::one() );
    // println!("Generator of Fq12: {}", TargetField::multiplicative_generator());
    // epgen.pow(r) should be TargetField::one()

    // Testing pairing's bilinearity:
    let epsum = Bls12Cheon::pairing(g1_rand2 + g1_rand, g2_rand);
    let ep0 = Bls12Cheon::pairing(g1_rand, g2_rand);
    let ep1 = Bls12Cheon::pairing(g1_rand2, g2_rand);
    println!("e(a0+a1,b): {}", epsum);
    println!("e(a0,b): {}", ep0);
    println!("e(a1,b): {}", ep1);
    println!("Testing bilinearity in G1: e(a0,b)+e(a1,b): {}, equal? {}", (ep0.0)*(ep1.0), (epsum.0 == (ep0.0 * ep1.0) ) );

    // Bilinearity Test-2:
    let g2_rand3 = G2::rand(&mut rng);
    let epsum2 = Bls12Cheon::pairing(g1_rand, g2_rand + g2_rand3);
    let ep20 = Bls12Cheon::pairing(g1_rand, g2_rand);
    let ep21 = Bls12Cheon::pairing(g1_rand, g2_rand3);
    println!("Testing bilinearity in G2: equal? {}", ((epsum2.0 == (ep20.0 * ep21.0) )) );

    // CHeon algo:
    // p in Cheon paper is r for us.
    // this is the d for bls12_cheon:
    let d : u32 = 702047372; // this is a 30bit factor of (r-1)
    let d1 = 11726539;
    let d2 = d - d1;
    let alpha  = Fr::from(973); // this is the secret
    
    let alpha_d = pow_sp(alpha, d as u128, 32);
    let alp_bigint128 : BigInteger128 = alpha.into();
    println!("##### GENERATING Cheon Puzzle Inputs: r: {}, d: {}, d1: {}, d2: {}, alpha: {}, alpha_d: {}, one: {}, alpha.0: {}, bigInt128: {}", r, d, d1, d2, alpha, alpha_d, Fr::one(), alpha.0, alp_bigint128.0[0] as u128);

    let alpha_d1 = pow_sp(alpha, d1 as u128, 32);
    let alpha_d2 = pow_sp(alpha, d2 as u128, 32);
    let g1_alpha = g1_gen.mul(Fr::from(alpha)); // g*alpha or g^alpha, alpha_P1
    let g1_alpha_d1 = g1_gen.mul(Fr::from(alpha_d1));
    let g2_alpha_d2 = g2_gen.mul(Fr::from(alpha_d2));
    let g_pair = Bls12Cheon::pairing(g1_gen, g2_gen);
    let g_1_pair = Bls12Cheon::pairing(g1_alpha, g2_gen);
    let g_d_pair = Bls12Cheon::pairing(g1_alpha_d1, g2_alpha_d2);
    println!("TESTING correctness (bilinearity) of pairings: {}, {}", g_1_pair.0 == pow_sp(g_pair.0, bigInt_to_u128(alpha.into()), 81 ), g_d_pair.0 == pow_sp(g_pair.0, bigInt_to_u128(alpha_d.into()), 81 ) );
    
    // Testing exponent timing:
    let tstart = Instant::now();
    let mut gdpow = g_pair.0;
    let num_exp = 100;
    let start = 1783264826459163; // 2s for 100, for start=1,1783264
    // 25s for 1000 for start=1783264826459163
    for i in start..num_exp+start
    {
        gdpow = pow_sp(g_pair.0, i, 81);
    }
    println!("JUST DID {} exponentiations in TargetField, Time taken: {:?}", num_exp, tstart.elapsed());

    // Testing bigInt_to_u128, .into() timing:
    let tstart1 = Instant::now();
    let num_ops = 100;
    let mut ans : u128 = 0;
    for i in 0..num_ops
    {
        if i % 2 == 0
        {
            ans = bigInt_to_u128(alpha_d.into()); 
        }
        else
        {
            ans = bigInt_to_u128(alpha_d1.into());
        }
    }
    println!("JUST DID {} Fr->u128 conversions, Time taken: {:?}", num_ops, tstart1.elapsed());

    cheon_attack(r, d, g_pair.0, g_1_pair.0, g_d_pair.0 );

    }
