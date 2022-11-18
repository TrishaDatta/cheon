use ark_ec::{Group, pairing::Pairing, CurveGroup};
use ark_ff::{Field, One, MontFp};
use ark_bls12_cheon::Bls12Cheon;
use  ark_bls12_cheon::{
     G1Projective as G1, G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12 as TargetField
};
// use ark_bls12_381::Bls12_381 as Bls12Cheon;
// use ark_bls12_381::{G1Projective as G1,  G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12 as TargetField};

use ark_std::{Zero, UniformRand, ops::Mul, start_timer, end_timer};

use std::time::{Instant};
use std::ops::Div;

use std::collections::HashMap;
use num_traits::pow;

/*
impl AsRef<[u64]> for u128
{
    fn as_ref(&self) -> &[u64]
    {
        // make an array of length 2
        let ans: [u64; 2] = [ self % (pow(2 as u128,64) ) , self/(pow(2 as u128,64) ) ];
        return ans;
    }
}
*/

// p is r (~80bit), d is a 30bit factor of r-1
// g_pair is a QuadExtField elem = e(g1_gen, g2_gen)
// g_1_pair is a QuadExtField elem = e(a*g1_gen, g2_gen) = e(g1_gen, g2_gen)^a
// g_d_pair is a QuadExtField/Fp12 elem = e((a^d1)*g1_gen, (a^d2)*g2_gen) = e(g1_gen, g2_gen)^(a^(d1+d2))
fn cheon_attack(p: u128, d: u32, g_pair: TargetField, g_1_pair: TargetField, g_d_pair: TargetField) -> i64
{
    /*
    let p_1 = p - 1;
    let fr_one = Fr::one();
    let zeta = Fr::from(2); // Fr::generator(); No generator() function, but we know 2 is a generator of Fr as per src/fr.rs
    // zeta^d can be done by pow() since d is ~30bits
    let zeta_hat = pow(zeta, d as usize);
    // Since d is a factor of p_1, p_1/d is actually a 51-bit whole number.
    // As per sage: m is a 26bit number
    let m = ( f64::sqrt( (p_1/(d as u128)) as f64 )).ceil() as u32;
    let m_hat =( ((p_1/(d as u128)) as f64)/(m as f64) ) .floor() as u32;
    println!("In cheon_attack: p_1: {}, zeta: {}, zeta_hat: {}, m: {}, m_hat: {}", p_1, zeta, zeta_hat, m, m_hat);

    // lookup dictionary: 
    // keys are points in target group, a quadExtField element [can stringify it if needed]
    let lookup1 = HashMap::new();
    let zeta_hat_inv = fr_one.div(Fr::from(zeta_hat));

    let mult = Fr::one();
    for u in 0..m
    {
        // TODO:1 convert mult to ?
        // mult might be 80bits since mult in Fr, so we cant use pow(.,.)
        // .pow requires mult to of class that implements AsRef<u64> 
        let pointu = g_d_pair.pow((mult as u64)); // pow(g_d_pair, mult as usize); // QuadExtField: has add & mul
        lookup1.insert(pointu, u);
        mult = mult.mul(zeta_hat_inv);
        if u % 100 == 7
        {
            println!("IN U Loop, u: {}, mult: {}, mult as usize: {}", u, mult, mult as usize);
        }
    }

    let k_0 : u128 = p - 1;
    let u_found : u32 = m+1;
    let zeta_hat_m = pow(zeta_hat, m as usize);
    let mult2 = Fr::one();
    for v in 0..m_hat
    {
        // TODO:2 
        let test = pow(g_pair, mult2 as usize);
        match lookup1.get(&test)
        {
            Some(&val) => {
                println!("YAYY!!! Key {} found in lookup1!! val: {}", test, val);
                // u_found is < m
                u_found = val;
                // k_0 < (m + p_1/d) < p_1
                k_0 = u_found as u128 + (m*v) as u128;
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
    let zeta_upside_down_hat = pow(zeta, (p_1/(d as u128)) as usize );

    let lookup2 = HashMap::new();
    for u_prime in 0..m_prime
    {
        // TODO: Need zeta inverse in Fr, raised to k0+uprime*p1/d.
        // Exp: zeta.inverse()^k0 * zeta_upside_down_hat.inverse()^u_prime
        // Exp can be 80-bit, so can't use pow().
        let exp = zeta.inverse().pow(k_0)  ;
        // let point = g_1_pair.pow();
        let point = pow(g_1_pair, zeta_upside_down_hat.pow( -1*u_prime)/zeta.pow(k_0) );
        lookup2.insert(point, u_prime);
    }

    // k_1 is < mprime + mprime*mprime_hat
    // i.e. < mprime + d < d+d.
    let k_1 : u32 = 2*d;
    // u_prime_found < m_prime, which is ~sqrt(d), so < d
    let u_prime_found : u64 = d;
    for v_prime in 0..m_prime_hat
    {
        let test = pow(g_pair, zeta_upside_down_hat.pow( m_prime*v_prime) );
        match lookup2.get(&test)
        {
            Some(&val) => {
                println!("YAYYY Key {} found in lookup2, val: {}", test, val);
                u_prime_found = val;
                k_1 = u_prime_found + m_prime * v_prime;
                break;
            },
            _ => println!("Key {} not found in lookup2 :(", test),
        }
    }

    if (k_0 == (p-1)) || (k_1 == 2*d)
    {
        println!("DAYUM, No k_0/k_1 FOUND!!!");
        return -1;
    }
    else
    {
        let alpha = k_0 + k_1*(p_1/(d as u128));
        println!("FOUND alpha!! {}", zeta.pow( alpha ) );
        return zeta.pow( alpha );
    }
    */
    return -1;
}

fn main() {
    println!("Hello, world!");
   
    // TESTING Fr:
    let mut alpha  = Fr::from(73); // this is the secret
    println!("alpha square_in_place: {}", alpha.square_in_place());
    let fr_11 = Fr::from(11);
    let fr_one = Fr::one();
    // let fq_11_mul_inv = Fq::from(4173972494597857957735299693145162369625330011647611405117074186456472896923437968505299267147688814383104883126907);
    let fr_11_m_inv = fr_one.div(fr_11); // .div(Fq::one());
    let fr_11_a_inv = Fr::from(-11);
    let fr_zero = Fr::zero();
    println!("Mult Identity in Fr: {}, Element in Fr: {}, Its inverse : {}, minv 11: {} , Fr(-11): {} , Fr(11)+Fr(-11): {}, alpha: {}, alpha^1: {}, alpha*alpha: {}", fr_one, fr_11, fr_11_m_inv, fr_11_m_inv*fr_11, fr_11_a_inv, (fr_11 + fr_11_a_inv) == fr_zero, alpha, alpha.pow(&[1 as u64]), alpha*alpha );

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
    const arr_len: usize = 71;
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

    // Testing degeneracy:
    let epgen = Bls12Cheon::pairing(g1_gen, g2_gen);
    println!("Testing degeneracy: e(g1,g2):  {}, Equal to One? {}", epgen.0, epgen.0 == TargetField::one());
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
    let r : u128 = 1114157594638178892192613;
    let d : u32 = 702047372; // this is a 30bit factor of (r-1)
    let d1 = 11726539;
    let d2 = d - d1;
    let alpha  = Fr::from(73); // this is the secret
    
    // TODO: 
    let alpha_d = alpha.pow(&[0 as u64]);
    println!("##### GENERATING Cheon Puzzle Inputs: r: {}, d: {}, d1: {}, d2: {}, alpha: {}, alpha_d: {}, one: {}, one.square_in_place(): {}", r, d, d1, d2, alpha, alpha_d, Fr::one(), Fr::one().square_in_place());

    /*
    let alpha_d1 = alpha.pow( d1);
    let alpha_d2 = alpha.pow( d2);
    let g1_alpha = g1_gen.mul(Fr::from(alpha)); // g*alpha or g^alpha, alpha_P1
    let g1_alpha_d1 = g1_gen.mul(Fr::from(alpha_d1));
    let g2_alpha_d2 = g2_gen.mul(Fr::from(alpha_d2));
    let g_pair = Bls12Cheon::pairing(g1_gen, g2_gen);
    let g_1_pair = Bls12Cheon::pairing(g1_alpha, g2_gen);
    let g_d_pair = Bls12Cheon::pairing(g1_alpha_d1, g2_alpha_d2);
    cheon_attack(r, d, g_pair.0, g_1_pair.0, g_d_pair.0 );
    */

    /*
    let mut g1_sum = G1::zero();

    let start =  Instant::now();
    let addn_time = start_timer!(|| format!("About to do {} CurveGrp additions", arr_len) );
    for i in 0..arr_len
    {
        g1_sum = g1_sum + g1_list_aff[i];
    }
    let duration = start.elapsed();
    end_timer!(addn_time);
    // println!("Time elapsed in {} additions via [CurveGrp] is {:?}", arr_len, duration);
    */
}
