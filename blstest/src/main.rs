use ark_ec::{Group, pairing::Pairing, CurveGroup};
use ark_ff::{Field, One, MontFp};
use ark_bls12_cheon::Bls12Cheon;
use  ark_bls12_cheon::{
     G1Projective as G1, G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12
};
// use ark_bls12_381::Bls12_381 as Bls12Cheon;
// use ark_bls12_381::{G1Projective as G1,  G2Projective as G2, Fr as Fr, Fq, Fq2, Fq6, Fq12};

use ark_std::{Zero, UniformRand, ops::Mul, start_timer, end_timer};

use std::time::{Instant};
use std::ops::Div;

fn main() {
    println!("Hello, world!");

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

    /*
    // Making an array of 10^6 random G1 elements
    const arr_len: usize = 100;
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
    */

    // println!("Time elapsed in sampling {} random pts via [random_scalar_mult] is: {:?}", arr_len, tduration);

    // Testing degeneracy:
    let epgen = Bls12Cheon::pairing(g1_gen, g2_gen);
    println!(" e(g1,g2):  {}", epgen.0);

    // Testing pairing's bilinearity:
    let epsum = Bls12Cheon::pairing(g1_gen + g1_rand, g2_rand);
    let ep0 = Bls12Cheon::pairing(g1_gen, g2_rand);
    let ep1 = Bls12Cheon::pairing(g1_rand, g2_rand);
    println!("e(a0+a1,b): {}", epsum);
    println!("e(a0,b): {}", ep0);
    println!("e(a1,b): {}", ep1);
    println!("Testing bilinearity in G1: e(a0,b)+e(a1,b): {}, equal? {}", (ep0.0)*(ep1.0), (epsum.0 == (ep0.0 * ep1.0) ) );

    // Bilinearity Test-2:
    let epsum2 = Bls12Cheon::pairing(g1_rand, g2_rand + g2_rand2);
    let ep20 = Bls12Cheon::pairing(g1_rand, g2_rand);
    let ep21 = Bls12Cheon::pairing(g1_rand, g2_rand2);
    println!("Testing bilinearity in G2: equal? {}", ((epsum2.0 == (ep20.0 * ep21.0) )) );

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
