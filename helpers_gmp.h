#include <stdio.h>

#ifndef HELPERS_GMP
#define HELPERS_GMP

void gcd(mpz_t a, mpz_t b, mpz_t gcd_value);
void extended_euclidean_algorithm(mpz_t a_original, mpz_t b_original, mpz_t x, mpz_t y, mpz_t gcd_value);
void eulers_totient_function(mpz_t number, mpz_t phi);
int inverse_under_mod_m(mpz_t a, mpz_t m, mpz_t a_inverse);
int chinese_remainder_theorem(int n, mpz_t *m, mpz_t *a, mpz_t x);
void fast_exponent_mod_m(mpz_t base, mpz_t original_power, mpz_t m, mpz_t result);
int order_a_mod_m(mpz_t a, mpz_t m, mpz_t h);

#endif