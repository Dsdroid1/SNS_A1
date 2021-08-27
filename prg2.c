// 2.	Extended Euclidean algorithm to output 𝑥, 𝑦 when 𝑎, 𝑏 is given, such 𝑎𝑥 + 𝑏𝑦 = 𝑔𝑐𝑑(𝑎, 𝑏)
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Invalid Arguments");
    }
    else
    {
        mpz_t a, b, x, y, gcd_value;
        mpz_init(a);
        mpz_init(b);
        mpz_init(x);
        mpz_init(y);
        mpz_init(gcd_value);
        mpz_set_str(a, argv[1], 10);
        mpz_set_str(b, argv[2], 10);
        // Call the extended euclidean algorithm function
        extended_euclidean_algorithm(a, b, x, y, gcd_value);
        gmp_printf("%Zd %Zd", x, y);

        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(gcd_value);
    }
}