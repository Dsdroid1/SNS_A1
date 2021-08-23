#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t a, m;
        mpz_init(a);
        mpz_init(m);
        mpz_set_str(a, argv[1], 10);
        mpz_set_str(m, argv[2], 10);
        mpz_t x, y, gcd_value;
        mpz_init(x);
        mpz_init(y);
        mpz_init(gcd_value);

        extended_euclidean_algorithm(a, m, x, y, gcd_value);
        if (mpz_cmp_ui(gcd_value, 1) == 0)
        {
            // Inverse of a exists under mod m
            // Here, just verify if x needs to be normalized under (0,m)
            mpz_mod(x, x, m);
            gmp_printf("Y %Zd", x);
        }
        else
        {
            printf("N");
        }
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(gcd_value);
        mpz_clear(a);
        mpz_clear(m);
    }
}