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
        mpz_t a_inv;
        mpz_init(a_inv);
        int exists;
        exists = inverse_under_mod_m(a, m, a_inv);

        if (exists == 1)
        {
            // Inverse of a exists under mod m
            gmp_printf("Y %Zd", a_inv);
        }
        else
        {
            printf("N");
        }

        mpz_clear(a_inv);
        mpz_clear(a);
        mpz_clear(m);
    }
}