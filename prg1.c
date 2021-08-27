// 1.	Given 𝑛, {𝑚𝑖}, 𝑖=1 𝑛 integers, print all common divisors of (𝑚1, 𝑚2,… , 𝑚𝑛 ).
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    // Arguments to this program are assumed to be in order - prg1 n m1 m2 m3 ... mn
    int n;
    mpz_t num, gcd_value;
    mpz_init(num);
    mpz_init(gcd_value);
    if (argc <= 2)
    {
        printf("Invalid Arguments!\n");
    }
    else
    {
        n = atoi(argv[1]);
        // G Calculate the running gcd of all these 'n' numbers
        if (argc != n + 2)
        {
            printf("Invalid Arguments!\n");
        }
        else
        {
            mpz_set_str(gcd_value, argv[2], 10);
            int i = 1;
            for (; i < n; i++)
            {
                mpz_set_str(num, argv[i + 2], 10);
                // Get the gcd of this number and the running gcd value
                gcd(num, gcd_value, gcd_value);
            }
            // Now get all factors of gcd_value itself
            mpz_t factor, root_n, remainder;
            mpz_init(factor);
            mpz_init(root_n);
            mpz_init(remainder);
            mpz_sqrt(root_n, gcd_value);
            mpz_set_ui(factor, 1);
            // For factors from 1 to root_n
            // // DEBUG
            // gmp_printf("%Zd is root of gcd\n", root_n);
            for (; mpz_cmp(factor, root_n) <= 0; mpz_add_ui(factor, factor, 1))
            {
                mpz_mod(remainder, gcd_value, factor);
                // // DEBUG
                // gmp_printf("%Zd leads to %Zd remainder\n", factor, remainder);
                if (mpz_cmp_ui(remainder, 0) == 0)
                {
                    // Factor found, print it
                    gmp_printf("%Zd ", factor);
                }
            }
            mpz_set(factor, root_n);
            mpz_t inverse_factor;
            mpz_init(inverse_factor);
            // For factors from root_n+1 to n
            for (; mpz_cmp_ui(factor, 0) > 0; mpz_sub_ui(factor, factor, 1))
            {
                mpz_mod(remainder, gcd_value, factor);
                // // DEBUG
                // gmp_printf("%Zd leads to %Zd remainder\n", factor, remainder);
                if (mpz_cmp_ui(remainder, 0) == 0)
                {
                    // Factor found, print it
                    mpz_divexact(inverse_factor, gcd_value, factor);
                    if (mpz_cmp(inverse_factor, root_n) != 0)
                    {
                        gmp_printf("%Zd ", inverse_factor);
                    }
                }
            }
            mpz_clear(inverse_factor);
            mpz_clear(factor);
            mpz_clear(root_n);
            mpz_clear(remainder);
        }
    }
    mpz_clear(gcd_value);
    mpz_clear(num);
}