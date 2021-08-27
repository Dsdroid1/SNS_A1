// 4.	[Reduced Residue System Modulo m] Given an integer m, output the RRSM_m set of integers. And also output the value of 𝜑(𝑚)
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t num, rrsm_size, i, gcd_value;
        mpz_init(num);
        mpz_init(rrsm_size);
        mpz_init(i);
        mpz_init(gcd_value);
        mpz_set_str(num, argv[1], 10);
        mpz_set_ui(i, 1);
        // 1 will always be a member
        mpz_set_ui(rrsm_size, 0);
        if (mpz_cmp_ui(num, 0) != 0)
        {
            for (; mpz_cmp(i, num) < 0; mpz_add_ui(i, i, 1))
            {
                gcd(num, i, gcd_value);
                if (mpz_cmp_ui(gcd_value, 1) == 0)
                {
                    // Belongs to RRSM_num
                    mpz_add_ui(rrsm_size, rrsm_size, 1);
                    gmp_printf("%Zd ", i);
                }
            }
            // gmp_printf("\nPhi = %Zd ", rrsm_size);
            gmp_printf("\n%Zd", rrsm_size);
        }
        mpz_clear(num);
        mpz_clear(rrsm_size);
        mpz_clear(i);
        mpz_clear(gcd_value);
    }
}