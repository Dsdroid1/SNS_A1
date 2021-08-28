// 4.	[Reduced Residue System Modulo m] Given an integer m, output the RRSM_m set of integers. And also output the value of 𝜑(𝑚)
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
// Print till INT_MAX specified
#include <limits.h>

// Function to calculate the gcd of 2 integers, using euclid's formula
void gcd(mpz_t a, mpz_t b, mpz_t gcd_value)
{
    // Assumes that all mpz_t passed are already initialized
    if (mpz_cmp_ui(b, 0) == 0)
    {
        mpz_set(gcd_value, a);
    }
    else
    {
        mpz_t a_mod_b;
        mpz_init(a_mod_b);
        // Get the remainder of 'a' divided by 'b'
        mpz_tdiv_r(a_mod_b, a, b);
        gcd(b, a_mod_b, gcd_value);
        mpz_clear(a_mod_b);
    }
}

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
                if (mpz_cmp_ui(gcd_value, 1) == 0 && /*Printing values upto INT_MAX*/ mpz_cmp_ui(i, INT_MAX) <= 0)
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