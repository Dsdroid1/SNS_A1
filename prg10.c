#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

int compare_for_qsort(const void *a, const void *b)
{
    mpz_t *mpz_a = (mpz_t *)a;
    mpz_t *mpz_b = (mpz_t *)b;
    return mpz_cmp(*mpz_a, *mpz_b);
}

void main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t m, phi, no_of_primitive_roots;
        mpz_init(m);
        mpz_init(phi);
        mpz_init(no_of_primitive_roots);
        mpz_set_str(m, argv[1], 10);
        // No.of primitive_roots is equal to phi(phi(m))
        eulers_totient_function(m, phi);
        eulers_totient_function(phi, no_of_primitive_roots);
        // Get all factors of phi(m)
        mpz_t *factors;
        int factor_count = 0;
        factors = factor_list(phi, &factor_count);
        // Now, for every member in RRSM_m, check if that raised to a factor of phi is never congruent to 1 mod m
        mpz_t rrsm_member, gcd_value;
        mpz_t power_mod_m, generator;
        mpz_init(power_mod_m);
        mpz_init(gcd_value);
        mpz_init(rrsm_member);
        mpz_init(generator);
        int found = 0;
        // DEBUG
        gmp_printf("No.of primitive_roots if atleast one exists:%Zd\n", no_of_primitive_roots);
        for (; mpz_cmp(rrsm_member, m) < 0 && found == 0; mpz_add_ui(rrsm_member, rrsm_member, 1))
        {
            gcd(rrsm_member, m, gcd_value);
            if (mpz_cmp_ui(gcd_value, 1) == 0)
            {
                // Belongs to rrsm, check the power mod m
                found = 1;
                // Get smallest primitive_root
                for (int i = 0; i < factor_count - 1 && found == 1; i++)
                {

                    fast_exponent_mod_m(rrsm_member, factors[i], m, power_mod_m);
                    if (mpz_cmp_ui(power_mod_m, 1) == 0)
                    {
                        // This cannot be a primitive_root
                        // gmp_printf("%Zd ", rrsm_member);
                        found = 0;
                    }
                }
                if (found == 1)
                {
                    mpz_set(generator, rrsm_member);
                    gmp_printf("Smallest primitive_root of %Zd is %Zd\n", m, rrsm_member);
                    // gmp_printf("%Zd ", rrsm_member);
                }
            }
        }

        if (found == 0)
        {
            printf("No primitive root exists!\n");
        }
        else
        {
            // Now generate all primitive_root by raising the generator to powers which are coprime with phi(m)
            mpz_t *primitive_roots, pr;
            primitive_roots = (mpz_t *)malloc(sizeof(mpz_t) * mpz_get_ui(no_of_primitive_roots));
            mpz_init(pr);
            int j = 0;
            for (int i = 1; i < mpz_get_ui(phi); i++)
            {
                mpz_set_ui(pr, i);
                // DEBUG
                gmp_printf("i=%Zd\n", pr);
                gcd(pr, phi, gcd_value);
                gmp_printf("GCD with phi:%Zd\n", gcd_value);
                if (mpz_cmp_ui(gcd_value, 1) == 0)
                {
                    // generator power pr is a primitive root
                    mpz_init(primitive_roots[j]);
                    fast_exponent_mod_m(generator, pr, m, primitive_roots[j]);
                    //DEBUG
                    gmp_printf("PR found as :%Zd\n", primitive_roots[j]);
                    j++;
                }
            }
            // Sort this using qsort
            qsort(primitive_roots, j, sizeof(mpz_t), compare_for_qsort);
            gmp_printf("No.of primitive_roots:%Zd\n", no_of_primitive_roots);
            for (int i = 0; i < j; i++)
            {
                gmp_printf("%Zd ", primitive_roots[i]);
            }
            mpz_clear(pr);
        }

        mpz_clear(power_mod_m);
        mpz_clear(generator);
        mpz_clear(rrsm_member);
        mpz_clear(gcd_value);
        mpz_clear(no_of_primitive_roots);
        mpz_clear(phi);
        mpz_clear(m);
    }
}