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
        mpz_t m, phi, no_of_primitive_roots, primitive_roots;
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
        gmp_printf("No.of primitive_roots:%Zd\n", no_of_primitive_roots);
        for (; mpz_cmp(rrsm_member, m) < 0 /*&& found == 0*/; mpz_add_ui(rrsm_member, rrsm_member, 1))
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
                    // mpz_set(generator, rrsm_member);
                    // gmp_printf("Smallest primitive_root of %Zd is %Zd", m, rrsm_member);
                    gmp_printf("%Zd ", rrsm_member);
                }
            }
        }
        // Now generate all primitive_root by raising the generator to powers which are coprime with phi(m)
        mpz_clear(power_mod_m);
        mpz_clear(generator);
        mpz_clear(rrsm_member);
        mpz_clear(gcd_value);
        mpz_clear(no_of_primitive_roots);
        mpz_clear(phi);
        mpz_clear(m);
    }
}