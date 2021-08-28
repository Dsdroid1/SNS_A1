// 10.	[Primitive roots] Given m, print the total number of primitive roots exist and print all the primitive roots, under modulo m.
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
// Print till INT_MAX specified
#include <limits.h>

int compare_for_qsort(const void *a, const void *b)
{
    mpz_t *mpz_a = (mpz_t *)a;
    mpz_t *mpz_b = (mpz_t *)b;
    return mpz_cmp(*mpz_a, *mpz_b);
}

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

// Function to calculate the size of RRSM_number, i.e the euler's totient function
void eulers_totient_function(mpz_t number, mpz_t phi)
{
    mpz_t n, factor, remainder, root_n, factor_power_multiplicity_minus_one, factor_minus_one;
    mpz_init(n);
    mpz_set(n, number);
    mpz_set_ui(phi, 1);
    mpz_init(factor);
    mpz_init(remainder);
    int factor_found = 0;
    mpz_init(root_n);
    mpz_set_ui(factor, 2);
    mpz_mod(remainder, n, factor);
    mpz_init(factor_power_multiplicity_minus_one);
    mpz_set_ui(factor_power_multiplicity_minus_one, 1);
    mpz_init(factor_minus_one);
    while (mpz_cmp_ui(remainder, 0) == 0)
    {

        if (factor_found == 1)
        {
            mpz_mul(factor_power_multiplicity_minus_one, factor_power_multiplicity_minus_one, factor);
        }
        factor_found = 1;
        mpz_divexact(n, n, factor);
        mpz_mod(remainder, n, factor);
    }
    if (factor_found)
    {
        mpz_mul(phi, phi, factor_power_multiplicity_minus_one);
        mpz_set(factor_minus_one, factor);
        mpz_sub_ui(factor_minus_one, factor_minus_one, 1);
        mpz_mul(phi, phi, factor_minus_one);
        factor_found = 0;
        mpz_set_ui(factor_power_multiplicity_minus_one, 1);
    }
    mpz_sqrt(root_n, n);
    mpz_set_ui(factor, 3);
    for (; mpz_cmp(factor, root_n) <= 0; mpz_add_ui(factor, factor, 2))
    {
        mpz_mod(remainder, n, factor);
        while (mpz_cmp_ui(remainder, 0) == 0)
        {

            if (factor_found == 1)
            {
                mpz_mul(factor_power_multiplicity_minus_one, factor_power_multiplicity_minus_one, factor);
            }
            factor_found = 1;
            mpz_divexact(n, n, factor);
            mpz_mod(remainder, n, factor);
        }
        if (factor_found)
        {
            mpz_mul(phi, phi, factor_power_multiplicity_minus_one);
            mpz_set(factor_minus_one, factor);
            mpz_sub_ui(factor_minus_one, factor_minus_one, 1);
            mpz_mul(phi, phi, factor_minus_one);
            factor_found = 0;
            mpz_set_ui(factor_power_multiplicity_minus_one, 1);
            mpz_sqrt(root_n, n);
        }
    }
    if (mpz_cmp_ui(n, 2) > 0)
    {
        mpz_set(factor, n);
        mpz_set(factor_minus_one, factor);
        mpz_sub_ui(factor_minus_one, factor_minus_one, 1);
        mpz_mul(phi, phi, factor_minus_one);
    }
    mpz_clear(root_n);
    mpz_clear(n);
    mpz_clear(factor);
    mpz_clear(factor_minus_one);
    mpz_clear(factor_power_multiplicity_minus_one);
    mpz_clear(remainder);
}

void fast_exponent_mod_m(mpz_t base, mpz_t original_power, mpz_t m, mpz_t result)
{
    // Works only for positive powers
    mpz_t remainder, base_power, two, power;
    mpz_init(remainder);
    mpz_init(two);
    mpz_init(base_power);
    mpz_init(power);
    mpz_set(power, original_power);

    mpz_set_ui(two, 2);
    mpz_set_ui(result, 1);

    // mpz_mod(remainder, power, two);
    // // DEBUG
    // gmp_printf("\nRemainder:%Zd", remainder);
    // if (mpz_cmp_ui(remainder, 1) == 0)
    // {
    //     mpz_add_ui(result, result, 1);
    //     // DEBUG
    //     gmp_printf("\nResult:%Zd ", result);
    // }
    // mpz_fdiv_q(power, power, two);
    // // DEBUG
    // gmp_printf("\nPower:%Zd", power);
    // mpz_set(base_power, base);
    mpz_mod(base_power, base, m);
    // // DEBUG
    // gmp_printf("\nBase raised value:%Zd", base_power);

    while (mpz_cmp_ui(power, 0) > 0)
    {
        mpz_mod(remainder, power, two);
        // // DEBUG
        // gmp_printf("\nRemainder:%Zd", remainder);
        if (mpz_cmp_ui(remainder, 1) == 0)
        {
            mpz_mul(result, result, base_power);
            mpz_mod(result, result, m);
            // // DEBUG
            // gmp_printf("\nResult:%Zd ", result);
        }
        mpz_mul(base_power, base_power, base_power);
        mpz_mod(base_power, base_power, m);
        mpz_fdiv_q(power, power, two);
        // // DEBUG
        // gmp_printf("\nPower:%Zd", power);
        // // DEBUG
        // gmp_printf("\nBase raised value:%Zd", base_power);
    }

    mpz_clear(remainder);
    mpz_clear(two);
    mpz_clear(base_power);
    mpz_clear(power);
}

mpz_t *factor_list(mpz_t n, long long int *no_of_factors)
{
    mpz_t *list;

    // Now get all factors of n itself
    mpz_t factor, root_n, remainder;
    mpz_init(factor);
    mpz_init(root_n);
    mpz_init(remainder);
    mpz_sqrt(root_n, n);
    mpz_set_ui(factor, 1);
    *no_of_factors = 0;
    mpz_t inverse_factor;
    mpz_init(inverse_factor);
    // First we get the count of factors, and then actually build the list
    // For count of factors
    for (; mpz_cmp(factor, root_n) <= 0; mpz_add_ui(factor, factor, 1))
    {
        mpz_mod(remainder, n, factor);
        if (mpz_cmp_ui(remainder, 0) == 0)
        {
            // Factor found, increase count
            if (mpz_cmp(factor, root_n) == 0)
            {
                mpz_divexact(inverse_factor, n, factor);
                if (mpz_cmp(inverse_factor, factor) == 0)
                {
                    *no_of_factors = *no_of_factors + 1;
                }
                else
                {
                    // factors are factor, n/factor
                    *no_of_factors = *no_of_factors + 2;
                }
            }
            else
            {
                // factors are factor, n/factor
                *no_of_factors = *no_of_factors + 2;
            }
        }
    }
    mpz_set_ui(factor, 1);
    int ptr = 0;

    list = (mpz_t *)malloc((*no_of_factors) * sizeof(mpz_t));
    // For factors from 1 to root_n
    for (; mpz_cmp(factor, root_n) <= 0; mpz_add_ui(factor, factor, 1))
    {
        mpz_mod(remainder, n, factor);
        if (mpz_cmp_ui(remainder, 0) == 0)
        {
            // Factor found, append it
            // gmp_printf("%Zd ", factor);
            mpz_init(list[ptr]);
            mpz_set(list[ptr], factor);
            ptr++;
        }
    }
    mpz_set(factor, root_n);

    // For factors from root_n+1 to n
    for (; mpz_cmp_ui(factor, 0) > 0; mpz_sub_ui(factor, factor, 1))
    {
        mpz_mod(remainder, n, factor);
        if (mpz_cmp_ui(remainder, 0) == 0)
        {
            // Factor found, append it
            mpz_divexact(inverse_factor, n, factor);
            if (mpz_cmp(inverse_factor, root_n) != 0)
            {
                mpz_init(list[ptr]);
                mpz_set(list[ptr], inverse_factor);
                ptr++;
            }
            // gmp_printf("%Zd ", inverse_factor);
        }
    }
    return list;
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
        long long int factor_count = 0;
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
        // gmp_printf("No.of primitive_roots if atleast one exists:%Zd\n", no_of_primitive_roots);
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
                    // gmp_printf("Smallest primitive_root of %Zd is %Zd\n", m, rrsm_member);
                    // gmp_printf("%Zd ", rrsm_member);
                }
            }
        }

        if (found == 0)
        {
            printf("N");
        }
        else
        {
            // Now generate all primitive_root by raising the generator to powers which are coprime with phi(m)
            mpz_t *primitive_roots, pr;
            primitive_roots = (mpz_t *)malloc(sizeof(mpz_t) * mpz_get_ui(no_of_primitive_roots));
            mpz_init(pr);
            long long int j = 0;
            for (long long int i = 1; i < mpz_get_ui(phi); i++)
            {
                mpz_set_ui(pr, i);
                // DEBUG
                // gmp_printf("i=%Zd\n", pr);
                gcd(pr, phi, gcd_value);
                // gmp_printf("GCD with phi:%Zd\n", gcd_value);
                if (mpz_cmp_ui(gcd_value, 1) == 0)
                {
                    // generator power pr is a primitive root
                    mpz_init(primitive_roots[j]);
                    fast_exponent_mod_m(generator, pr, m, primitive_roots[j]);
                    // //DEBUG
                    // gmp_printf("PR found as :%Zd\n", primitive_roots[j]);
                    j++;
                }
            }
            // Sort this using qsort
            qsort(primitive_roots, j, sizeof(mpz_t), compare_for_qsort);
            gmp_printf("%Zd\n", no_of_primitive_roots);
            for (int i = 0; i < j; i++)
            {
                /*Print till INT_MAX specified*/
                if (mpz_cmp_ui(primitive_roots[i], INT_MAX) > 0)
                {
                    break;
                }
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