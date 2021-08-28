// 9.	[Order] Given a and m, print order of a under modulo m.
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

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

int order_a_mod_m(mpz_t a, mpz_t m, mpz_t h)
{
    mpz_t phi, gcd_value;
    mpz_init(gcd_value);
    int exists = 1;
    gcd(a, m, gcd_value);
    if (mpz_cmp_ui(gcd_value, 1) == 0)
    {
        // Get the phi of m
        mpz_init(phi);
        eulers_totient_function(m, phi);
        // DEBUG
        // gmp_printf("Phi of %Zd is %Zd\n", m, phi);
        // We need to only check all factors of phi for order equation
        mpz_t factor, root_n, remainder, order_expression_value;
        mpz_init(factor);
        mpz_init(root_n);
        mpz_init(remainder);
        mpz_init(order_expression_value);

        mpz_set_ui(factor, 1);
        mpz_sqrt(root_n, phi);
        int found = 0;
        // For factors from 1 to root_n
        for (; mpz_cmp(factor, root_n) <= 0 && found == 0; mpz_add_ui(factor, factor, 1))
        {

            mpz_mod(remainder, phi, factor);
            if (mpz_cmp_ui(remainder, 0) == 0)
            {
                // Factor found, check it
                // DEBUG
                // gmp_printf("Factor:%Zd\n", factor);
                fast_exponent_mod_m(a, factor, m, order_expression_value);
                // DEBUG
                // gmp_printf("Order expression value:%Zd\n", order_expression_value);
                if (mpz_cmp_ui(order_expression_value, 1) == 0)
                {
                    found = 1;
                    exists = 1;
                    mpz_set(h, factor);
                }
            }
        }
        // DEBUG
        // printf("First part done\n");
        mpz_set(factor, root_n);
        mpz_t inverse_factor;
        mpz_init(inverse_factor);
        // For factors from root_n+1 to n
        for (; mpz_cmp_ui(factor, 0) > 0 && found == 0; mpz_sub_ui(factor, factor, 1))
        {
            mpz_mod(remainder, phi, factor);
            if (mpz_cmp_ui(remainder, 0) == 0)
            {
                // Factor found, check it
                mpz_divexact(inverse_factor, phi, factor);
                // gmp_printf("%Zd ", inverse_factor);
                // DEBUG
                // gmp_printf("Factor:%Zd\n", inverse_factor);
                fast_exponent_mod_m(a, inverse_factor, m, order_expression_value);
                // DEBUG
                // gmp_printf("Order expression value:%Zd\n", order_expression_value);
                if (mpz_cmp_ui(order_expression_value, 1) == 0)
                {
                    found = 1;
                    exists = 1;
                    mpz_set(h, inverse_factor);
                }
            }
        }
        mpz_clear(inverse_factor);

        mpz_clear(order_expression_value);
        mpz_clear(factor);
        mpz_clear(root_n);
        mpz_clear(remainder);
        mpz_clear(phi);
    }
    else
    {
        // Order will not exist as gcd is not 1
        exists = 0;
    }

    mpz_clear(gcd_value);
    return exists;
}

void main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t a, m, h;
        mpz_init(a);
        mpz_init(m);
        mpz_init(h);

        mpz_set_str(a, argv[1], 10);
        mpz_set_str(m, argv[2], 10);

        int exists = order_a_mod_m(a, m, h);
        if (exists == 1)
        {
            gmp_printf("%Zd", h);
        }
        else
        {
            printf("N");
        }

        mpz_clear(a);
        mpz_clear(m);
        mpz_clear(h);
    }
}
