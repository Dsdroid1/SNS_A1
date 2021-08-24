#include <stdio.h>
#include <gmp.h>
#include "helpers_gmp.h"

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

// Function to calculate the gcd of 2 numbers, and find 'x' and 'y' such that ax+by = gcd(a,b)
void extended_euclidean_algorithm(mpz_t a_original, mpz_t b_original, mpz_t x, mpz_t y, mpz_t gcd_value)
{
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);
    mpz_set(a, a_original);
    mpz_set(b, b_original);
    if (mpz_cmp(a, b) < 0)
    {
        extended_euclidean_algorithm(b, a, y, x, gcd_value);
    }
    else
    {
        mpz_t x1, x2, x_tmp;
        mpz_t y1, y2, y_tmp;
        mpz_t q, r;
        mpz_init(x1);
        mpz_init(x2);
        mpz_init(x_tmp);
        mpz_init(y1);
        mpz_init(y2);
        mpz_init(y_tmp);
        mpz_init(q);
        mpz_init(r);
        mpz_set_ui(x1, 1);
        mpz_set_ui(x2, 0);
        mpz_set_ui(y1, 0);
        mpz_set_ui(y2, 1);
        mpz_set_ui(x_tmp, 0);
        mpz_set_ui(y_tmp, 0);
        do
        {
            // Expressing a = qb + r => r = a - qb
            mpz_tdiv_qr(q, r, a, b);
            // Updating the coefficients to get the representation of 'r' in terms of 'a_original' and 'b_original'
            // x_tmp = x1 - q*x2
            mpz_mul(x_tmp, q, x2);
            mpz_neg(x_tmp, x_tmp);
            mpz_add(x_tmp, x1, x_tmp);
            // y_tmp = y1 - q*y2
            mpz_mul(y_tmp, q, y2);
            mpz_neg(y_tmp, y_tmp);
            mpz_add(y_tmp, y1, y_tmp);

            mpz_set(x1, x2);
            mpz_set(y1, y2);
            mpz_set(x2, x_tmp);
            mpz_set(y2, y_tmp);
            mpz_set(a, b);
            mpz_set(b, r);
        } while (mpz_cmp_ui(r, 0) != 0);
        mpz_set(x, x1);
        mpz_set(y, y1);
        mpz_set(gcd_value, a);
        // Clear all allocated resources
        mpz_clear(x1);
        mpz_clear(x2);
        mpz_clear(x_tmp);
        mpz_clear(y1);
        mpz_clear(y2);
        mpz_clear(y_tmp);
        mpz_clear(q);
        mpz_clear(r);
    }
    mpz_clear(a);
    mpz_clear(b);
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

int inverse_under_mod_m(mpz_t a, mpz_t m, mpz_t a_inverse)
{
    // Inverse only exists if gcd(a,m)=1
    int exists = 0;
    mpz_t g, x, y;
    mpz_init(g);
    mpz_init(x);
    mpz_init(y);

    extended_euclidean_algorithm(a, m, x, y, g);
    if (mpz_cmp_ui(g, 1) == 0)
    {
        // x is the inverse_under_mod_m
        mpz_mod(x, x, m);
        mpz_set(a_inverse, x);
        exists = 1;
    }
    else
    {
        // Inverse does not exist
        exists = 0;
    }

    mpz_clear(g);
    mpz_clear(x);
    mpz_clear(y);
    return exists;
}

int chinese_remainder_theorem(int n, mpz_t *m, mpz_t *a, mpz_t x)
{
    // Check if all mi,mj are relatively coprime when i!=j
    int i, j, exists = 1;
    mpz_t g, M;
    mpz_init(g);
    mpz_init(M);
    mpz_set(M, m[0]);
    for (i = 0; i < n - 1 && exists == 1; i++)
    {
        mpz_mul(M, M, m[i + 1]);
        for (j = i + 1; j < n; j++)
        {
            gcd(m[i], m[j], g);
            if (mpz_cmp_ui(g, 1) != 0)
            {
                exists = 0;
            }
        }
    }
    mpz_clear(g);
    if (exists == 1)
    {
        mpz_t b, M_by_mj, mult_term;
        // b is the inverse of M/mj under mod mj
        mpz_init(b);
        mpz_init(M_by_mj);
        mpz_init(mult_term);
        mpz_set_ui(x, 0);

        for (int j = 0; j < n; j++)
        {
            mpz_divexact(M_by_mj, M, m[j]);
            inverse_under_mod_m(M_by_mj, m[j], b);
            // //DEBUG
            // gmp_printf("M_by_mj is %Zd\n", M_by_mj);
            // gmp_printf("mj is %Zd\n", m[j]);
            // gmp_printf("ai is %Zd\n", a[i]);
            // gmp_printf("ai_inverse is %Zd\n", b);
            mpz_set(mult_term, M_by_mj);
            mpz_mul(mult_term, mult_term, b);
            mpz_mul(mult_term, mult_term, a[j]);
            mpz_add(x, x, mult_term);
        }

        mpz_clear(mult_term);

        mpz_clear(M_by_mj);
        mpz_clear(b);
        // Normalize the answer under mod M
        mpz_mod(x, x, M);
    }
    mpz_clear(M);
    return exists;
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
        gmp_printf("Phi of %Zd is %Zd\n", m, phi);
        // We need to only check all factors of phi for order equation
        mpz_t factor, root_n, remainder, order_expression_value;
        mpz_init(factor);
        mpz_init(root_n);
        mpz_init(remainder);
        mpz_init(order_expression_value);

        mpz_set_ui(factor, 1);
        mpz_sqrt(root_n, m);
        int found = 0;
        // For factors from 1 to root_n
        for (; mpz_cmp(factor, root_n) < 0 && found == 0; mpz_add_ui(factor, factor, 1))
        {

            mpz_mod(remainder, phi, factor);
            if (mpz_cmp_ui(remainder, 0) == 0)
            {
                // Factor found, check it
                // DEBUG
                gmp_printf("Factor:%Zd\n", factor);
                fast_exponent_mod_m(a, factor, m, order_expression_value);
                // DEBUG
                gmp_printf("Order expression value:%Zd\n", order_expression_value);
                if (mpz_cmp_ui(order_expression_value, 1) == 0)
                {
                    found = 1;
                    exists = 1;
                    mpz_set(h, factor);
                }
            }
        }
        // DEBUG
        printf("First part done\n");
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
                gmp_printf("Factor:%Zd\n", inverse_factor);
                fast_exponent_mod_m(a, inverse_factor, m, order_expression_value);
                // DEBUG
                gmp_printf("Order expression value:%Zd\n", order_expression_value);
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