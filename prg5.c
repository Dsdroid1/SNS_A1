// 5.	Given 𝑎, 𝑥 and 𝑛, output 𝑎^𝑥(𝑚𝑜𝑑 𝑛). Use Fermat’s theorem concept. Also print intermediate equations while computing, in any readable form.
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

void main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t a, x, n;
        mpz_init(a);
        mpz_init(x);
        mpz_init(n);

        mpz_set_str(a, argv[1], 10);
        mpz_set_str(x, argv[2], 10);
        mpz_set_str(n, argv[3], 10);

        mpz_t phi_of_n, r, q, ans, gcd_value;
        mpz_init(phi_of_n);
        mpz_init(r);
        mpz_init(q);
        mpz_init(ans);
        mpz_init(gcd_value);
        // Only valid if gcd(a,m) is 1
        gcd(a, n, gcd_value);
        if (mpz_cmp_ui(gcd_value, 1) == 0)
        {
            eulers_totient_function(n, phi_of_n);
            // Represent x as q*phi(n)+r
            // We know a^phi(n) mod n is 1
            mpz_fdiv_qr(q, r, x, phi_of_n);
            // Hence calculate only a^r mod n
            fast_exponent_mod_m(a, x, n, ans);
            // gmp_printf("%Zd ^ %Zd (mod %Zd) = %Zd ^ (%Zd*[%Zd] + %Zd) (mod %Zd)", a, x, n, a, q, phi_of_n, r, n);
            // gmp_printf(" = [%Zd ^ (%Zd*[%Zd]) (mod %Zd) * %Zd ^ %Zd (mod %Zd)] (mod %Zd)", a, q, phi_of_n, n, a, r, n, n);
            // gmp_printf(" = 1 * %Zd = %Zd", ans, ans);
            gmp_printf("%Zd", ans);
        }
        else
        {
            // Fermats theorem is not applicable here
            //mpz_powm(ans, a, x, n);
            fast_exponent_mod_m(a, x, n, ans);
            gmp_printf("%Zd", ans);
        }
        mpz_clear(q);
        mpz_clear(r);
        mpz_clear(phi_of_n);
        mpz_clear(ans);
        mpz_clear(gcd_value);
        mpz_clear(a);
        mpz_clear(x);
        mpz_clear(n);
    }
}