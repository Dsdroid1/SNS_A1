// 6.	Given a and m, first print whether multiplicative inverse of 𝑎 (𝑚𝑜𝑑 𝑚) exist Y/N then output its inverse, if exist.
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

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