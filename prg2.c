// 2.	Extended Euclidean algorithm to output 𝑥, 𝑦 when 𝑎, 𝑏 is given, such 𝑎𝑥 + 𝑏𝑦 = 𝑔𝑐𝑑(𝑎, 𝑏)
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

void main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Invalid Arguments");
    }
    else
    {
        mpz_t a, b, x, y, gcd_value;
        mpz_init(a);
        mpz_init(b);
        mpz_init(x);
        mpz_init(y);
        mpz_init(gcd_value);
        mpz_set_str(a, argv[1], 10);
        mpz_set_str(b, argv[2], 10);
        // Call the extended euclidean algorithm function
        extended_euclidean_algorithm(a, b, x, y, gcd_value);
        gmp_printf("%Zd %Zd", x, y);

        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(gcd_value);
    }
}