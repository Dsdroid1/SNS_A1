// 7.	[Solutions of congruence] Given a, b and m, print whether solution to the congruence 𝑎𝑥 ≡ 𝑏(𝑚𝑜𝑑𝑚) exist Y/N. If Yes then output the number of solutions and all the solutions.
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
    if (argc != 4)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t a, x, b, m, gcd_value;
        mpz_init(a);
        mpz_init(x);
        mpz_init(gcd_value);
        mpz_init(b);
        mpz_init(m);

        mpz_set_str(a, argv[1], 10);
        mpz_set_str(b, argv[2], 10);
        mpz_set_str(m, argv[3], 10);

        // Solution exists only if gcd(a,m)|b
        gcd(a, m, gcd_value);
        mpz_t remainder;
        mpz_init(remainder);
        mpz_mod(remainder, b, gcd_value);
        if (mpz_cmp_ui(remainder, 0) == 0)
        {
            // Solution exists
            // Check about alpha and beta
            // b = beta*g, m= mu*g, a= alpha*g where g is gcd(a,m)
            // gcd(alpha,mu)=1
            // The equation reduces to x = beta* alpha_inverse mod mu
            mpz_t beta, mu, alpha;
            mpz_init(beta);
            mpz_init(mu);
            mpz_init(alpha);

            mpz_divexact(beta, b, gcd_value);
            mpz_divexact(alpha, a, gcd_value);
            mpz_divexact(mu, m, gcd_value);

            mpz_t alpha_inverse, solutions, i;
            mpz_init(alpha_inverse);
            mpz_init(solutions);

            mpz_init(i);
            // Here gcd(mu,alpha)=1=> inverse exists
            inverse_under_mod_m(alpha, mu, alpha_inverse);
            printf("Y ");
            gmp_printf("%Zd\n", gcd_value);
            // The first solution is beta*alpha_inverse
            mpz_mul(solutions, beta, alpha_inverse);
            mpz_mod(solutions, solutions, mu);
            mpz_set_ui(i, 1);
            /* Print solution only if less than INT_MAX */
            if (mpz_cmp_ui(solutions, INT_MAX) <= 0)
            {
                gmp_printf("%Zd ", solutions);
            }

            for (; mpz_cmp(i, gcd_value) < 0 && /*Print only if solutions is less than INT_MAX*/ mpz_cmp_ui(solutions, INT_MAX) <= 0; mpz_add_ui(i, i, 1))
            {
                mpz_add(solutions, solutions, mu);
                gmp_printf("%Zd ", solutions);
            }

            mpz_clear(i);
            mpz_clear(alpha_inverse);
            mpz_clear(solutions);

            mpz_clear(beta);
            mpz_clear(alpha);
            mpz_clear(mu);
        }
        else
        {
            // No solution can exist
            printf("N");
        }

        mpz_clear(a);
        mpz_clear(x);
        mpz_clear(gcd_value);
        mpz_clear(b);
        mpz_clear(m);
    }
}