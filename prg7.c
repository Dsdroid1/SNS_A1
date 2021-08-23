#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

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

            // The first solution is beta*alpha_inverse
            mpz_mul(solutions, beta, alpha_inverse);
            mpz_mod(solutions, solutions, mu);
            mpz_set_ui(i, 1);
            gmp_printf("%Zd ", solutions);
            for (; mpz_cmp(i, gcd_value) < 0; mpz_add_ui(i, i, 1))
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