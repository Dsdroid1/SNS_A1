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
            mpz_powm(ans, a, r, n);
            gmp_printf("%Zd ^ %Zd (mod %Zd) = %Zd ^ (%Zd*[%Zd] + %Zd) (mod %Zd)", a, x, n, a, q, phi_of_n, r, n);
            gmp_printf(" = [%Zd ^ (%Zd*[%Zd]) (mod %Zd) * %Zd ^ %Zd (mod %Zd)] (mod %Zd)", a, q, phi_of_n, n, a, r, n, n);
            gmp_printf(" = 1 * %Zd = %Zd", ans, ans);
            gmp_printf("\nFinal Ans: %Zd", ans);
        }
        else
        {
            // Fermats theorem is not applicable here
            //mpz_powm(ans, a, x, n);
            fast_exponent_mod_m(a, x, n, ans);
            gmp_printf("\nFinal Ans(Normal calculation as fermats theorem is not applicable here): %Zd", ans);
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