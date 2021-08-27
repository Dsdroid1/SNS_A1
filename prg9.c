// 9.	[Order] Given a and m, print order of a under modulo m.
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

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
