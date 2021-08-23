#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        mpz_t num, factor, root_n, remainder;
        mpz_init(num);
        mpz_init(factor);
        mpz_init(root_n);
        mpz_init(remainder);
        mpz_set_str(num, argv[1], 10);
        mpz_set_ui(factor, 2);

        mpz_mod(remainder, num, factor);
        while (mpz_cmp_ui(remainder, 0) == 0)
        {
            gmp_printf("%Zd ", factor);
            mpz_divexact(num, num, factor);
            mpz_mod(remainder, num, factor);
        }
        mpz_sqrt(root_n, num);
        mpz_set_ui(factor, 3);
        for (; mpz_cmp(factor, root_n) <= 0; mpz_add_ui(factor, factor, 2))
        {
            mpz_mod(remainder, num, factor);
            while (mpz_cmp_ui(remainder, 0) == 0)
            {
                gmp_printf("%Zd ", factor);
                mpz_divexact(num, num, factor);
                mpz_sqrt(root_n, num);
                mpz_mod(remainder, num, factor);
            }
        }
        if (mpz_cmp_ui(num, 2) > 0)
        {
            gmp_printf("%Zd", num);
        }

        mpz_clear(num);
        mpz_clear(factor);
        mpz_clear(root_n);
        mpz_clear(remainder);
    }
}