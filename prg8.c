#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main(int argc, char *argv[])
{
    int n;
    mpz_t *a, *m, *b;
    if (argc <= 2)
    {
        printf("Invalid Arguments\n");
    }
    else
    {
        n = atoi(argv[1]);
        if (argc != 2 + n * 3)
        {
            printf("Invalid Arguments\n");
        }
        else
        {
            // We assume order of arguments as ai,bi,mi and such repeated for all n congruences
            a = (mpz_t *)malloc(n * sizeof(mpz_t));
            b = (mpz_t *)malloc(n * sizeof(mpz_t));
            m = (mpz_t *)malloc(n * sizeof(mpz_t));
            mpz_t ai_inverse_under_mod_mi, M;
            mpz_init(ai_inverse_under_mod_mi);
            mpz_init(M);
            mpz_set_ui(M, 1);
            int idx, exists = 1;
            for (int i = 0; i < n && exists == 1; i++)
            {
                mpz_init(a[i]);
                mpz_init(b[i]);
                mpz_init(m[i]);
                idx = 3 * i;
                mpz_set_str(a[i], argv[idx + 2], 10);
                mpz_set_str(b[i], argv[idx + 3], 10);
                mpz_set_str(m[i], argv[idx + 4], 10);
                mpz_mul(M, M, m[i]);
                // Check if gcd(ai,mi) = 1 for eqns to be reduced to CRT format,i.e inverse of ai to exist
                exists = inverse_under_mod_m(a[i], m[i], ai_inverse_under_mod_mi);
                if (exists == 1)
                {
                    mpz_mul(b[i], b[i], ai_inverse_under_mod_mi);
                    mpz_mod(b[i], b[i], m[i]);
                    // Modify eqn to form of CRT
                    // //DEBUG
                    // gmp_printf("Inverse of %Zd under %Zd is %Zd:\n", a[i], m[i], ai_inverse_under_mod_mi);
                    // gmp_printf("b[i] updates to :%Zd\n", b[i]);
                }
            }
            mpz_clear(ai_inverse_under_mod_mi);
            if (exists == 0)
            {
                printf("N");
            }
            else
            {
                mpz_t x;
                mpz_init(x);
                exists = chinese_remainder_theorem(n, m, b, x);
                if (exists == 1)
                {
                    gmp_printf("%Zd(mod %Zd)", x, M);
                }
                else
                {
                    printf("N");
                }
                mpz_clear(x);
            }
            // Remember to free memory for a,b,m;
            mpz_clear(M);
        }
    }
}