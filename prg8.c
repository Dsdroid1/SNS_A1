// 8.	[Solution to system of congruences] Given set of integers {(𝑎𝑖 , 𝑏𝑖 , 𝑚𝑖 )}𝑖=1 to 𝑛 , print whether there exist common solution 𝑥 which satisfy the system of congruences of the form 𝑎𝑖𝑥 ≡ 𝑏𝑖 (𝑚𝑜𝑑 𝑚𝑖 ). If exist then print all the solutions. Use user defined CRT function
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

int chinese_remainder_theorem(int n, mpz_t *m, mpz_t *a, mpz_t x)
{
    // Check if all mi,mj are relatively coprime when i!=j
    int i, j, exists = 1;
    mpz_t g, M;
    mpz_init(g);
    mpz_init(M);
    // mpz_set(M, m[0]);
    mpz_set_ui(M, 1);
    // for (i = 0; i < n - 1 && exists == 1; i++)
    // {
    //     mpz_mul(M, M, m[i + 1]);
    //     for (j = i + 1; j < n; j++)
    //     {
    //         gcd(m[i], m[j], g);
    //         if (mpz_cmp_ui(g, 1) != 0)
    //         {
    //             exists = 0;
    //         }
    //     }
    // }
    /* Optimized way for pairwise coprimes */
    mpz_t LCM;
    mpz_init(LCM);
    mpz_set_ui(LCM, 1);
    for (i = 0; i < n; i++)
    {
        mpz_mul(M, M, m[i]);
        gcd(LCM, m[i], g);
        mpz_mul(LCM, LCM, m[i]);
        mpz_divexact(LCM, LCM, g);
    }

    if (mpz_cmp(LCM, M) == 0)
    {
        exists = 1;
    }
    else
    {
        exists = 0;
    }
    mpz_clear(LCM);
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
                    gmp_printf("Y %Zd(mod %Zd)", x, M);
                }
                else
                {
                    printf("N");
                }
                mpz_clear(x);
            }
            // Remember to free memory for a,b,m;
            mpz_clear(M);
            for (int i = 0; i <= idx / 3; i++)
            {
                mpz_clear(a[i]);
                mpz_clear(b[i]);
                mpz_clear(m[i]);
            }
            free(a);
            free(b);
            free(m);
        }
    }
}