#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "helpers_gmp.h"

void main()
{
    // // 1. Test for working of gcd functions
    // mpz_t a, b, gcd_val;
    // mpz_init(a);
    // mpz_init(b);
    // mpz_init(gcd_val);
    // char str[1024]; // To scan the big integer as a string
    // printf("1. Woking of gcd function:");
    // printf("\nEnter the 1st number(in base-10):");
    // scanf("%1023s", str); // To avoid buffer overflow
    // mpz_set_str(a, str, 10);
    // printf("\nEnter the 2nd number(in base-10):");
    // scanf("%1023s", str);
    // mpz_set_str(b, str, 10);
    // gcd(a, b, gcd_val);
    // gmp_printf("\nThe value gcd is:%Zd\n", gcd_val);
    // mpz_clear(a);
    // mpz_clear(b);
    // mpz_clear(gcd_val);

    // 2. Test for extended_euclidean_algorithm
    // mpz_t a, b, x, y, gcd_val;
    // mpz_init(a);
    // mpz_init(b);
    // mpz_init(x);
    // mpz_init(y);
    // mpz_init(gcd_val);
    // char str[1024]; // To scan the big integer as a string
    // printf("1. Woking of extended euclidean algorithm function:");
    // printf("\nEnter the 1st number(in base-10):");
    // scanf("%1023s", str); // To avoid buffer overflow
    // mpz_set_str(a, str, 10);
    // printf("\nEnter the 2nd number(in base-10):");
    // scanf("%1023s", str);
    // mpz_set_str(b, str, 10);
    // extended_euclidean_algorithm(a, b, x, y, gcd_val);
    // mpz_t ax, by, linear_combo;
    // mpz_init(ax);
    // mpz_init(by);
    // mpz_mul(ax, a, x);
    // mpz_mul(by, b, y);
    // mpz_add(linear_combo, ax, by);
    // gmp_printf("The gcd is %Zd\n", gcd_val);
    // gmp_printf("Linear combo: %Zd*(%Zd) + %Zd*(%Zd) = %Zd\n", a, x, b, y, linear_combo);
    // mpz_clear(ax);
    // mpz_clear(by);
    // mpz_clear(linear_combo);
    // mpz_clear(a);
    // mpz_clear(b);
    // mpz_clear(x);
    // mpz_clear(y);
    // mpz_clear(gcd_val);

    // // 3. Test for euler's totient function
    // mpz_t n, phi;
    // mpz_init(n);
    // mpz_init(phi);
    // char str[1024]; // To scan the big integer as a string
    // printf("1. Woking of euler's totient function:");
    // printf("\nEnter the number(in base-10):");
    // scanf("%1023s", str); // To avoid buffer overflow
    // mpz_set_str(n, str, 10);
    // eulers_totient_function(n, phi);
    // gmp_printf("The value of phi is :%Zd\n", phi);
    // mpz_clear(n);
    // mpz_clear(phi);

    // // 4. Test for CRT
    // int n;
    // printf("Enter n:");
    // scanf("%d", &n);
    // mpz_t *a, *m;
    // a = (mpz_t *)malloc(n * sizeof(mpz_t));
    // m = (mpz_t *)malloc(n * sizeof(mpz_t));
    // printf("Enter all a-m pairs(one by one): ");
    // char str[1024];
    // for (int i = 0; i < n; i++)
    // {
    //     mpz_init(a[i]);
    //     scanf("%1023s", str);
    //     mpz_set_str(a[i], str, 10);
    //     mpz_init(m[i]);
    //     scanf("%1023s", str);
    //     mpz_set_str(m[i], str, 10);
    // }
    // mpz_t x;
    // mpz_init(x);
    // chinese_remainder_theorem(n, m, a, x);
    // gmp_printf("%Zd", x);

    // 5. Test for fast_exponent_mod_m
    mpz_t base, exp, m, results;
    mpz_init(base);
    mpz_init(exp);
    mpz_init(m);
    mpz_init(results);
    char str[1024];
    printf("Enter base:");
    scanf("%1023s", str);
    mpz_set_str(base, str, 10);
    printf("Enter exp:");
    scanf("%1023s", str);
    mpz_set_str(exp, str, 10);
    printf("Enter m:");
    scanf("%1023s", str);
    mpz_set_str(m, str, 10);
    fast_exponent_mod_m(base, exp, m, results);
    gmp_printf("\n The exponent is %Zd", results);
    mpz_clear(results);
    mpz_clear(m);
    mpz_clear(exp);
    mpz_clear(base);
}