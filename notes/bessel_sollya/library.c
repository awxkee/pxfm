#include "library.h"
#include <mpfi.h>

__attribute__((visibility("default")))
int bessel_y1(mpfi_t result, mpfi_t x, int n)
{
    mpfr_t a, b;
    mpfr_init2(a, mpfi_get_prec(result));
    mpfr_init2(b, mpfi_get_prec(result));

    mpfi_get_left(a, x);
    mpfr_y1(a, a, GMP_RNDD);
    mpfi_get_right(b, x);
    mpfr_y1(b, b,GMP_RNDU);
    mpfi_interv_fr(result, a, b);

    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <mpfi.h>
#include <mpfr.h>

#include <mpfr.h>

void bessel_i1(mpfr_t result, const mpfr_t x, mpfr_prec_t prec, int max_terms, const mpfr_t epsilon)
{
    if (mpfr_sgn(x) == 0)
    {
        mpfr_set_d(result, 0, MPFR_RNDN);
        return;
    }
    mpfr_t sum, term, x_half, x_pow, k_fact, kp1_fact, denom;

    mpfr_inits2(prec, sum, term, x_half, x_pow, k_fact, kp1_fact, denom, (mpfr_ptr)0);
    mpfr_set_ui(sum, 0, MPFR_RNDN);

    // x_half = x / 2
    mpfr_div_ui(x_half, x, 2, MPFR_RNDN);

    // x_pow = x_half^1 = x/2
    mpfr_set(x_pow, x_half, MPFR_RNDN);

    // k_fact = 1, kp1_fact = 1
    mpfr_set_ui(k_fact, 1, MPFR_RNDN);
    mpfr_set_ui(kp1_fact, 1, MPFR_RNDN);

    for (int k = 0; k < max_terms; ++k)
    {
        if (k > 0)
        {
            mpfr_mul_ui(k_fact, k_fact, k, MPFR_RNDN); // k!
            mpfr_mul_ui(kp1_fact, kp1_fact, k + 1, MPFR_RNDN); // (k+1)!
            mpfr_mul(x_pow, x_pow, x_half, MPFR_RNDN);
            mpfr_mul(x_pow, x_pow, x_half, MPFR_RNDN); // (x/2)^{2k+1}
        }

        mpfr_mul(denom, k_fact, kp1_fact, MPFR_RNDN); // denom = k!(k+1)!
        mpfr_div(term, x_pow, denom, MPFR_RNDN); // term = x_pow / denom

        if (mpfr_cmpabs(term, epsilon) < 0)
            break;

        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    mpfr_set(result, sum, MPFR_RNDN);

    mpfr_clears(sum, term, x_half, x_pow, k_fact, kp1_fact, denom, (mpfr_ptr)0);
}

void besseli0(mpfr_t result, const mpfr_t x, mpfr_prec_t prec)
{
    mpfr_t term, sum, k_fact, x_half_pow, x_half;
    mpfr_inits2(prec, term, sum, k_fact, x_half_pow, x_half, (mpfr_ptr)0);

    mpfr_set_ui(sum, 1, MPFR_RNDN); // sum = 1
    mpfr_set_ui(k_fact, 1, MPFR_RNDN); // k! = 1
    mpfr_set_ui(x_half_pow, 1, MPFR_RNDN); // (x/2)^0 = 1

    mpfr_div_ui(x_half, x, 2, MPFR_RNDN); // x/2

    for (int k = 1; k < 1500; ++k)
    {
        // x_half_pow *= (x/2)^2
        mpfr_mul(x_half_pow, x_half_pow, x_half, MPFR_RNDN);
        mpfr_mul(x_half_pow, x_half_pow, x_half, MPFR_RNDN);

        // k_fact *= k
        mpfr_mul_ui(k_fact, k_fact, k, MPFR_RNDN);

        // term = x_half_pow / (k_fact)^2
        mpfr_mul(term, k_fact, k_fact, MPFR_RNDN);
        mpfr_div(term, x_half_pow, term, MPFR_RNDN);

        // sum += term
        mpfr_add(sum, sum, term, MPFR_RNDN);

        if (mpfr_cmp_d(term, 1e-40) < 0) break;
    }

    mpfr_set(result, sum, MPFR_RNDN);
    mpfr_clears(term, k_fact, sum, x_half_pow, x_half, (mpfr_ptr)0);
}

void compute_i0_approximant_asympt(mpfr_t result, const mpfr_t x, mpfr_prec_t prec)
{
    mpfr_t sqrt_x, exp_mx, bessi, recip, ones;
    mpfr_inits2(prec, sqrt_x, exp_mx, bessi, recip, ones, (mpfr_ptr)0);

    mpfr_set_d(ones, 1, MPFR_RNDN);

    // mpfr_div(recip, ones, x, MPFR_RNDN);
    mpfr_sqrt(sqrt_x, x, MPFR_RNDN);
    mpfr_neg(exp_mx, x, MPFR_RNDN);
    mpfr_exp(exp_mx, exp_mx, MPFR_RNDN);

    besseli0(bessi, x, prec);

    mpfr_mul(result, sqrt_x, exp_mx, MPFR_RNDN);
    mpfr_mul(result, result, bessi, MPFR_RNDN);

    mpfr_clears(sqrt_x, exp_mx, bessi, recip, ones, (mpfr_ptr)0);
}

void compute_i1_approximant_asympt_big(mpfr_t result, const mpfr_t x, mpfr_prec_t prec)
{
    mpfr_t sqrt_x, exp_mx, bessi, recip, ones, eps;
    mpfr_inits2(prec, sqrt_x, exp_mx, bessi, recip, ones, eps, (mpfr_ptr)0);

    mpfr_set_d(ones, 1, MPFR_RNDN);
    mpfr_set_d(eps, 1e-41, MPFR_RNDN);

    // mpfr_div(recip, ones, x, MPFR_RNDN);
    mpfr_sqrt(sqrt_x, x, MPFR_RNDN);
    mpfr_neg(exp_mx, x, MPFR_RNDN);
    mpfr_exp(exp_mx, exp_mx, MPFR_RNDN);

    bessel_i1(bessi, x, prec, 1500, eps);

    mpfr_mul(result, sqrt_x, exp_mx, MPFR_RNDN);
    mpfr_mul(result, result, bessi, MPFR_RNDN);

    mpfr_clears(sqrt_x, exp_mx, bessi, recip, ones, eps, (mpfr_ptr)0);
}

__attribute__((visibility("default")))
int bessel_i0_approximant(mpfi_t result, mpfi_t x, int n)
{
    mpfr_t a, b;
    mpfr_init2(a, mpfi_get_prec(result));
    mpfr_init2(b, mpfi_get_prec(result));

    mpfi_get_left(a, x);
    compute_i0_approximant_asympt(a, a, mpfi_get_prec(result));

    mpfi_get_right(b, x);
    compute_i0_approximant_asympt(b, b, mpfi_get_prec(result));
    mpfi_interv_fr(result, a, b);

    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
}

__attribute__((visibility("default")))
int bessel_i0(mpfi_t result, mpfi_t x, int n)
{
    mpfr_t a, b;
    mpfr_init2(a, mpfi_get_prec(result));
    mpfr_init2(b, mpfi_get_prec(result));

    mpfi_get_left(a, x);
    besseli0(a, a, mpfi_get_prec(result));

    mpfi_get_right(b, x);
    besseli0(b, b, mpfi_get_prec(result));
    mpfi_interv_fr(result, a, b);

    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
}

// void compute_i1_approximant_asympt_small(mpfr_t result, const mpfr_t x, mpfr_prec_t prec)
// {
//     if (mpfr_sgn(x) == 0)
//     {
//         mpfr_set_d(result, 0, MPFR_RNDN);
//         return;
//     }
//     mpfr_t p1, eps, p2, p3, two_over_2;
//     mpfr_inits2(prec, p1, eps, p2, p3, two_over_2, (mpfr_ptr)0);
//
//     mpfr_set_d(eps, 1e-41, MPFR_RNDN);
//
//     bessel_i1(p1, x, prec, 1500, eps);
//
//     mpfr_div(p1, p1, x, MPFR_RNDN);
//     mpfr_mul_ui(p1, p1, 2, MPFR_RNDN);
//
//     mpfr_mul_d(two_over_2, x, 0.5, MPFR_RNDN);
//
//     mpfr_mul(p2, two_over_2, two_over_2, MPFR_RNDN);
//
//     mpfr_mul_d(p2, p2, 0.5, MPFR_RNDN);
//
//     mpfr_mul(p3, two_over_2, two_over_2, MPFR_RNDN);
//     mpfr_mul(p3, p3, p3, MPFR_RNDN);
//
//     mpfr_add_d(result, p1, -1, MPFR_RNDN);
//     mpfr_sub(result, result, p2, MPFR_RNDN);
//
//     mpfr_clears(p1, eps, p2, p3,two_over_2, (mpfr_ptr)0);
// }

void compute_i1_approximant_asympt_small(mpfr_t result, const mpfr_t x, mpfr_prec_t prec)
{
    if (mpfr_sgn(x) == 0) {
        mpfr_set_d(result, 0.0, MPFR_RNDN);
        return;
    }

    mpfr_t i1x, two_i1x_over_x, y, num ,eps, denom, tmp;
    mpfr_inits2(prec, i1x, two_i1x_over_x, y, num, denom, tmp, eps, (mpfr_ptr)0);
    mpfr_set_d(eps, 1e-41, MPFR_RNDN);
    // Compute I1(x)
    bessel_i1(i1x, x, prec, 1500, eps);

    // Compute 2 * I1(x) / x
    mpfr_mul_ui(two_i1x_over_x, i1x, 2, MPFR_RNDN);
    mpfr_div(two_i1x_over_x, two_i1x_over_x, x, MPFR_RNDN);

    // y = (x/2)^2 = (x^2) / 4
    mpfr_sqr(y, x, MPFR_RNDN);
    mpfr_div_ui(y, y, 4, MPFR_RNDN);

    // Numerator = 2*I1(x)/x - 1 - 0.5 * y
    mpfr_set(num, two_i1x_over_x, MPFR_RNDN);
    mpfr_sub_ui(num, num, 1, MPFR_RNDN);            // num -= 1
    mpfr_mul_d(tmp, y, 0.5, MPFR_RNDN);              // tmp = 0.5 * y
    mpfr_sub(num, num, tmp, MPFR_RNDN);              // num -= tmp

    // Denominator = y^2
    mpfr_sqr(denom, y, MPFR_RNDN);

    // Final result: result = num / denom
    mpfr_div(result, num, denom, MPFR_RNDN);

    mpfr_clears(i1x, two_i1x_over_x, y, num, eps, denom, tmp, (mpfr_ptr)0);
}

__attribute__((visibility("default")))
int bessel_i1_approximant_small(mpfi_t result, mpfi_t x, int n)
{
    mpfr_t a, b;
    mpfr_init2(a, mpfi_get_prec(result));
    mpfr_init2(b, mpfi_get_prec(result));

    mpfi_get_left(a, x);
    compute_i1_approximant_asympt_small(a, a, mpfi_get_prec(result));

    mpfi_get_right(b, x);
    compute_i1_approximant_asympt_small(b, b, mpfi_get_prec(result));
    mpfi_interv_fr(result, a, b);

    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
}

__attribute__((visibility("default")))
int bessel_i1_approximant_big(mpfi_t result, mpfi_t x, int n)
{
    mpfr_t a, b;
    mpfr_init2(a, mpfi_get_prec(result));
    mpfr_init2(b, mpfi_get_prec(result));

    mpfr_t eps;
    mpfr_init2(eps, mpfi_get_prec(result));

    mpfr_set_d(eps, 1e-41, MPFR_RNDN);

    mpfi_get_left(a, x);
    compute_i1_approximant_asympt_big(a, a, mpfi_get_prec(result));

    mpfi_get_right(b, x);
    compute_i1_approximant_asympt_big(b, b, mpfi_get_prec(result));
    mpfi_interv_fr(result, a, b);

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(eps);
    return 0;
}
