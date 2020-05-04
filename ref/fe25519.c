#include "fe25519.h"

/* Whether or not `a` is equal to `b`. Inputs must be representable by 2 bytes. */
static crypto_uint32 equal(
    crypto_uint32 a,
    crypto_uint32 b
)
{
    /* If a = b, then x = 0. */
    crypto_uint32 x = a ^ b;

    /* If a = b, then x = 4294967295. Otherwise, x = [0,65534]. */
    x -= 1;

    /* If a = b, then x = 1. Otherwise, x = 0. */
    x >>= 31;

    return x;
}

/* Whether or not `a` is greater or equal to `b`. Inputs must be representable
 * by 2 bytes. */
static crypto_uint32 ge(
    crypto_uint32 a,
    crypto_uint32 b
)
{
    unsigned int x = a;

    /* If a >= b, then x = [0,65535]. Otherwise, x = [4294901761,4294967295]. */
    x -= (unsigned int) b;

    /* If a >= b, then x = 0. Otherwise, x = 1. */
    x >>= 31;

    /* If a >= b, then x = 1. Otherwise, x = 0. */
    x ^= 1;

    return x;
}

/* Calculates 19 * a. */
static crypto_uint32 times19(
    crypto_uint32 a
)
{
    return (a << 4) + (a << 1) + a;
}

/* Calculates 38 * a. */
static crypto_uint32 times38(
    crypto_uint32 a
)
{
    return (a << 5) + (a << 2) + (a << 1);
}

/* Reduce `r` modulo 2^255 after an addition or subtraction. */
/* TODO: Why do we iterate four times here? */
static void reduce_add_sub(
    fe25519 *r
)
{
    crypto_uint32 temporary;
    int i, rep;

    for (rep = 0; rep < 4; rep++) {
        /* When representing a number modulo 2^255 in 32 buckets, the last
         * bucket is capped to seven bits. Everything that does not fit into
         * those seven bits of the last byte can be seen as belonging to a 33rd
         * bucket, which of course does not exist. The content of the 33rd
         * bucket would represent the multiples of 2^255. Now, since 2^255 is
         * congruent to 19 modulo 2^255-19, we can also represent it as a
         * multiple of 19. */
        temporary = r->v[31] >> 7;
        r->v[31] &= 127;
        temporary = times19(temporary);
        r->v[0] += temporary;

        /* After the above transformation, the first bucket could be could have
         * overflowed. We iterate over all buckets, and shift the carry along
         * the whole number. */
        for (i = 0; i < 31; i++) {
            temporary = r->v[i] >> 8;
            r->v[i + 1] += temporary;
            r->v[i] &= 255;
        }
    }
}

/* Reduce `r` modulo 2^255 after a multiplication. See reduce_add_sub() for
 * details on this function. The only difference is the number of iterations. */
/* TODO: Why do we iterate twice here? */
static void reduce_mul(
    fe25519 *r
)
{
    crypto_uint32 temporary;
    int i, rep;

    for (rep = 0; rep < 2; rep++) {
        temporary = r->v[31] >> 7;
        r->v[31] &= 127;
        temporary = times19(temporary);
        r->v[0] += temporary;

        for (i = 0; i < 31; i++) {
            temporary = r->v[i] >> 8;
            r->v[i + 1] += temporary;
            r->v[i] &= 255;
        }
    }
}

/* Reduce `r` modulo 2^255-19. This assumes that `r` has been reduced modulo
 * 2^255 and that the digits are already reduced. */
void fe25519_freeze(
    fe25519 *r
)
{
    int i;

    /* First we check if r is in [2^255-19,2^255-1]. Note that with reduced
     * digits, we can hence check if all bits are set for all digits but the
     * first. The first digit must be greater or equal to 2^8-19 = 237.
     * Remember that we reduce modulo 2^255-19. The last byte of 2^255-19 in
     * our representation is 2^8-19. */
    crypto_uint32 m = equal(r->v[31], 127);
    for (i = 30; i > 0; i--)
        m &= equal(r->v[i], 255);
    m &= ge(r->v[0], 237);

    /* If r is in [2^255-19,2^255-1], then m = 1. Otherwise, m = 0. From that
     * we build a mask. If m = 1, then -m is a mask with all ones (due to the
     * representation as Two's complement). If m = 0, m is a mask with all
     * zeros. */
    m = -m;

    /* If the number is in [2^255-19,2^255-1], we now proceed by setting all
     * the digits but the first to zero. To make the last byte "wrap around"
     * when reaching the modulus, we subtract the modulus 2^8-19 = 237 from the
     * first digit. */
    r->v[31] -= m & 127;
    for (i = 30; i > 0; i--)
        r->v[i] -= m & 255;
    r->v[0] -= m & 237;
}

/* Convert `x` to our bucket representation. Note that the most significant bit
 * of the last byte is disregarded, i.e., the function only accepts numbers in
 * [0,2^255-1]. */
void fe25519_unpack(
    fe25519 *r,
    const unsigned char x[32]
)
{
    int i;

    for (i = 0; i < 32; i++)
        r->v[i] = x[i];

    r->v[31] &= 127;
}

/* Converts x to the packed format. This assumes input x has been reduced
 * modulo 2^255. */
void fe25519_pack(
    unsigned char r[32],
    const fe25519 *x
)
{
    int i;

    fe25519 y = *x;

    /* Reduce the number modulo p. */
    fe25519_freeze(&y);

    for (i = 0; i < 32; i++)
        r[i] = y.v[i];
}

/* Returns 1 if x = y, otherwise returns 0. Note that this function does not
 * run in constant time. */
int fe25519_iseq_vartime(
    const fe25519 *x,
    const fe25519 *y
)
{
    int i;

    fe25519 temporary1 = *x;
    fe25519 t2 = *y;
    fe25519_freeze(&temporary1);
    fe25519_freeze(&t2);

    for (i = 0; i < 32; i++)
        if (temporary1.v[i] != t2.v[i])
            return 0;

    return 1;
}

/* If `b` is `1`, moves `x` into `r`. For `b` equal to `0`, `r` is left
 * unchanged. Hence the name "conditional move". This function runs in constant
 * time. */
void fe25519_cmov(
    fe25519 *r,
    const fe25519 *x,
    unsigned char b
)
{
    int i;

    /* The char is written into an integer so we can invert it in the next step. */
    crypto_uint32 mask = b;

    /* Next, we build a mask. For b = 1, the resulting mask is all ones (due to
     * the representation as Two's complement). Otherwise, the mask is all
     * zeros. */
    mask = -mask;

    /* The XOR here can be seen as the "difference" between x and r. Using the
     * mask, we conditionally apply that "difference". */
    for (i = 0; i < 32; i++)
        r->v[i] ^= mask & (x->v[i] ^ r->v[i]);
}

/* Sets `r` to one. */
void fe25519_setone(
    fe25519 *r
)
{
    int i;

    r->v[0] = 1;

    for (i = 1; i < 32; i++)
        r->v[i] = 0;
}

/* Sets `r` to zero. */
void fe25519_setzero(
    fe25519 *r
)
{
    int i;

    for (i = 0; i < 32; i++)
        r->v[i] = 0;
}

/* Uses fe25519_sub() to negate `x` by subtracting it from zero. */
void fe25519_neg(
    fe25519 *r,
    const fe25519 *x
)
{
    fe25519 temporary;
    int i;

    for (i = 0; i < 32; i++)
        temporary.v[i] = x->v[i];

    fe25519_setzero(r);
    fe25519_sub(r, r, &temporary);
}

/* Returns `1` if `x` is odd, `0` otherwise. */
unsigned char fe25519_getparity(
    const fe25519 *x
)
{
    fe25519 temporary = *x;
    fe25519_freeze(&temporary);
    return temporary.v[0] & 1;
}

/* Adds `x` and `y`, and writes the result into `r`. */
void fe25519_add(
    fe25519 *r,
    const fe25519 *x,
    const fe25519 *y
)
{
    int i;

    for (i = 0; i < 32; i++)
        r->v[i] = x->v[i] + y->v[i];

    /* The digits could have overflowed, so we need to reduce them. */
    reduce_add_sub(r);
}

/* Subtracts `y` from `x`, and writes the result into `r`. */
void fe25519_sub(
    fe25519 *r,
    const fe25519 *x,
    const fe25519 *y
)
{
    int i;
    crypto_uint32 temporary[32];

    /* To prevent the individual buckets from underflowing, we need to make
     * sure that the buckets are "full enough" so we can subtract both numbers
     * (x and y) without causing an underflow. For the first bucket, we could
     * have x = y = 2^8-19 in the "worst" case. This means, we need to
     * initialize the bucket with 2*(2^8-19). When adding x and subtracting y
     * in this bucket, the bucket is guaranteed not to overflow. For the last
     * bucket, we assume the operands not to have the last bit set, which is
     * why we initialize it with 2*(2^7-1). All the buckets inbetween deal with
     * full-byte numbers, so they need to be filled with 2*(2^8-1) before the
     * subtraction. But why can we just "fill" the buckets as a preparation,
     * doesn't this yield a wrong result for the subtraction? Since we operate
     * in a finite field, we can add p arbitrarily often. In this case, we
     * would actually calculate 2*p-x-y. If you write down the prepared
     * buckets, this would look like
     * [0x00fe 0x01fe 0x01fe ... 0x01fe 0x01fe 0x01da] (big-endian ordering).
     * We also know that 2*p = 2*(2^255-19) = [0xff 0xff ... 0xff 0xff 0xda].
     * One can easily see that when overlapping the buckets, that they sum up
     * to 2*p. */

    temporary[0] = x->v[0] + 0x1da; /* 0x1da = 2*(2^8-19) */
    temporary[31] = x->v[31] + 0xfe; /* 0x1fe = 2*(2^8-1) */

    for (i = 1; i < 31; i++)
        temporary[i] = x->v[i] + 0x1fe; /* 0xfe = 2*(2^7-1) */

    for (i = 0; i < 32; i++)
        r->v[i] = temporary[i] - y->v[i];

    /* The digits could have overflowed, so we need to reduce them. */
    reduce_add_sub(r);
}

/* Multiplies `y` and `x`, and writes the result into `r`. */
void fe25519_mul(
    fe25519 *r,
    const fe25519 *x,
    const fe25519 *y
)
{
    int i, j;
    crypto_uint32 temporary[63];

    /* Initialize `temporary` to zero. If we were calculating with full 32-byte
     * numbers, the multiplication could not be bound by 63 byte. But since our
     * elements are in a field of order 2^255-19, the maximum value in the
     * following multiplication is (2^255 - 18)^2. */
    for (i = 0; i < 63; i++)
        temporary[i] = 0;

    /* Multiply the two operands using "long multiplication". */
    for (i = 0; i < 32; i++)
        for (j = 0; j < 32; j++)
            temporary[i + j] += x->v[i] * y->v[j];

    /* Since the buffer is too large for us to process, we need to reduce it
     * into 32 bytes again. The second half of the buffer is used for the parts
     * exceeding 2^256-1. This is because in the first 32 bytes, we can
     * represent all numbers in [0,2^256-1]. Now, let's illustrate how we can
     * reduce the bytes of the second half. The 33rd byte can be regarded as
     * "shifted" by 32 bytes, i.e., 256 bits, to the right. Luckily, 2^256 is
     * congruent to 38, so we can multiply the 33rd byte with 38, and add it to
     * the first byte. Note that a multiplication with 38 is equal to one with
     * 2^256 in our arithmetic. */
    for (i = 32; i < 63; i++)
        r->v[i - 32] = temporary[i - 32] + times38(temporary[i]);
    r->v[31] = temporary[31];           /* result now in r[0]...r[31] */

    /* The digits could have overflowed, so we need to reduce them. */
    reduce_mul(r);
}

/* Squares `x`, and writes the result into `r`. */
void fe25519_square(
    fe25519 *r,
    const fe25519 *x
)
{
    fe25519_mul(r, x, x);
}

/* Inverts `x`, and writes the result into `r`. Since in our arithmetic we have
 * x * x^(-1) = 1 = x^(p-1), it must hold that x^(-1) = x^(p-2). Thus, we can
 * calculate x^(p-2) to derive the inverse. In our case, p is 2^255-19, so we
 * want to have x^(2^255-21). To calculate this number, we use the
 * square-and-multiply method. */
void fe25519_invert(
    fe25519 *r,
    const fe25519 *x
)
{
    int i;

    fe25519 z2;
    fe25519 z9;
    fe25519 z11;

    fe25519 z2_5_0;
    fe25519 z2_10_0;
    fe25519 z2_20_0;
    fe25519 z2_50_0;
    fe25519 z2_100_0;

    fe25519 temporary0;
    fe25519 temporary1;

    fe25519_square(&z2, x); /* 2 */

    fe25519_square(&temporary1, &z2); /* 4 */
    fe25519_square(&temporary0, &temporary1); /* 8 */
    fe25519_mul(&z9, &temporary0, x); /* 9 */

    fe25519_mul(&z11, &z9, &z2); /* 11 */

    fe25519_square(&temporary0, &z11); /* 22 */
    fe25519_mul(&z2_5_0, &temporary0, &z9); /* 31 = 2^5 - 2^0 */

    /* We now have x^31. Since the exponents will grow quickly, we'll write
     * x^(2^5 - 2^0) from now instead. When squaring this number, we have
     * x^(2*(2^5 - 2^0)) = x^(2^6 - 2^1). */

    fe25519_square(&temporary0, &z2_5_0); /* 2^6 - 2^1 */
    fe25519_square(&temporary1, &temporary0); /* 2^7 - 2^2 */
    fe25519_square(&temporary0, &temporary1); /* 2^8 - 2^3 */
    fe25519_square(&temporary1, &temporary0); /* 2^9 - 2^4 */
    fe25519_square(&temporary0, &temporary1); /* 2^10 - 2^5 */
    fe25519_mul(&z2_10_0, &temporary0, &z2_5_0); /* 2^10 - 2^0 */

    fe25519_square(&temporary0, &z2_10_0); /* 2^11 - 2^1 */
    fe25519_square(&temporary1, &temporary0); /* 2^12 - 2^2 */
    for (i = 2; i < 10; i += 2) {
        fe25519_square(&temporary0, &temporary1);
        fe25519_square(&temporary1, &temporary0);
    } /* 2^20 - 2^10 */
    fe25519_mul(&z2_20_0, &temporary1, &z2_10_0); /* 2^20 - 2^0 */

    fe25519_square(&temporary0, &z2_20_0); /* 2^21 - 2^1 */
    fe25519_square(&temporary1, &temporary0); /* 2^22 - 2^2 */
    for (i = 2; i < 20; i += 2) {
        fe25519_square(&temporary0, &temporary1);
        fe25519_square(&temporary1, &temporary0);
    } /* 2^40 - 2^20 */
    fe25519_mul(&temporary0, &temporary1, &z2_20_0); /* 2^40 - 2^0 */
    fe25519_square(&temporary1, &temporary0); /* 2^41 - 2^1 */
    fe25519_square(&temporary0, &temporary1); /* 2^42 - 2^2 */
    for (i = 2; i < 10; i += 2) {
        fe25519_square(&temporary1, &temporary0);
        fe25519_square(&temporary0, &temporary1);
    } /* 2^50 - 2^10 */
    fe25519_mul(&z2_50_0, &temporary0, &z2_10_0); /* 2^50 - 2^0 */

    fe25519_square(&temporary0, &z2_50_0); /* 2^51 - 2^1 */
    fe25519_square(&temporary1, &temporary0); /* 2^52 - 2^2 */
    for (i = 2; i < 50; i += 2) {
        fe25519_square(&temporary0, &temporary1);
        fe25519_square(&temporary1, &temporary0);
    } /* 2^100 - 2^50 */
    fe25519_mul(&z2_100_0, &temporary1, &z2_50_0); /* 2^100 - 2^0 */

    fe25519_square(&temporary1, &z2_100_0); /* 2^101 - 2^1 */
    fe25519_square(&temporary0, &temporary1); /* 2^102 - 2^2 */
    for (i = 2; i < 100; i += 2) {
        fe25519_square(&temporary1, &temporary0);
        fe25519_square(&temporary0, &temporary1);
    } /* 2^200 - 2^100 */
    fe25519_mul(&temporary1, &temporary0, &z2_100_0); /* 2^200 - 2^0 */
    fe25519_square(&temporary0, &temporary1); /* 2^201 - 2^1 */
    fe25519_square(&temporary1, &temporary0); /* 2^202 - 2^2 */
    for (i = 2; i < 50; i += 2) {
        fe25519_square(&temporary0, &temporary1);
        fe25519_square(&temporary1, &temporary0);
    } /* 2^250 - 2^50 */
    fe25519_mul(&temporary0, &temporary1, &z2_50_0); /* 2^250 - 2^0 */
    fe25519_square(&temporary1, &temporary0); /* 2^251 - 2^1 */
    fe25519_square(&temporary0, &temporary1); /* 2^252 - 2^2 */
    fe25519_square(&temporary1, &temporary0); /* 2^253 - 2^3 */
    fe25519_square(&temporary0, &temporary1); /* 2^254 - 2^4 */
    fe25519_square(&temporary1, &temporary0); /* 2^255 - 2^5 */
    fe25519_mul(r, &temporary1, &z11); /* 2^255 - 21 */
}

/* Calculate x^2523 and store the result in `r`. This uses the
 * square-and-multiply approach. See fe25519_invert() for a bit of guidance on
 * the calculation. */
void fe25519_pow2523( fe25519 *r, const fe25519 *x) { int i;

    fe25519 z2;
    fe25519 z9;
    fe25519 z11;

    fe25519 z2_5_0;
    fe25519 z2_10_0;
    fe25519 z2_20_0;
    fe25519 z2_50_0;
    fe25519 z2_100_0;

    fe25519 temporary;

    fe25519_square(&z2, x); /* 2 */

    fe25519_square(&temporary, &z2); /* 4 */
    fe25519_square(&temporary, &temporary); /* 8 */
    fe25519_mul(&z9, &temporary, x); /* 9 */

    fe25519_mul(&z11, &z9, &z2); /* 11 */

    fe25519_square(&temporary, &z11); /* 22 */
    fe25519_mul(&z2_5_0, &temporary, &z9); /* 31 = 2^5 - 2^0 */

    fe25519_square(&temporary, &z2_5_0); /* 2^6 - 2^1 */
    for (i = 1; i < 5; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^10 - 2^5 */
    fe25519_mul(&z2_10_0, &temporary, &z2_5_0); /* 2^10 - 2^0 */

    fe25519_square(&temporary, &z2_10_0); /* 2^11 - 2^1 */
    for (i = 1; i < 10; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^20 - 2^10 */
    fe25519_mul(&z2_20_0, &temporary, &z2_10_0); /* 2^20 - 2^0 */

    fe25519_square(&temporary, &z2_20_0); /* 2^21 - 2^1 */
    for (i = 1; i < 20; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^40 - 2^20 */
    fe25519_mul(&temporary, &temporary, &z2_20_0); /* 2^40 - 2^0 */
    fe25519_square(&temporary, &temporary); /* 2^41 - 2^1 */
    for (i = 1; i < 10; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^50 - 2^10 */
    fe25519_mul(&z2_50_0, &temporary, &z2_10_0); /* 2^50 - 2^0 */

    fe25519_square(&temporary, &z2_50_0); /* 2^51 - 2^1 */
    for (i = 1; i < 50; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^100 - 2^50 */
    fe25519_mul(&z2_100_0, &temporary, &z2_50_0); /* 2^100 - 2^0 */

    fe25519_square(&temporary, &z2_100_0); /* 2^101 - 2^1 */
    for (i = 1; i < 100; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^200 - 2^100 */
    fe25519_mul(&temporary, &temporary, &z2_100_0); /* 2^200 - 2^0 */
    fe25519_square(&temporary, &temporary); /* 2^201 - 2^1 */
    for (i = 1; i < 50; i++) {
        fe25519_square(&temporary, &temporary);
    } /* 2^250 - 2^50 */
    fe25519_mul(&temporary, &temporary, &z2_50_0); /* 2^250 - 2^0 */
    fe25519_square(&temporary, &temporary); /* 2^251 - 2^1 */
    fe25519_square(&temporary, &temporary); /* 2^252 - 2^2 */
    fe25519_mul(r, &temporary, x); /* 2^252 - 3 */
}
