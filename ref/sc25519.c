#include "sc25519.h"

/*
 * Arithmetic modulo the group order L.
 * L = 2^252 + 27742317777372353535851937790883648493
 */
static const crypto_uint32 L[32] =
    { 0xED, 0xD3, 0xF5, 0x5C, 0x1A, 0x63, 0x12, 0x58, 0xD6, 0x9C, 0xF7, 0xA2, 0xDE, 0xF9, 0xDE, 0x14,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
};

/*
 * Precomputation of μ from HAC, Alg. 14.42.
 * μ = floor(b^(2k)/m)
 * Radix b = 2^64, k = 4.
 */
static const crypto_uint32 mu[33] =
    { 0x1B, 0x13, 0x2C, 0x0A, 0xA3, 0xE5, 0x9C, 0xED, 0xA7, 0x29, 0x63, 0x08, 0x5D, 0x21, 0x06, 0x21,
    0xEB, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x0F
};

/* Whether or not `a` is less than `b`. Input must be representable in 16-bit. */
static crypto_uint32 lt(
    crypto_uint32 a,
    crypto_uint32 b
)
{
    unsigned int x = a;

    /* Subtract b from a. If the result is in [0,65535], then a is bigger than
     * b. Otherwise, the result is in [4294901761,4294967295]. */
    x -= (unsigned int) b;

    /* Shift the result. Now if x=0, then a is bigger than b, so return false. */
    x >>= 31;

    return x;
}

/* Add or subtract L once. The result is an element in a group with order L,
 * given that the input is already "close" to the reduced element. Assumes that
 * the digits of result are already reduced. */
static void reduce_add_sub(
    sc25519 *result
)
{
    /* TODO: Rename these variables to something more useful. */
    crypto_uint32 pb = 0;
    crypto_uint32 b;
    int i;
    unsigned char subtraction[32];

    /* Calculate result - L. See the calculation of r1 - r2 in barrett_reduce()
     * for details. In essence, the loop handles the case where
     * result[i] < L[i]. */
    for (i = 0; i < 32; i++) {
        pb += L[i];
        b = lt(result->v[i], pb);
        subtraction[i] = result->v[i] - pb + (b << 8);
        pb = b;
    }

    /* If for the last pair it holds that result[i] < L[i], then b will be 1, so this mask will
     * be all zeros. Otherwise, the mask will be all ones. */
    crypto_uint32 mask = b - 1;

    /* If we have result < L, then the subtraction underflowed. Note that we
     * luckily did not change the result buffer yet! So if the subtraction
     * underflowed, we can actually leave `result` as is, given our assumption
     * that `result` is already "close" to the reduced element. Note that in
     * this case, the mask is set to all zeros. Hence, `result` will not be
     * changed. */
    /* TODO: What exactly does the loop do if mask is all ones? */
    for (i = 0; i < 32; i++)
        result->v[i] ^= mask & (result->v[i] ^ subtraction[i]);
}

/* Reduce coefficients of x before calling barrett_reduce */
/* See HAC, Alg. 14.42. */
static void barrett_reduce(
    sc25519 *r,
    const crypto_uint32 x[64]
)
{
    int i, j;

    crypto_uint32 q2[66];
    crypto_uint32 r2[33];

    /* Initialize q2 to zero. */
    for (i = 0; i < 66; ++i)
        q2[i] = 0;

    /* Initialize r2 to zero. */
    for (i = 0; i < 33; ++i)
        r2[i] = 0;

    /* Step 1 */

    crypto_uint32 *q3 = q2 + 33;
    crypto_uint32 carry;

    /* Calculate q3. Note that we do not explicitly calculate q2, i.e., the
     * result in q2 is not equivalent to Step 1 from the book. Instead, we
     * directly divide by b^(k+1) and take the floor of the result. */
    for (i = 0; i < 33; i++)
        for (j = 0; j < 33; j++)
            if (i + j >= 31)
                q2[i + j] += mu[i] * x[j + 31];

    /* TODO: Why is this needed, considering that we want to round down? */
    carry = q2[31] >> 8;
    q2[32] += carry;
    carry = q2[32] >> 8;
    q2[33] += carry;

    /* Step 2 */

    crypto_uint32 r1[33];

    /* TODO: Where did the reduction go? */
    for (i = 0; i < 33; i++)
        r1[i] = x[i];

    /* Calculate r2. This is similar to the step in sc25519_mul(). The
     * difference is that the second operand is bigger in size. */
    for (i = 0; i < 32; i++)
        for (j = 0; j < 33; j++)
            if (i + j < 33)
                r2[i + j] += L[i] * q3[j];

    /* Just as in sc25519_mul(), the digits of the result must be reduced. */
    for (i = 0; i < 32; i++) {
        carry = r2[i] >> 8;
        r2[i + 1] += carry;
        r2[i] &= 0xff;
    }

    crypto_uint32 b;
    crypto_uint32 pb = 0;

    /* Calculate r1 - r2. Note that we could have r1[i] < r2[i]. In that case,
     * we need to add the radix to our digit, and subtract the next digit in
     * round i+1 by one more. Look up the Austrian subtraction method for
     * details. Note that for the last pair, we do not handle a carry at all.
     * See the next step for details. */
    for (i = 0; i < 32; i++) {
        pb += r2[i];
        b = lt(r1[i], pb);
        r->v[i] = r1[i] - pb + (b << 8);
        pb = b;
    }

    /* Step 3 */

    /* This step was left out by the reference. I'm not sure myself if it can
     * be that r1 < r2. If you express r2 in terms of x, you will have roughly
     * r2 ~~ x/(b^2) mod b^(k+1). */

    /* Step 4 */

    /* Read HAC, Note 14.44 (ii). The loop in Step 4 is repeated at most twice.
     * Hence, we always call it twice and save checking the loop condition. */

    /* Since r1 and r2 were already reduced before subtraction, the digits of r
     * are also already reduced. */

    reduce_add_sub(r);
    reduce_add_sub(r);
}

/* Converts a 32-byte integer into a scalar of a group with order L. Note that
 * since we store bytes, the digits are already reduced to the chosen radix. */
void sc25519_from32bytes(
    sc25519 *result,
    const unsigned char x[32]
)
{
    int i;
    crypto_uint32 t[64];

    for (i = 0; i < 32; i++)
        t[i] = x[i];

    for (i = 32; i < 64; ++i)
        t[i] = 0;

    barrett_reduce(result, t);
}

/* Converts a 64-byte integer into a scalar of a group with order L. Note that
 * since we store bytes, the digits are already reduced to the chosen radix. */
void sc25519_from64bytes(
    sc25519 *result,
    const unsigned char x[64]
)
{
    int i;
    crypto_uint32 t[64];

    for (i = 0; i < 64; i++)
        t[i] = x[i];

    barrett_reduce(result, t);
}

/* Converts a scalar of a group with order L into a 64-byte integer. */
void sc25519_to32bytes(
    unsigned char result[32],
    const sc25519 *x
)
{
    int i;

    /* We only need to copy the digits. */
    for (i = 0; i < 32; i++)
        result[i] = x->v[i];
}

/* Adds two elements of a group with order L. */
void sc25519_add(
    sc25519 *result,
    const sc25519 *x,
    const sc25519 *y
)
{
    int i, carry;

    /* Multiply digits pairwise, store the carry for each pair. */
    for (i = 0; i < 32; i++)
        result->v[i] = x->v[i] + y->v[i];

    /* Reduce the digits of the result: Apply the carry on the next digit. */
    for (i = 0; i < 31; i++) {
        carry = result->v[i] >> 8;
        result->v[i + 1] += carry;
        result->v[i] &= 0xff;
    }

    /* Add or subtract L once. We don't need a full reduction here, since the
     * result can be at most one L off. */
    reduce_add_sub(result);
}

/* Multiplies two elements of a group with order L. */
void sc25519_mul(
    sc25519 *result,
    const sc25519 *x,
    const sc25519 *y
)
{
    int i, j, carry;
    crypto_uint32 result_with_carry[64];

    /* Initialize result_with_carry to zero. */
    for (i = 0; i < 64; i++)
        result_with_carry[i] = 0;

    /* Multiply the two operands using "long multiplication". */
    for (i = 0; i < 32; i++)
        for (j = 0; j < 32; j++)
            result_with_carry[i + j] += x->v[i] * y->v[j];

    /* Reduce the digits of the result: Apply the carry on the next digit. */
    for (i = 0; i < 63; i++) {
        carry = result_with_carry[i] >> 8;
        result_with_carry[i + 1] += carry;
        result_with_carry[i] &= 0xff;
    }

    /* Here we need to do a fully reduction. This is because we don't have
     * useful bounds for the result. */
    barrett_reduce(result, result_with_carry);
}

/* TODO: Understand and document sc25519_window3(). */
/* Convert s into a representation of the form \sum_{i=0}^{84}r[i]2^3 with r[i]
 * in {-4,...,3} */
void sc25519_window3(
    signed char r[85],
    const sc25519 *s
)
{
    char carry;
    int i;

    for (i = 0; i < 10; i++) {
        r[8 * i + 0] = s->v[3 * i + 0] & 7;
        r[8 * i + 1] = (s->v[3 * i + 0] >> 3) & 7;
        r[8 * i + 2] = (s->v[3 * i + 0] >> 6) & 7;
        r[8 * i + 2] ^= (s->v[3 * i + 1] << 2) & 7;
        r[8 * i + 3] = (s->v[3 * i + 1] >> 1) & 7;
        r[8 * i + 4] = (s->v[3 * i + 1] >> 4) & 7;
        r[8 * i + 5] = (s->v[3 * i + 1] >> 7) & 7;
        r[8 * i + 5] ^= (s->v[3 * i + 2] << 1) & 7;
        r[8 * i + 6] = (s->v[3 * i + 2] >> 2) & 7;
        r[8 * i + 7] = (s->v[3 * i + 2] >> 5) & 7;
    }

    r[8 * i + 0] = s->v[3 * i + 0] & 7;
    r[8 * i + 1] = (s->v[3 * i + 0] >> 3) & 7;
    r[8 * i + 2] = (s->v[3 * i + 0] >> 6) & 7;
    r[8 * i + 2] ^= (s->v[3 * i + 1] << 2) & 7;
    r[8 * i + 3] = (s->v[3 * i + 1] >> 1) & 7;
    r[8 * i + 4] = (s->v[3 * i + 1] >> 4) & 7;

    /* Making it signed */
    carry = 0;
    for (i = 0; i < 84; i++) {
        r[i] += carry;
        r[i + 1] += r[i] >> 3;
        r[i] &= 7;
        carry = r[i] >> 2;
        r[i] -= carry << 3;
    }

    r[84] += carry;
}

/* TODO: Understand and document sc25519_2interleave2(). */
void sc25519_2interleave2(
    unsigned char r[127],
    const sc25519 *s1,
    const sc25519 *s2
)
{
    int i;

    for (i = 0; i < 31; i++) {
        r[4 * i] = (s1->v[i] & 3) ^ ((s2->v[i] & 3) << 2);
        r[4 * i + 1] = ((s1->v[i] >> 2) & 3) ^ (((s2->v[i] >> 2) & 3) << 2);
        r[4 * i + 2] = ((s1->v[i] >> 4) & 3) ^ (((s2->v[i] >> 4) & 3) << 2);
        r[4 * i + 3] = ((s1->v[i] >> 6) & 3) ^ (((s2->v[i] >> 6) & 3) << 2);
    }

    r[124] = (s1->v[31] & 3) ^ ((s2->v[31] & 3) << 2);
    r[125] = ((s1->v[31] >> 2) & 3) ^ (((s2->v[31] >> 2) & 3) << 2);
    r[126] = ((s1->v[31] >> 4) & 3) ^ (((s2->v[31] >> 4) & 3) << 2);
}
