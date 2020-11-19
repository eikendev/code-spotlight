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
 * Precomputation of μ from HAC, algorithm. 14.42. Here we have μ =
 * floor(b^(2k)/m) with radix b = 2^8 and k = 32. In Python, we can write
 * ```python3
 * def get_mu(b, k, m):
 *     return hex((b**(2 * k)) // m)
 *
 * L = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
 * print(get_mu(2**8, 32, L))
 * ```
 */
static const crypto_uint32 mu[33] =
    { 0x1B, 0x13, 0x2C, 0x0A, 0xA3, 0xE5, 0x9C, 0xED, 0xA7, 0x29, 0x63, 0x08, 0x5D, 0x21, 0x06, 0x21,
    0xEB, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x0F
};

/* Whether or not `a` is less than `b`. Inputs must be representable by 2 bytes. */
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
     * changed. If mask is all ones, the inner XOR here can be seen as
     * calculating the "difference" between subtraction and result. Using the
     * mask, we conditionally apply that "difference", i.e., we write
     * subtraction into result. */
    for (i = 0; i < 32; i++)
        result->v[i] ^= mask & (result->v[i] ^ subtraction[i]);
}

/* Reduce coefficients of x before calling barrett_reduce */
/* See HAC, algorithm 14.42. */
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

    /* Next we calculate q3. First, you might be wondering why there is no q1.
     * Instead of calculating q1 explicitly, we only consider x[31:], which is
     * equivalent to using the q1 from the book. Further note that the q2 here
     * does not strictly equal the q2 from the book. Here, q2 is used as a
     * larger buffer, where the least-significant digits are always zero. */
    for (i = 0; i < 33; i++)
        for (j = 0; j < 33; j++)
            if (i + j >= 31)
                q2[i + j] += mu[i] * x[j + 31];

    /* From step 2 onwards, we will only operate on q3. Since we have
     * q3 = q2 + 33 (reflecting q3 <- ⌊q2 / b^(k+1)⌋ from the book), we need to
     * carry any overflowed digits into q3. */
    carry = q2[31] >> 8;
    q2[32] += carry;
    carry = q2[32] >> 8;
    q2[33] += carry;

    /* Step 2 */

    crypto_uint32 r1[33];

    /* The reduction mod b^(k+1) is equivalent to disregarding x[33:63]. */
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

    /* Add digits pairwise, store the carry for each pair. */
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

/* Convert `s` into a representation similar to radix 2^3. The difference is
 * that the digits will be stored as numbers in [-4,...,3], so that the sum
 * over all r[i] * 8^i equals `s`. Since a number is represented by 255 bits,
 * and this representation uses 3 bits per digit, we need 255 / 3 = 85 digits
 * to represent a number. Note that this function runs in constant time.*/
void sc25519_window3(
    signed char r[85],
    const sc25519 *s
)
{
    int i;

    /* As a first step, we convert to the radix 2^3, but without signed digits.
     * Before we had 32 digits in radix 2^8, now there are 85 in radix 2^3.
     * Note that we can fit 3 digits from s (which are 3 * 8 = 24 bits) into 8
     * digits in r (which are 8 * 3 = 24 bits). To illustrate, the bit
     * transformation can be described like this:
     * [ 000 000 100 111 111 221 222 222 ], where `0` indicates a bit from the
     * first digits of s, `1` a bit from the second, and `2` a bit from the
     * third digit respectively. Also note that the individual digits are
     * written in a way so the most-significant bit is on the very left. Taking
     * the third digit, `100`, as an example: it contains the two most
     * significant bits of the first digit of s, and the least significant bit
     * of the second digit of s. The following loop will process the first 10
     * pairs, i.e., the first 30 bytes of s. Afterwards, there are only 15 bits
     * left. These can be fit into five digits in r. The last bit of the last
     * digit in s is not used. */

    for (i = 0; i < 10; i++) {
        r[8 * i + 0] = s->v[3 * i + 0] & 7; // `000`
        r[8 * i + 1] = (s->v[3 * i + 0] >> 3) & 7; // `000`
        r[8 * i + 2] = (s->v[3 * i + 0] >> 6) & 7; // `X00`
        r[8 * i + 2] ^= (s->v[3 * i + 1] << 2) & 7; // `100`
        r[8 * i + 3] = (s->v[3 * i + 1] >> 1) & 7; // `111`
        r[8 * i + 4] = (s->v[3 * i + 1] >> 4) & 7; // `111`
        r[8 * i + 5] = (s->v[3 * i + 1] >> 7) & 7; // `XX1`
        r[8 * i + 5] ^= (s->v[3 * i + 2] << 1) & 7; // `221`
        r[8 * i + 6] = (s->v[3 * i + 2] >> 2) & 7; // `222`
        r[8 * i + 7] = (s->v[3 * i + 2] >> 5) & 7; // `222`
    }

    r[8 * i + 0] = s->v[3 * i + 0] & 7; // `000`
    r[8 * i + 1] = (s->v[3 * i + 0] >> 3) & 7; // `000`
    r[8 * i + 2] = (s->v[3 * i + 0] >> 6) & 7; // `X00`
    r[8 * i + 2] ^= (s->v[3 * i + 1] << 2) & 7; // `100`
    r[8 * i + 3] = (s->v[3 * i + 1] >> 1) & 7; // `111`
    r[8 * i + 4] = (s->v[3 * i + 1] >> 4) & 7; // `111`

    /* Next, we need to adjust the range of the individual digits. Since we
     * cannot simply subtract the digits by some value, as that would change
     * the total value of the number, we need to take care of the carry.
     * Specifically, we check whether the digit has the most significant bit
     * set (the third bit, due to radix 2^8). If that is the case, we subtract
     * 8 from it. For instance, a digit with the value 4 will be transformed
     * into -4. Let's show why this works by example: Assume that before this
     * transformation, we only had the digits [ 0 4 0 ]. As a decimal number,
     * that is 0*8^0 + 4*8^1 + 0*8^2 = 32. After our transformation, we have
     * the digits [ 0 -4 1 ], because the onderflow was carried to the next
     * digit. The number is still the same, because 0*8^0 - 4*8^1 + 1*8^2 = 32.
     * Why does this work? Let's replace the value of the digit with the symbol
     * n, and the index 1 with i. So, before we had n*8^i = 4*8^1. After the
     * transformation, we are left with (n-8)*8^i + 8^(i+1) = -4*8^1 + 1*8^2.
     * It holds that n*8^i = (n-8)*8^i + 8^(i+1), because 0 = (-8)*8^i +
     * 8^(i+1). */

    char carry = 0;

    for (i = 0; i < 84; i++) {
        r[i] += carry;
        r[i + 1] += r[i] >> 3;
        r[i] &= 7;
        carry = r[i] >> 2;
        r[i] -= carry << 3;
    }

    r[84] += carry;
}

/* TODO: Add explanation for this function. */
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
