#include "fe25519.h"
#include "sc25519.h"
#include "ge25519.h"

/* Arithmetic on the twisted Edwards curve -x^2 + y^2 = 1 + dx^2y^2. */

/* d = -(121665/121666)
 *   = 37095705934669439343138083508754565189542113879843219016388785533085940283555 */
static const fe25519 ge25519_ecd = {{
    0xA3, 0x78, 0x59, 0x13, 0xCA, 0x4D, 0xEB, 0x75,
    0xAB, 0xD8, 0x41, 0x41, 0x4D, 0x0A, 0x70, 0x00,
    0x98, 0xE8, 0x79, 0x77, 0x79, 0x40, 0xC7, 0x8C,
    0x73, 0xFE, 0x6F, 0x2B, 0xEE, 0x6C, 0x03, 0x52
}};

/* k = 2*d */
static const fe25519 ge25519_ec2d = {{
    0x59, 0xF1, 0xB2, 0x26, 0x94, 0x9B, 0xD6, 0xEB,
    0x56, 0xB1, 0x83, 0x82, 0x9A, 0x14, 0xE0, 0x00,
    0x30, 0xD1, 0xF3, 0xEE, 0xF2, 0x80, 0x8E, 0x19,
    0xE7, 0xFC, 0xDF, 0x56, 0xDC, 0xD9, 0x06, 0x24
}};

/* sqrt(-1) */
static const fe25519 ge25519_sqrtm1 = {{
    0xB0, 0xA0, 0x0E, 0x4A, 0x27, 0x1B, 0xEE, 0xC4,
    0x78, 0xE4, 0x2F, 0xAD, 0x06, 0x18, 0x43, 0x2F,
    0xA7, 0xD7, 0xFB, 0x3D, 0x99, 0x00, 0x4D, 0x2B,
    0x0B, 0xDF, 0xC1, 0x4F, 0x80, 0x24, 0x83, 0x2B
}};

#define ge25519_p3 ge25519

/* This is sometimes called the completed point representation. See [0].
 * [0] https://doc-internal.dalek.rs/curve25519_dalek/backend/serial/curve_models/index.html */
typedef struct {
    fe25519 x;
    fe25519 z;
    fe25519 y;
    fe25519 t;
} ge25519_p1p1;

/* The projective representation. It must hold that x = X/Z and y = Y/Z. See
 * [0].
 * [0] https://www.hyperelliptic.org/EFD/g1p/auto-twisted-projective.html */
typedef struct {
    fe25519 x;
    fe25519 y;
    fe25519 z;
} ge25519_p2;

/* The affine representation. */
typedef struct {
    fe25519 x;
    fe25519 y;
} ge25519_aff;

/* Packed coordinates of the base point B = (B_x, B_y).
 * B_x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
 * B_y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
 * Note that this is in extended homogeneous coordinates (X, Y, Z, T), with
 * x = X/Z, y = Y/Z, x * y = T/Z, and Z = 1. */
const ge25519 ge25519_base = {
    {{
        0x1A, 0xD5, 0x25, 0x8F, 0x60, 0x2D, 0x56, 0xC9,
        0xB2, 0xA7, 0x25, 0x95, 0x60, 0xC7, 0x2C, 0x69,
        0x5C, 0xDC, 0xD6, 0xFD, 0x31, 0xE2, 0xA4, 0xC0,
        0xFE, 0x53, 0x6E, 0xCD, 0xD3, 0x36, 0x69, 0x21
    }},
    {{
        0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
    }},
    {{
        0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    }},
    {{
        0xA3, 0xDD, 0xB7, 0xA5, 0xB3, 0x8A, 0xDE, 0x6D,
        0xF5, 0x52, 0x51, 0x77, 0x80, 0x9F, 0xF0, 0x20,
        0x7D, 0xE3, 0xAB, 0x64, 0x8E, 0x4E, 0xEA, 0x66,
        0x65, 0x76, 0x8B, 0xD7, 0x0F, 0x5F, 0x87, 0x67
    }}
};

/* Multiples of the base point in affine representation. For each of the 85
 * indices, there are 5 possible multiples (disregarding the signs), and
 * 5 * 85 = 425. */
static const ge25519_aff ge25519_base_multiples_affine[425] = {
#include "ge25519_base.data"
};

/* TODO: Add explanation for this function. */
static void p1p1_to_p2(
    ge25519_p2 *r,
    const ge25519_p1p1 *p
)
{
    fe25519_mul(&r->x, &p->x, &p->t);
    fe25519_mul(&r->y, &p->y, &p->z);
    fe25519_mul(&r->z, &p->z, &p->t);
}

/* TODO: Add explanation for this function. */
static void p1p1_to_p3(
    ge25519_p3 *r,
    const ge25519_p1p1 *p
)
{
    p1p1_to_p2((ge25519_p2 *) r, p);
    fe25519_mul(&r->t, &p->x, &p->y);
}

/* Add `q` and `r`, and store the result in `r`. See [0].
 * [0] https://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-madd-2008-hwcd-3 */
static void ge25519_mixadd2(
    ge25519_p3 *r,
    const ge25519_aff *q
)
{
    fe25519 a, b, t1, t2, c, d, e, f, g, h, qt;

    /* The addition assumes that Z2=1. We can simply calculate T2 = X2*Y2. For
     * details, have a look at the representation using extended coordinates. */
    fe25519_mul(&qt, &q->x, &q->y);

    /* A = (Y1-X1)*(Y2-X2) */
    /* B = (Y1+X1)*(Y2+X2) */
    fe25519_sub(&a, &r->y, &r->x);
    fe25519_add(&b, &r->y, &r->x);
    fe25519_sub(&t1, &q->y, &q->x);
    fe25519_add(&t2, &q->y, &q->x);
    fe25519_mul(&a, &a, &t1);
    fe25519_mul(&b, &b, &t2);

    fe25519_sub(&e, &b, &a); /* E = B-A */
    fe25519_add(&h, &b, &a); /* H = B+A */

    /* C = T1*k*T2 */
    fe25519_mul(&c, &r->t, &qt);
    fe25519_mul(&c, &c, &ge25519_ec2d);

    fe25519_add(&d, &r->z, &r->z); /* D = Z1*2 */
    fe25519_sub(&f, &d, &c); /* F = D-C */
    fe25519_add(&g, &d, &c); /* G = D+C */

    fe25519_mul(&r->x, &e, &f); /* X3 = E*F */
    fe25519_mul(&r->y, &h, &g); /* Y3 = G*H */
    fe25519_mul(&r->z, &g, &f); /* Z3 = F*G */
    fe25519_mul(&r->t, &e, &h); /* T3 = E*H */
}

/* Add `p` and `q`, and store the result in `r`. See [0]. Note how this does
 * not explicitly calculate (X, Y, Z, T), because the result should be a
 * "completed point".
 * [0] https://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-add-2008-hwcd-3
 */
static void add_p1p1(
    ge25519_p1p1 *r,
    const ge25519_p3 *p,
    const ge25519_p3 *q
)
{
    fe25519 a, b, c, d, t;

    /* A = (Y1-X1)*(Y2-X2) */
    fe25519_sub(&a, &p->y, &p->x); /* t0 = Y1-X1 */
    fe25519_sub(&t, &q->y, &q->x); /* t1 = Y2-X2 */
    fe25519_mul(&a, &a, &t); /* A = t0*t1 */

    /* B = (Y1+X1)*(Y2+X2) */
    fe25519_add(&b, &p->x, &p->y); /* t2 = Y1+X1 */
    fe25519_add(&t, &q->x, &q->y); /* t3 = Y2+X2 */
    fe25519_mul(&b, &b, &t); /* B = t2*t3 */

    /* C = T1*k*T2 */
    fe25519_mul(&c, &p->t, &q->t); /* t4 = T1*T2 (different from EFD) */
    fe25519_mul(&c, &c, &ge25519_ec2d); /* C = t4*k (different from EFD) */

    /* D = Z1*2*Z2 */
    fe25519_mul(&d, &p->z, &q->z); /* t5 = Z1*Z2 (different from EFD) */
    fe25519_add(&d, &d, &d); /* D = 2*t5 (different from EFD) */

    fe25519_sub(&r->x, &b, &a); /* E = B-A */
    fe25519_sub(&r->t, &d, &c); /* F = D-C */
    fe25519_add(&r->z, &d, &c); /* G = D+C */
    fe25519_add(&r->y, &b, &a); /* H = B+A */
}

/* Douple `p`, and store the result in `r`. See [0]. Note that a=-1, so
 * multiplying with a is a negation. Also note how this does not explicitly
 * calculate (X, Y, T, Z), because the result should be a "completed point".
 * [0] http://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#doubling-dbl-2008-hwcd */
static void dbl_p1p1(
    ge25519_p1p1 *r,
    const ge25519_p2 *p
)
{
    fe25519 a, b, c, d;

    fe25519_square(&a, &p->x); /* A = X1^2 */
    fe25519_square(&b, &p->y); /* B = Y1^2 */

    /* C = 2*Z1^2 */
    fe25519_square(&c, &p->z); /* t0 = Z1^2 */
    fe25519_add(&c, &c, &c); /* C = 2*t0 */

    fe25519_neg(&d, &a); /* D = a*A */

    /* E = (X1+Y1)^2-A-B */
    fe25519_add(&r->x, &p->x, &p->y); /* t1 = X1+Y1 */
    fe25519_square(&r->x, &r->x); /* t2 = t1^2 */
    fe25519_sub(&r->x, &r->x, &a); /* t3 = t2-A */
    fe25519_sub(&r->x, &r->x, &b); /* E = t3-B */

    fe25519_add(&r->z, &d, &b); /* G = D+B */
    fe25519_sub(&r->t, &r->z, &c); /* F = G-C */
    fe25519_sub(&r->y, &d, &b); /* H = D-B */
}

/* If `b` is `1`, moves `p` into `r`. For `b` equal to `0`, `r` is left
 * unchanged. Hence the name "conditional move". This function runs in constant
 * time. */
static void cmov_aff(
    ge25519_aff *r,
    const ge25519_aff *p,
    unsigned char b
)
{
    fe25519_cmov(&r->x, &p->x, b);
    fe25519_cmov(&r->y, &p->y, b);
}

/* Whether or not `b` is equal to `c`. */
static unsigned char equal(
    signed char b,
    signed char c
)
{
    unsigned char ub = b;
    unsigned char uc = c;

    /* If b = c, then x = 0. Otherwise, x is in [1,255]. */
    unsigned char x = ub ^ uc;

    /* Next we want to treat x as an integer. */
    crypto_uint32 y = x;

    /* If b = c, then y = 4294967295. Otherwise, y is in [0,254]. */
    y -= 1;

    /* If b = c, then y = 1. Otherwise, y = 0. */
    y >>= 31;

    return y;
}

/* Whether or not `b` is negative. */
static unsigned char negative(
    signed char b
)
{
    /* If b < 0, then x is in [18446744073709551361,18446744073709551615]. Otherwise, x is in [0,255]. */
    /* TODO: Why do we use a long long here? */
    unsigned long long x = b;

    /* If b < 0, then x = 1. Otherwise, x = 0. */
    x >>= 63;

    return x;
}

/* Pick the `b`th base multiple for position `pos` from a pre-computed lookup
 * table. Here, the position refers to the index of the digit that we currently
 * process. See ge25519_scalarmult_base() for details. This function runs in
 * constant time. */
static void choose_t(
    ge25519_aff *t,
    unsigned long long pos,
    signed char b
)
{
    /* Since we want this function to run in constant time, we cannot use a
     * conditional here. Instead, we always initialize the result buffer with
     * the 0th multiple, and conditionally overwrite this with another value in
     * the next step. Note that the 0th multiple will always be (0, 1) (the
     * neutral element), but we cannot skip this step due to the constant time
     * requirement. */
    *t = ge25519_base_multiples_affine[5 * pos + 0];

    /* These calls conditionally overwrite the 0th multiple. The condition here
     * is b matches the multiplicity of a certain digit. We can disregard the
     * sign in this step, and adjust the sign in the next step. Under the hood,
     * cmov_aff() uses a bitmask to ensure the constant time requirement. */
    cmov_aff(t, &ge25519_base_multiples_affine[5 * pos + 1], equal(b, 1) | equal(b, -1));
    cmov_aff(t, &ge25519_base_multiples_affine[5 * pos + 2], equal(b, 2) | equal(b, -2));
    cmov_aff(t, &ge25519_base_multiples_affine[5 * pos + 3], equal(b, 3) | equal(b, -3));
    cmov_aff(t, &ge25519_base_multiples_affine[5 * pos + 4], equal(b, -4));

    /* In the last step, we disregarded the sign, so the result for b = 1 and
     * b = -1 will be the same. We now have to negate the x coordinate if b is
     * negative. */
    fe25519 v;
    fe25519_neg(&v, &t->x);
    fe25519_cmov(&t->x, &v, negative(b));
}

/* Set `r` to the neutral element of the group. The neutral element for Ed25519
 * is (0,1). Since we use extended homogeneous coordinates (X, Y, Z, T), with
 * x = X/Z, y = Y/Z, x * y = T/Z, we can set Z = 1. */
static void setneutral(
    ge25519 *r
)
{
    /* Set r = (0, 1, 1, 0). */
    fe25519_setzero(&r->x);
    fe25519_setone(&r->y);
    fe25519_setone(&r->z);
    fe25519_setzero(&r->t);
}

/* return 0 on success, -1 otherwise */
/* TODO: Add explanation for this function. */
int ge25519_unpackneg_vartime(
    ge25519_p3 *r,
    const unsigned char p[32]
)
{
    unsigned char par;
    fe25519 t, chk, num, den, den2, den4, den6;
    fe25519_setone(&r->z);
    par = p[31] >> 7;

    /* fe25519_unpack() disregards the most-significant bit. In packed
     * representation, that bit is used to store the sign of the x coordinate. */
    fe25519_unpack(&r->y, p);
    fe25519_square(&num, &r->y);        /* x = y^2 */
    fe25519_mul(&den, &num, &ge25519_ecd);      /* den = dy^2 */
    fe25519_sub(&num, &num, &r->z);     /* x = y^2-1 */
    fe25519_add(&den, &r->z, &den);     /* den = dy^2+1 */

    /* Computation of sqrt(num/den) */
    /* 1.: computation of num^((p-5)/8)*den^((7p-35)/8) = (num*den^7)^((p-5)/8) */
    fe25519_square(&den2, &den);
    fe25519_square(&den4, &den2);
    fe25519_mul(&den6, &den4, &den2);
    fe25519_mul(&t, &den6, &num);
    fe25519_mul(&t, &t, &den);

    fe25519_pow2523(&t, &t);
    /* 2. computation of r->x = t * num * den^3 */
    fe25519_mul(&t, &t, &num);
    fe25519_mul(&t, &t, &den);
    fe25519_mul(&t, &t, &den);
    fe25519_mul(&r->x, &t, &den);

    /* 3. Check whether sqrt computation gave correct result, multiply by sqrt(-1) if not: */
    fe25519_square(&chk, &r->x);
    fe25519_mul(&chk, &chk, &den);
    if (!fe25519_iseq_vartime(&chk, &num))
        fe25519_mul(&r->x, &r->x, &ge25519_sqrtm1);

    /* 4. Now we have one of the two square roots, except if input was not a square */
    fe25519_square(&chk, &r->x);
    fe25519_mul(&chk, &chk, &den);
    if (!fe25519_iseq_vartime(&chk, &num))
        return -1;

    /* 5. Choose the desired square root according to parity: */
    if (fe25519_getparity(&r->x) != (1 - par))
        fe25519_neg(&r->x, &r->x);

    fe25519_mul(&r->t, &r->x, &r->y);
    return 0;
}

/* Encode a point as described in RFC 8032, section 5.1.2. */
void ge25519_pack(
    unsigned char r[32],
    const ge25519_p3 *p
)
{
    fe25519 tx, ty, zi;

    /* Convert the extended coordinates back into affine representation as
     * described in [0].
     * [0] https://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html */
    fe25519_invert(&zi, &p->z);
    fe25519_mul(&tx, &p->x, &zi);
    fe25519_mul(&ty, &p->y, &zi);

    /* Encode the y-coordinate as a 32-octet string. */
    fe25519_pack(r, &ty);

    /* Copy the least significant bit of the x-coordinate to the most
     * significant bit of the result. Given only the y-coordinate, there would
     * be two possible points on the curve. Therefore, we have to store the
     * sign of the x-coordinate in the encoding of the point. */
    r[31] ^= fe25519_getparity(&tx) << 7;
}

/* Compute [s1]p1 + [s2]p2.
 * TODO: Is this the Bos–Coster method? */
/* TODO: Add explanation for this function. */
void ge25519_double_scalarmult_vartime(
    ge25519_p3 *r,
    const ge25519_p3 *p1,
    const sc25519 *s1,
    const ge25519_p3 *p2,
    const sc25519 *s2
)
{
    ge25519_p1p1 tp1p1;
    ge25519_p3 pre[16];
    unsigned char b[127];
    int i;

    /* precomputation s2 s1 */
    setneutral(pre);            /* 00 00 */
    pre[1] = *p1;               /* 00 01 */
    dbl_p1p1(&tp1p1, (ge25519_p2 *) p1);
    p1p1_to_p3(&pre[2], &tp1p1);        /* 00 10 */
    add_p1p1(&tp1p1, &pre[1], &pre[2]);
    p1p1_to_p3(&pre[3], &tp1p1);        /* 00 11 */
    pre[4] = *p2;               /* 01 00 */
    add_p1p1(&tp1p1, &pre[1], &pre[4]);
    p1p1_to_p3(&pre[5], &tp1p1);        /* 01 01 */
    add_p1p1(&tp1p1, &pre[2], &pre[4]);
    p1p1_to_p3(&pre[6], &tp1p1);        /* 01 10 */
    add_p1p1(&tp1p1, &pre[3], &pre[4]);
    p1p1_to_p3(&pre[7], &tp1p1);        /* 01 11 */
    dbl_p1p1(&tp1p1, (ge25519_p2 *) p2);
    p1p1_to_p3(&pre[8], &tp1p1);        /* 10 00 */
    add_p1p1(&tp1p1, &pre[1], &pre[8]);
    p1p1_to_p3(&pre[9], &tp1p1);        /* 10 01 */
    dbl_p1p1(&tp1p1, (ge25519_p2 *) & pre[5]);
    p1p1_to_p3(&pre[10], &tp1p1);       /* 10 10 */
    add_p1p1(&tp1p1, &pre[3], &pre[8]);
    p1p1_to_p3(&pre[11], &tp1p1);       /* 10 11 */
    add_p1p1(&tp1p1, &pre[4], &pre[8]);
    p1p1_to_p3(&pre[12], &tp1p1);       /* 11 00 */
    add_p1p1(&tp1p1, &pre[1], &pre[12]);
    p1p1_to_p3(&pre[13], &tp1p1);       /* 11 01 */
    add_p1p1(&tp1p1, &pre[2], &pre[12]);
    p1p1_to_p3(&pre[14], &tp1p1);       /* 11 10 */
    add_p1p1(&tp1p1, &pre[3], &pre[12]);
    p1p1_to_p3(&pre[15], &tp1p1);       /* 11 11 */

    sc25519_2interleave2(b, s1, s2);

    /* scalar multiplication */
    *r = pre[b[126]];
    for (i = 125; i >= 0; i--) {
        dbl_p1p1(&tp1p1, (ge25519_p2 *) r);
        p1p1_to_p2((ge25519_p2 *) r, &tp1p1);
        dbl_p1p1(&tp1p1, (ge25519_p2 *) r);
        if (b[i] != 0) {
            p1p1_to_p3(r, &tp1p1);
            add_p1p1(&tp1p1, r, &pre[b[i]]);
        }
        if (i != 0)
            p1p1_to_p2((ge25519_p2 *) r, &tp1p1);
        else
            p1p1_to_p3(r, &tp1p1);
    }
}

/* Compute [s]g, where g is the base point. Note that this function runs in
 * constant time. */
void ge25519_scalarmult_base(
    ge25519_p3 *r,
    const sc25519 *s
)
{
    signed char b[85];
    int i;
    ge25519_aff t;

    /* First, we transform s into a new representation. In this representation,
     * we use radix 8, but our digits are signed. Since our numbers are
     * representable in 255 bits, we need 255 / 3 = 85 digits. */
    sc25519_window3(b, s);

    /* Using the pre-computed lookup table, we initialize r with the multiple
     * for the first digit. */
    choose_t((ge25519_aff *) r, 0, b[0]);

    /* choose_t() does not set z and t coordinates. Be reminded that we use
     * extended homogeneous coordinates (X, Y, Z, T), with x = X/Z, y = Y/Z,
     * x * y = T/Z, which allows us to set Z = 1. */
    fe25519_setone(&r->z);
    fe25519_mul(&r->t, &r->x, &r->y);

    /* Since we processed the first digit, there are now 84 digits remaining.
     * Note that our desired result is simply the sum over all multiples for
     * the digits in b. For details, see section 4 in "High-speed high-security
     * signatures". */
    for (i = 1; i < 85; i++) {
        /* Pick the next base multiple from the lookup table. */
        choose_t(&t, (unsigned long long)i, b[i]);

        /* Add it to the accumulating result. */
        ge25519_mixadd2(r, &t);
    }
}
