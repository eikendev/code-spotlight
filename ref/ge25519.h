#ifndef GE25519_H
#define GE25519_H

#include "fe25519.h"
#include "sc25519.h"

/* The representation using extended coordinates. It must hold that x = X/Z,
 * y = Y/Z, and x*y = T/Z. Note that for our curve, we have a=-1, which allows
 * for faster arithmetic than in the general case. See
 * https://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html. */
typedef struct {
    fe25519 x;
    fe25519 y;
    fe25519 z;
    fe25519 t;
} ge25519;

const ge25519 ge25519_base;

int ge25519_unpackneg_vartime(
    ge25519 *r,
    const unsigned char p[32]
);

void ge25519_pack(
    unsigned char r[32],
    const ge25519 *p
);

void ge25519_double_scalarmult_vartime(
    ge25519 *r,
    const ge25519 *p1,
    const sc25519 *s1,
    const ge25519 *p2,
    const sc25519 *s2
);

void ge25519_scalarmult_base(
    ge25519 *r,
    const sc25519 *s
);

#endif
