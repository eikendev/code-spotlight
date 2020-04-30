#include <string.h>

#include "crypto_sign.h"
#include "crypto_hash_sha512.h"
#include "randombytes.h"
#include "ge25519.h"

/* Implements Key Generation as described in RFC 8032, section 5.1.5. */

int crypto_sign_keypair(
    unsigned char *public_key,
    unsigned char *secret_key
)
{
    /* Generate random bytes for the secret key. This is not specified in the
     * RFC, but is needed in the measure function of this package. */
    randombytes(secret_key, 32);

    /* Key Generation: Step 1 */

    unsigned char h[64];
    crypto_hash_sha512(h, secret_key, 32);

    /* Key Generation: Step 2 */

    h[0] &= 248;
    h[31] &= 127;
    h[31] |= 64;

    /* Key Generation: Step 3 */

    sc25519 secret_scalar; // s, 32 bytes
    sc25519_from32bytes(&secret_scalar, h);

    ge25519 sB; // [s]B
    ge25519_scalarmult_base(&sB, &secret_scalar);

    /* Key Generation: Step 4 */

    ge25519_pack(public_key, &sB);

    /* Finally, we append the public key to the private key. This is not
     * specified in the RFC, but is needed in the measure function of this
     * package. */
    memmove(secret_key + 32, public_key, 32);

    return 0;
}
