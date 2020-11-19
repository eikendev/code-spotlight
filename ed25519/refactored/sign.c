#include <string.h>

#include "crypto_sign.h"
#include "crypto_hash_sha512.h"
#include "ge25519.h"

/* Implements Sign as described in RFC 8032, section 5.1.6. */

int crypto_sign(
    unsigned char *signed_message,
    unsigned long long *signed_message_len,
    const unsigned char *message,
    const unsigned long long message_len,
    const unsigned char *private_key
)
{
    /* Sign: Step 1 */
    /* Key Generation: Step 1 */

    unsigned char h[64];
    crypto_hash_sha512(h, private_key, 32);

    /* Key Generation: Step 2 */

    h[0] &= 248;
    h[31] &= 127;
    h[31] |= 64;

    sc25519 secret_scalar; // 32 bytes
    sc25519_from32bytes(&secret_scalar, h);

    unsigned char *prefix = h + 32; // 32 bytes

    /* Key Generation: Step 3 and Step 4 */

    /* The public key A can be copied from the arguments. Instead of
     * calculating this each time, it's safe to assume the user has some place
     * to store this. */
    unsigned char public_key[32];
    memmove(public_key, private_key + 32, 32);

    /* Sign: Step 2 */

    /* We use the signed_message buffer to calculate the nonce. Note that at
     * this point, the name signed_message is misleading, since the prefix will
     * not be part of the signed message. */
    memmove(signed_message + 64, message, message_len);
    memmove(signed_message + 32, prefix, 32);

    unsigned char r_digest[64];
    crypto_hash_sha512(r_digest, signed_message + 32, message_len + 32);

    sc25519 r; // 64 bytes
    sc25519_from64bytes(&r, r_digest);

    /* Sign: Step 3 */

    ge25519 rB; // [r]B
    ge25519_scalarmult_base(&rB, &r);

    /* Encode [r]B as R and move R into signed_message. Note that the moving is
     * not actually part of Step 3, but done here for efficiency reasons. */
    ge25519_pack(signed_message, &rB);

    /* Sign: Step 4 */

    /* We want to calculate SHA512(dom2(F, C) || R || A || PH(M)). dom2(F, C)
     * is empty for plain Ed25519, and PH is the identity function. Thus, we
     * still need to move the public key A into the buffer, so we have
     * (R || A || M). */
    memmove(signed_message + 32, public_key, 32);

    unsigned char hram[64];
    crypto_hash_sha512(hram, signed_message, message_len + 64);

    sc25519 k; // 64 bytes
    sc25519_from64bytes(&k, hram);

    /* Sign: Step 5 */

    sc25519 *S = &k; // 64 bytes
    sc25519_mul(S, &k, &secret_scalar);
    sc25519_add(S, S, &r);

    /* Sign: Step 6 */

    /* The signed message is (R || S || M). We are left with writing S into the
     * buffer. */
    sc25519_to32bytes(signed_message + 32, S);

    /* Write the output length. This is not described in the RFC, since it's
     * implementation specific. */
    *signed_message_len = message_len + 64;

    return 0;
}
