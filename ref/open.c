#include <string.h>

#include "crypto_sign.h"
#include "crypto_verify_32.h"
#include "crypto_hash_sha512.h"
#include "ge25519.h"

/* Implements Verify as described in RFC 8032, section 5.1.7. */

int crypto_sign_open(
    unsigned char *message,
    unsigned long long *message_length,
    const unsigned char *signed_message, // (R || S || M)
    unsigned long long signed_message_len,
    const unsigned char *public_key // A
)
{
    /* Verify: Step 1 */

    if (signed_message_len < 64)
        goto bad_signature;

    if (signed_message[63] & 224)
        goto bad_signature;

    ge25519 public_key_decoded_neg; // A'
    if (ge25519_unpackneg_vartime(&public_key_decoded_neg, public_key))
        goto bad_signature;

    unsigned char public_key_copy[32];
    memmove(public_key_copy, public_key, 32);

    unsigned char R[32];
    memmove(R, signed_message, 32);

    sc25519 S; // 32 bytes
    sc25519_from32bytes(&S, signed_message + 32);

    /* Verify: Step 2 */

    /* message = (R || S || M) */
    memmove(message, signed_message, signed_message_len);
    /* message = (R || A || M) */
    memmove(message + 32, public_key_copy, 32);

    /* h = SHA512(dom2(F, C) || R || A || PH(M)) */
    /*   = SHA512(R || A || M) */
    unsigned char h[64];
    crypto_hash_sha512(h, message, signed_message_len);

    sc25519 k; // 64 bytes
    sc25519_from64bytes(&k, h);

    /* Verify: Step 3 */

    /* We now check if [S]B = R + [k]A' as specified in the RFC. At the time of
     * writing, I'm not sure if this potentially rejects correct solutions, as
     * the RFC states it as if it's a stronger requirement. Note that we can
     * write the equation as -[k]A' + [S]B = R. */

    /* Compute [k]public_key_decoded_neg + [S]ge25519_base. */
    ge25519 R_check;
    ge25519_double_scalarmult_vartime(&R_check, &public_key_decoded_neg, &k, &ge25519_base, &S);

    unsigned char R_check_encoded[32];
    ge25519_pack(R_check_encoded, &R_check);

    /* Iterate all bytes, and check for equality. */
    int differentbits = crypto_verify_32(R, R_check_encoded);

    if (differentbits != 0)
        goto bad_signature;

    /* Write the output buffer. This is implementation-specific. */
    memmove(message, message + 64, signed_message_len - 64);
    memset(message + signed_message_len - 64, 0, 64);
    *message_length = signed_message_len - 64;

    return 0;

bad_signature:
    *message_length = (unsigned long long) -1;
    memset(message, 0, signed_message_len);

    return -1;
}
