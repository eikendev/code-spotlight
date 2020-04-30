# ref

## About

This is the reference implementation of Ed25519.
It was originally written by Daniel J. Bernstein, Niels Duif, Tanja Lange, Peter Schwabe, and Bo-Yin Yang.

## Implementation

The Barrett reduction is implemented as explained in the [Handbook of Applied Cryptography](http://cacr.uwaterloo.ca/hac/about/chap14.pdf).
You will find references to the algorithm in the code, so if you are interested in the details, it makes sense to open code and the book side-by-side.

The high-level functions for signing, verifying and key generation are explained in [RFC 8032](https://tools.ietf.org/html/rfc8032#section-5.1).
Specifically, the code references sections 5.1.5 to 5.1.7.

## Changes

- All functions that were not used at all were removed.
- The defines for integration with SUPERCOP were removed.
- The code was reformatted so it has a consistent styling.
- Some variables were renamed to give them a more meaningful name.
- Comments were added or refined to explain each individual step.
