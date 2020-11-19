# Reference Implementation of Ed25519

I've recently started reading the reference implementation of [Ed25519](https://ed25519.cr.yp.to/), which is included in the [SUPERCOP](https://bench.cr.yp.to/supercop.html) benchmark toolkit.
It was written by Daniel J. Bernstein, Niels Duif, Tanja Lange, Peter Schwabe, and Bo-Yin Yang.

Even if I think that the primitive itself is beautiful, the code was not easy to read at all!

This is why I went ahead and tried to understand it.
Throughout this process, I've commented the code so it's easier for me if I want to revisit it later.
Also, I want to share it so others can benefit, too!

Note that I'm not fully done documenting the reference implementation.

_This project is only meant to provide help in understanding the implementation.
The intention is by no means to semantically improve it in any way.
I'm not an expert in this field, so the content in this repository should not be seen as professional advice._

## Implementation Background

The high-level functions for signing, verifying and key generation are explained in [RFC 8032](https://tools.ietf.org/html/rfc8032#section-5.1).
Specifically, the code references sections 5.1.5 to 5.1.7.

The Barrett reduction is implemented as explained in the [Handbook of Applied Cryptography](http://cacr.uwaterloo.ca/hac/about/chap14.pdf).
You will find references to the algorithm in the code, so if you are interested in the details, it makes sense to open code and the book side-by-side.

## Changes in Refactored Version

- All functions that were not used at all were removed.
- The defines for integration with SUPERCOP were removed.
- The code was reformatted so it has a consistent styling.
- Some variables were renamed to give them a more meaningful name.
- Comments were added or refined to explain each individual step.
