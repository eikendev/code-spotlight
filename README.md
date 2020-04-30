# Explaining implementations of Ed25519

I've recently started reading the reference implementation of Ed25519, which is included in the [SUPERCOP](https://bench.cr.yp.to/supercop.html) benchmark toolkit.
Even if I think that the primitive itself is beautiful, the code was not easy to read at all!

This is why I went ahead and tried to understand it.
Throughout this process, I've commented the code so it's easier for me if I want to revisit it later.
Also, I want to share it so others can benefit, too!

Note that I'm not fully done documenting the reference implementation.
Maybe I'll include other implementations of Ed25519 in the future. Stay tuned! :wink:

_This project is only meant to provide help in understanding other implementations.
The intention is by no means to semantically improve existing implementations.
I'm not an expert in this field, so the content in this repository should not be seen as professional advice._
