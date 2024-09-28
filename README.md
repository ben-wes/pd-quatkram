# quatkram Pd library
arbitrary collection of externals and abstractions originated around quaternion transformation for Pd (Pure Data)

| object | description |
| ---- | ----- |
| `[qacc~]` | Quaternion accumulator (quaternion multiplication and normalization of current state with incoming quaternion for each sample) |
| `[qmul~]` | Quaternion multiplication - right quaternion can be defined with creation args. left falls back to identity quaternion if no input provided |
| `[qdiv~]` | Quaternion division - similar to multiplication with inverse quaternion |
| `[qconj~]` | Conjugate quaternion (w -x -y -z) |
| `[qinv~]` | Inverse quaternion (similar to conjugate divided by squared magnitude) |
| `[qnorm~]` | Quaternion norm (= squared magnitude) |
| `[qnormalize~]` | Normalize input quaternion to unit quaternion |
| `[qmag~]` | Magnitude of quaternion |
| `[qfromgyroaccel~]` | Derive quaternion from sensor input (needs more parameters and more efficient implementation) |
| `[qvtrans~]` | Quaternion-based vector transformation - expects 3-channel vector on left inlet and quaternion on right inlet, outputs transformed vector |
| `[faccwrap~]` | Float accumulator and wrapper - adds incoming float to current value for each sample and wraps around `-1..1` (avoiding the problem of losing precision when using vanilla's [rpole~] and [wrap~] for this kind of integration and wrapping) |
| `[faccbounce~]` | Float accumulator with boundary bouncing - similar to [faccwrap~], but bounces back from boundaries `-1 / 1` instead of wrapping |
| `[atan2~]` | Signal rate atan2 implementation |
| `[noisen~]` | Outputs normally distributed values (Gaussian noise) at signal rate. Can be seeded with non-zero values |
| `[mc_conv~]` | Applies convolution across channels of a multichannel input signal using a user-defined kernel (channel-domain, not time-domain) |
| `[mc_conv2d~]` | Applies convolution across channels representing a 2d square grid with input signal (also representing a square grid) |
| `[mc_route~]` | Routing multichannel input with another multichannel routing signal (of same or smaller channel count - in which case the routing pattern gets repeated for the other input channels) |
| `[urn~]` | Urn model random number generator on signal rate, triggered by impulses (or just signal `1`) and seedable with signal inlet (seed applied for each cycle) |
