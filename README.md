# quatkram Pd library
arbitrary collection of externals and abstractions originated around quaternion transformation for Pd (Pure Data)

| object | description |
| ---- | ----- |
| `[qacc~]` | Quaternion accumulator (quaternion multiplication and normalization of current state with incoming quaternion for each sample) |
| `[qmul~]` | Quaternion multiplication - right quaternion can be defined with creation args. both inputs fall back to identity quaternion if no input provided |
| `[qdiv~]` | Quaternion division - similar to multiplication with inverse quaternion |
| `[qconj~]` <sup>abs<sup> | Conjugate quaternion (w -x -y -z) |
| `[qinv~]` <sup>abs<sup> | Inverse quaternion (similar to conjugate divided by squared magnitude) |
| `[qnorm~]` <sup>abs<sup> | Quaternion norm (= squared magnitude) |
| `[qnormalize~]` <sup>abs<sup> | Normalize input quaternion to unit quaternion |
| `[qmag~]` <sup>abs<sup> | Magnitude of quaternion |
| `[qfromgyroaccel~]` <sup>abs<sup> | Derive quaternion from sensor input (needs more parameters and more efficient implementation) |
| `[qvtrans~]` | Quaternion-based vector transformation - expects 3-channel vector on left inlet and quaternion on right inlet, outputs transformed vector |
| `[faccwrap~]` | Float accumulator and wrapper with multichannel support - adds incoming float to current value for each sample and wraps around a defined range (default `-1..1`) to avoid the problem of losing precision when using vanilla's [rpole~] and [wrap~] for this kind of integration and wrapping |
| `[faccbounce~]` | Float accumulator with boundary bouncing with multichannel support - similar to [faccwrap~], but bounces back from defined boundaries (default `-1 / 1`) instead of wrapping, actual direction state change can be achieved with `-r` arg on initization |
| `[faccleak~]` | Float accumulator with leaky behavior (leaky integrator) |
| `[atan2~]` | Signal rate atan2 implementation |
| `[noisen~]` | Outputs normally distributed values (Gaussian noise) at signal rate. Can be seeded with non-zero values |
| `[nchans~]` | Outputs channels of multichannel signal on dsp activation and `bang` input |
| `[mc_conv~]` | Applies convolution across channels of a multichannel input signal using a user-defined kernel (channel-domain, not time-domain) |
| `[mc_conv2d~]` | Applies convolution across channels representing a 2d square grid with input signal (also representing a square grid) |
| `[mc_route~]` | Routing multichannel input with another multichannel routing signal (of same or smaller channel count - in which case the routing pattern gets repeated for the other input channels) |
| `[urn~]` | Urn model random number generator on signal rate, triggered by impulses (or just signal `1`) and seedable with signal inlet (seed applied for each cycle) |
| `[tetra2pos]` | derive location of point based on 4 distance measurements from 4 tetrahedral vertex coordinates |
| `[tetra2pos_abs]` <sup>abs<sup> | see above, vanilla abstraction |
| `[tabsmear~]` | smear into table at arbitrary speed |
| `[tabredraw~]` | refresh garray display |
| `[tabloop~]` | loop table in arbitrary segment (interpolated) |
| `[sampdel~]` | multichannel delay signal by n samples |
| `[frft~]` | fractional fourier transform, applied on block |
| `[frft]` | fractional fourier transform, applied on list input |
| `[zcflip~]` | zero-crossing signal flipper (flips on right inlet impulse after zero-crossing was found for silent signal and zero-crossing is met for audible signal |
| `[zc~]` | multichannel capable zerocrossing detection (outputs impulses on left outlet for upward and on right outlet for downward zero crossings) |
| `[zcflip~]` | zero-crossing signal flipper (flips on right inlet impulse after zero-crossing was found for silent signal and zero-crossing is met for audible signal |
