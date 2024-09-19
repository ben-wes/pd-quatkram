# quatkram Pd library
arbitrary collection of externals originated around quaternion transformation for Pd (Pure Data)

* `[faccwrap~]` Float accumulator and wrapper - adds incoming float to current value for each sample and wraps around `-1..1` (avoiding the problem of losing precision when using vanilla's [rpole~] and [wrap~] for this kind of integration and wrapping)
* `[faccbounce~]` Float accumulator with boundary bouncing - similar to [faccwrap~], but bounces back from boundaries `-1 / 1` instead of wrapping
* `[qacc~]` Quaternion accumulator (quaternion multiplication and normalization of current state with incoming quaternion for each sample)
* `[qvrot~]` Quaternion-based vector rotation - expects 3-channel vector on left inlet and quaternion on right inlet, outputs transformed vector
* `[noisen~]` Outputs normally distributed values (Gaussian noise) at signal rate. Can be seeded with non-zero values
* `[mc_conv~]` Applies convolution across channels of a multichannel input signal using a user-defined kernel (channel-domain, not time-domain)
