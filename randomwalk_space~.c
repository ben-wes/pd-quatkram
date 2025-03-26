/*
randomwalk_space~ - Continuous random walk in 3D space at signal rate
2024, Ben Wesch

Functionality:
- Performs continuous random walk in 3D space at signal rate
- Step size and angle range can be controlled via signals
- Outputs both Cartesian (x,y,z) and spherical (theta,phi) coordinates as signals

Inlets:
1. Signal controls step size (in turns/sample, where 1.0 = full rotation)
2. Signal controls angle range variation (in turns, where 1.0 = full rotation)

Outlets:
1. 3-channel signal [x y z] - Cartesian coordinates
2. Signal theta (longitude, degrees)
3. Signal phi (latitude, degrees)

Creation arguments: 
1. step size (turns/sample), default = 0.001 (0.36° per sample)
2. angle range (turns), default = 0.5 (±180 degrees)
3. start_x, default = 0.0
4. start_y, default = 0.0
5. start_z, default = 1.0 (north pole)

A step size of 0.001 means:
- At 44.1kHz: ~15.876 full rotations per second
- At 48kHz: ~17.28 full rotations per second
*/

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Define UNUSED macro for parameters that are intentionally not used
#define UNUSED(x) ((void)(x))

// Define PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Vector3 structure (now using double)
typedef struct _Vector3 {
    double x;
    double y;
    double z;
} Vector3;

// Quaternion structure (now using double)
typedef struct _Quaternion {
    double w;
    double x;
    double y;
    double z;
} Quaternion;

// Main PD object structure
static t_class *randomwalk_space_tilde_class;

typedef struct _randomwalk_space_tilde {
    t_object x_obj;
    t_sample f;              // For first signal inlet (step size)
    Vector3 position;        // Current position in space
    Vector3 step_axis;       // Current step axis
    double angle_offset;     // Angle offset in turns (wrapped to [-0.5, 0.5])
    t_sample **vec_out;      // 3-channel position output
    t_sample *zero_buffer;   // Dynamic zero buffer
    int buffer_size;         // Size of zero buffer
    Vector3 init_position;   // Store initial position for reset
    unsigned int x_state;    // Random state (renamed to match PD convention)
} t_randomwalk_space_tilde;

// Function prototypes
void randomwalk_space_tilde_seed(t_randomwalk_space_tilde *x, t_floatarg f);
void *randomwalk_space_tilde_new(t_symbol *s, int argc, t_atom *argv);
void randomwalk_space_tilde_free(t_randomwalk_space_tilde *x);
static t_int *randomwalk_space_tilde_perform(t_int *w);
static void randomwalk_space_tilde_dsp(t_randomwalk_space_tilde *x, t_signal **sp);
static double random_float(t_randomwalk_space_tilde *x);
static Vector3 generate_random_perpendicular_to_direction(t_randomwalk_space_tilde *x, Vector3 *direction);

// Vector operations
float vector3_magnitude(Vector3 *v) {
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

void vector3_normalize(Vector3 *v) {
    float mag = vector3_magnitude(v);
    if (mag > 0.0001f) {
        v->x /= mag;
        v->y /= mag;
        v->z /= mag;
    }
}

Vector3 vector3_cross(Vector3 *a, Vector3 *b) {
    Vector3 result;
    result.x = a->y * b->z - a->z * b->y;
    result.y = a->z * b->x - a->x * b->z;
    result.z = a->x * b->y - a->y * b->x;
    return result;
}

float vector3_dot(Vector3 *a, Vector3 *b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

Vector3 vector3_scale(Vector3 *v, float scalar) {
    Vector3 result;
    result.x = v->x * scalar;
    result.y = v->y * scalar;
    result.z = v->z * scalar;
    return result;
}

Vector3 vector3_add(Vector3 *a, Vector3 *b) {
    Vector3 result;
    result.x = a->x + b->x;
    result.y = a->y + b->y;
    result.z = a->z + b->z;
    return result;
}

// Quaternion operations
Quaternion quaternion_from_axis_angle(Vector3 *axis, float angle_radians) {
    Quaternion q;
    float half_angle = angle_radians * 0.5f;
    float s = sin(half_angle);
    
    q.w = cos(half_angle);
    q.x = axis->x * s;
    q.y = axis->y * s;
    q.z = axis->z * s;
    
    return q;
}

Quaternion quaternion_multiply(Quaternion *a, Quaternion *b) {
    Quaternion result;
    
    result.w = a->w * b->w - a->x * b->x - a->y * b->y - a->z * b->z;
    result.x = a->w * b->x + a->x * b->w + a->y * b->z - a->z * b->y;
    result.y = a->w * b->y - a->x * b->z + a->y * b->w + a->z * b->x;
    result.z = a->w * b->z + a->x * b->y - a->y * b->x + a->z * b->w;
    
    return result;
}

Quaternion quaternion_conjugate(Quaternion *q) {
    Quaternion result;
    result.w = q->w;
    result.x = -q->x;
    result.y = -q->y;
    result.z = -q->z;
    return result;
}

// Rotate a vector using a quaternion
Vector3 quaternion_rotate_vector(Quaternion *q, Vector3 *v) {
    // Create a quaternion from the vector (with w=0)
    Quaternion v_quat;
    v_quat.w = 0;
    v_quat.x = v->x;
    v_quat.y = v->y;
    v_quat.z = v->z;
    
    // Perform the quaternion rotation: q * v * q^-1
    Quaternion q_conj = quaternion_conjugate(q);
    Quaternion temp = quaternion_multiply(q, &v_quat);
    Quaternion result = quaternion_multiply(&temp, &q_conj);
    
    // Extract the vector part
    Vector3 rotated;
    rotated.x = result.x;
    rotated.y = result.y;
    rotated.z = result.z;
    
    return rotated;
}

// Generate a random perpendicular vector using permutation trick
Vector3 generate_perpendicular_vector(Vector3 *position) {
    Vector3 norm_pos = *position;
    vector3_normalize(&norm_pos);
    
    // Simple permutation trick: (x,y,z) -> (-y,x,0)
    Vector3 perp;
    perp.x = -norm_pos.y;
    perp.y = norm_pos.x;
    perp.z = 0;
    
    // Handle edge case if we're near the poles
    if (vector3_magnitude(&perp) < 0.001f) {
        perp.x = 1.0f;
        perp.y = 0.0f;
        perp.z = 0.0f;
    }
    
    vector3_normalize(&perp);
    return perp;
}

// Move the implementation up before generate_random_perpendicular_to_direction
static double random_float(t_randomwalk_space_tilde *x)
{
    unsigned int randval = x->x_state;
    x->x_state = randval = randval * 472940017 + 832416023;
    return ((double)randval * (1./4294967296.) * 2.0 - 1.0);
}

// Then our new function
static Vector3 generate_random_perpendicular_to_direction(t_randomwalk_space_tilde *x, Vector3 *direction) {
    // Get a base perpendicular vector using the same permutation trick
    Vector3 perp;
    perp.x = -direction->y;
    perp.y = direction->x;
    perp.z = 0;
    
    // Handle edge case if we're near vertical
    if (vector3_magnitude(&perp) < 0.001f) {
        perp.x = 1.0f;
        perp.y = 0.0f;
        perp.z = 0.0f;
    }
    
    vector3_normalize(&perp);
    
    // Rotate this perpendicular vector around the direction by a random angle
    // to get a random perpendicular direction
    double random_angle = random_float(x) * 2.0 * M_PI;
    Quaternion q = quaternion_from_axis_angle(direction, random_angle);
    return quaternion_rotate_vector(&q, &perp);
}

// Convert position to spherical coordinates
void cartesian_to_spherical(Vector3 *pos, float *theta, float *phi) {
    float r = vector3_magnitude(pos);
    *theta = atan2(pos->y, pos->x);
    *phi = acos(pos->z / r);
    
    // Convert to degrees
    *theta = *theta * 180.0f / M_PI;
    *phi = *phi * 180.0f / M_PI;
    
    // Ensure theta is in [0, 360)
    if (*theta < 0) {
        *theta += 360.0f;
    }
}

// Utility to wrap a value to [-0.5, 0.5]
static double wrap_half(double x) {
    x = fmod(x + 0.5, 1.0) - 0.5;
    return x < -0.5 ? x + 1.0 : x;
}

// Add the angle offset method
static void randomwalk_space_tilde_offset(t_randomwalk_space_tilde *x, t_floatarg f)
{
    x->angle_offset = wrap_half((double)f);
}

// Add PD's seed generation function
static int makeseed(void)
{
    static unsigned int random_nextseed = 1489853723;
    random_nextseed = random_nextseed * 435898247 + 938284287;
    return (random_nextseed & 0x7fffffff);
}

// Then modify the perform routine
static t_int *randomwalk_space_tilde_perform(t_int *w)
{
    t_randomwalk_space_tilde *x = (t_randomwalk_space_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *stepsize_in = (t_sample *)(w[3]);
    t_sample *anglerange_in = (t_sample *)(w[4]);
    t_sample *ox = x->vec_out[0], *oy = x->vec_out[1], *oz = x->vec_out[2];
    
    while (n--) {
        double step_size = *stepsize_in++;  // Now in actual units, not turns
        double angle_range = *anglerange_in++ * (2.0 * M_PI); // Still in turns for rotation
        
        // Calculate total angle including random variation and offset
        double total_angle = x->angle_offset * 2.0 * M_PI;
        if (angle_range > 0.0) {
            double variation = random_float(x) * angle_range;
            total_angle += variation;
        }

        // Apply rotation to step_axis if we have any angle to apply
        if (total_angle != 0.0) {
            // Generate a random perpendicular to our current direction
            Vector3 rot_axis = generate_random_perpendicular_to_direction(x, &x->step_axis);
            
            // Rotate step_axis around this perpendicular
            Quaternion q = quaternion_from_axis_angle(&rot_axis, total_angle);
            x->step_axis = quaternion_rotate_vector(&q, &x->step_axis);
            vector3_normalize(&x->step_axis);
        }
        
        // Move position along step_axis by step_size
        Vector3 step = vector3_scale(&x->step_axis, step_size);
        x->position = vector3_add(&x->position, &step);
        
        // Output position
        *ox++ = (t_sample)x->position.x;
        *oy++ = (t_sample)x->position.y;
        *oz++ = (t_sample)x->position.z;
    }
    
    return (w+5);
}

static void randomwalk_space_tilde_dsp(t_randomwalk_space_tilde *x, t_signal **sp)
{
    int vec_size = sp[0]->s_n;

    // Reallocate zero buffer if necessary
    if (x->buffer_size != vec_size) {
        if (x->zero_buffer) {
            freebytes(x->zero_buffer, x->buffer_size * sizeof(t_sample));
        }
        x->zero_buffer = (t_sample *)getbytes(vec_size * sizeof(t_sample));
        if (!x->zero_buffer) {
            pd_error(x, "randomwalk_space~: out of memory");
            return;
        }
        x->buffer_size = vec_size;
        for (int i = 0; i < vec_size; i++) {
            x->zero_buffer[i] = 0;
        }
    }

    // Set up multichannel output
    signal_setmultiout(&sp[2], 3);  // 3-channel output
    for (int i = 0; i < 3; i++) {
        x->vec_out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }

    dsp_add(randomwalk_space_tilde_perform, 4,
            x,
            sp[0]->s_n,
            sp[0]->s_vec,      // stepsize signal
            sp[1]->s_vec);     // anglerange signal
}

// Update new to initialize the RNG state
void *randomwalk_space_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    UNUSED(s);
    t_randomwalk_space_tilde *x = (t_randomwalk_space_tilde *)pd_new(randomwalk_space_tilde_class);
    
    // Initialize buffers
    x->vec_out = (t_sample **)getbytes(3 * sizeof(t_sample *));
    x->zero_buffer = NULL;
    x->buffer_size = 0;
    
    // Default values
    float stepsize = (argc > 0) ? atom_getfloatarg(0, argc, argv) : 0.001f;
    float anglerange = (argc > 1) ? atom_getfloatarg(1, argc, argv) : 0.5f;
    
    // Initialize position
    x->position.x = 0.0;
    x->position.y = 0.0;
    x->position.z = 0.0;  // Start at origin by default
    
    // Store initial position
    x->init_position = x->position;
    
    if (argc > 4) {
        x->position.x = atom_getfloatarg(2, argc, argv);
        x->position.y = atom_getfloatarg(3, argc, argv);
        x->position.z = atom_getfloatarg(4, argc, argv);
        x->init_position = x->position;
    }
    
    // Initialize step_axis as a normalized vector (any direction)
    x->step_axis.x = 0.0;
    x->step_axis.y = 0.0;
    x->step_axis.z = 1.0;  // Start moving upward by default
    
    // Initialize angle offset
    x->angle_offset = 0.0;
    
    // Create inlets (left to right)
    x->f = stepsize;  // First inlet default value
    signalinlet_new(&x->x_obj, anglerange);  // Second inlet with default value
    
    // Create outlet
    outlet_new(&x->x_obj, &s_signal);  // 3-channel position output
    
    // Initialize RNG state using PD's method
    x->x_state = makeseed();
    
    return (void *)x;
}

void randomwalk_space_tilde_free(t_randomwalk_space_tilde *x)
{
    freebytes(x->vec_out, 3 * sizeof(t_sample *));
    if (x->zero_buffer) {
        freebytes(x->zero_buffer, x->buffer_size * sizeof(t_sample));
    }
}

// Reset method
static void randomwalk_space_tilde_reset(t_randomwalk_space_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    UNUSED(s);
    
    if (argc >= 3) {
        // Set new position from arguments
        x->position.x = atom_getfloatarg(0, argc, argv);
        x->position.y = atom_getfloatarg(1, argc, argv);
        x->position.z = atom_getfloatarg(2, argc, argv);
    } else {
        // Reset to initial position
        x->position = x->init_position;
    }
    
    // Reset step_axis to default direction
    x->step_axis.x = 0.0;
    x->step_axis.y = 0.0;
    x->step_axis.z = 1.0;
}


// Implementation
void randomwalk_space_tilde_seed(t_randomwalk_space_tilde *x, t_floatarg f)
{
    x->x_state = (unsigned int)f;
}

void randomwalk_space_tilde_setup(void)
{
    randomwalk_space_tilde_class = class_new(gensym("randomwalk_space~"),
        (t_newmethod)randomwalk_space_tilde_new,
        (t_method)randomwalk_space_tilde_free,
        sizeof(t_randomwalk_space_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);
    
    class_addmethod(randomwalk_space_tilde_class, 
                   (t_method)randomwalk_space_tilde_dsp, 
                   gensym("dsp"), 
                   A_CANT, 0);
    class_addmethod(randomwalk_space_tilde_class,
                   (t_method)randomwalk_space_tilde_seed,
                   gensym("seed"),
                   A_FLOAT, 0);
    class_addmethod(randomwalk_space_tilde_class,
                   (t_method)randomwalk_space_tilde_offset,
                   gensym("offset"),
                   A_FLOAT, 0);
    class_addmethod(randomwalk_space_tilde_class,
                   (t_method)randomwalk_space_tilde_reset,
                   gensym("reset"),
                   A_GIMME, 0);
    CLASS_MAINSIGNALIN(randomwalk_space_tilde_class, t_randomwalk_space_tilde, f);
}