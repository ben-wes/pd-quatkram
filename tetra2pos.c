#include "m_pd.h"
#include <string.h>
#include <math.h>

static t_class *tetra2pos_class;

typedef struct _tetra2pos {
    t_object x_obj;
    t_float edge_length;  // in mm
    t_float positions[4][3]; // actual positions in mm
    int debug;
    t_outlet *position_out;  // position in mm
} t_tetra2pos;

static void calculate_positions(t_tetra2pos *x) {
    t_float a = x->edge_length / 2;  // half edge length
    t_float sqrt3 = sqrt(3.0f);
    t_float sqrt6 = sqrt(6.0f);
    
    // Front mic
    x->positions[0][0] = 0;
    x->positions[0][1] = a * 2 / sqrt3;
    x->positions[0][2] = -a / sqrt6;
    
    // Left back mic
    x->positions[1][0] = -a;
    x->positions[1][1] = -a / sqrt3;
    x->positions[1][2] = -a / sqrt6;
    
    // Right back mic
    x->positions[2][0] = a;
    x->positions[2][1] = -a / sqrt3;
    x->positions[2][2] = -a / sqrt6;
    
    // Top mic
    x->positions[3][0] = 0;
    x->positions[3][1] = 0;
    x->positions[3][2] = a * 3 / sqrt6;
}

static void tetra2pos_debug(t_tetra2pos *x, t_floatarg f) {
    x->debug = (int)f;
}

static void tetra2pos_edge(t_tetra2pos *x, t_floatarg f) {
    if (f <= 0) {
        pd_error(x, "tetra2pos: edge length must be positive");
        return;
    }
    x->edge_length = f;
    calculate_positions(x);
    
    if (x->debug) {
        post("tetra2pos: edge length set to %.1f mm", x->edge_length);
        post("tetra2pos: mic positions (mm):");
        for (int i = 0; i < 4; i++) {
            post("  %d: %.1f %.1f %.1f", i, x->positions[i][0], x->positions[i][1], x->positions[i][2]);
        }
    }
}

static void solve_linear_system(t_float A[3][3], t_float b[3], t_float result[3]) {
    t_float det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
                
    if (fabs(det) < 0.0001f) {
        result[0] = result[1] = result[2] = 0;
        return;
    }
    
    for (int i = 0; i < 3; i++) {
        t_float temp[3][3];
        memcpy(temp, A, sizeof(temp));
        for (int j = 0; j < 3; j++) {
            temp[j][i] = b[j];
        }
        t_float det_i = temp[0][0]*(temp[1][1]*temp[2][2] - temp[1][2]*temp[2][1])
                      - temp[0][1]*(temp[1][0]*temp[2][2] - temp[1][2]*temp[2][0])
                      + temp[0][2]*(temp[1][0]*temp[2][1] - temp[1][1]*temp[2][0]);
        result[i] = det_i/det;
    }
}

static void solve_position_toa(t_float distances[4], t_float positions[4][3], t_float result[3]) {
    t_float A[3][3] = {{0}};
    t_float b[3] = {0};
    
    for (int i = 0; i < 3; i++) {
        t_float dx = positions[i+1][0] - positions[0][0];
        t_float dy = positions[i+1][1] - positions[0][1];
        t_float dz = positions[i+1][2] - positions[0][2];
        
        t_float d0_sq = distances[0] * distances[0];
        t_float di_sq = distances[i+1] * distances[i+1];
        
        t_float p0_sq = positions[0][0]*positions[0][0] + 
                       positions[0][1]*positions[0][1] + 
                       positions[0][2]*positions[0][2];
        t_float pi_sq = positions[i+1][0]*positions[i+1][0] + 
                       positions[i+1][1]*positions[i+1][1] + 
                       positions[i+1][2]*positions[i+1][2];
        
        A[i][0] = 2 * dx;
        A[i][1] = 2 * dy;
        A[i][2] = 2 * dz;
        b[i] = d0_sq - di_sq - p0_sq + pi_sq;
    }
    
    solve_linear_system(A, b, result);
}

static void solve_position_tdoa(t_float distances[4], t_float positions[4][3], t_float result[3]) {
    t_float tdoa[3];
    for (int i = 0; i < 3; i++) {
        tdoa[i] = distances[i + 1] - distances[0];
    }
    
    t_float min_error = 1e10;
    t_float best_result[3] = {0, 0, 0};
    t_float calc_distances[4] = {0};
    
    // Find maximum absolute TDOA to help with range estimation
    t_float max_abs_tdoa = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(tdoa[i]) > max_abs_tdoa) max_abs_tdoa = fabs(tdoa[i]);
    }
    
    // Start with a reasonable minimum distance that's at least max_abs_tdoa
    t_float r1_start = fmax(1000, max_abs_tdoa);
    t_float r1_end = 50000;      // Maximum reasonable distance
    t_float r1_step = r1_start * 0.1;  // Initial step proportional to start distance
    
    // Multiple passes with decreasing step sizes
    for (int phase = 0; phase < 5; phase++) {
        // Adjust angular resolution based on phase
        t_float angle_step = (phase < 2) ? M_PI / 8 :
                            (phase < 3) ? M_PI / 16 : M_PI / 32;
        
        for (t_float r1 = r1_start; r1 < r1_end; r1 += r1_step) {
            // Try multiple base angles for each radius
            for (t_float angle = 0; angle < M_PI * 2; angle += angle_step) {
                t_float test_distances[4];
                test_distances[0] = r1;
                test_distances[1] = r1 + tdoa[0];
                test_distances[2] = r1 + tdoa[1];
                test_distances[3] = r1 + tdoa[2];
                
                // Skip if any distance is negative or too small
                if (test_distances[1] < max_abs_tdoa || 
                    test_distances[2] < max_abs_tdoa || 
                    test_distances[3] < max_abs_tdoa) {
                    continue;
                }
                
                t_float test_result[3];
                solve_position_toa(test_distances, positions, test_result);
                
                t_float error = 0;
                t_float max_calc_dist = 0;
                
                for (int i = 0; i < 4; i++) {
                    t_float delta_x = test_result[0] - positions[i][0];
                    t_float delta_y = test_result[1] - positions[i][1];
                    t_float delta_z = test_result[2] - positions[i][2];
                    calc_distances[i] = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
                    if (calc_distances[i] > max_calc_dist) max_calc_dist = calc_distances[i];
                }
                
                // Calculate error with both absolute and relative components
                for (int i = 0; i < 3; i++) {
                    t_float calc_tdoa = calc_distances[i+1] - calc_distances[0];
                    t_float abs_error = fabs(calc_tdoa - tdoa[i]);
                    error += abs_error + abs_error / (max_calc_dist + 1.0f);
                }
                
                if (error < min_error) {
                    min_error = error;
                    best_result[0] = test_result[0];
                    best_result[1] = test_result[1];
                    best_result[2] = test_result[2];
                }
            }
        }
        
        if (phase < 3) {
            t_float best_r1 = calc_distances[0];
            t_float window_factor = (phase < 2) ? 5.0f : 3.0f;
            r1_start = best_r1 - r1_step * window_factor;
            r1_end = best_r1 + r1_step * window_factor;
            r1_step *= 0.1;  // Each phase reduces step size by factor of 10
        }
    }
    
    result[0] = best_result[0];
    result[1] = best_result[1];
    result[2] = best_result[2];
}

static void tetra2pos_relative(t_tetra2pos *x, t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    if (argc != 4) {
        pd_error(x, "tetra2pos: relative message expects 4 distances (mm)");
        return;
    }

    t_float distances[4];
    for (int i = 0; i < 4; i++) {
        distances[i] = atom_getfloat(argv + i);
    }

    t_float position[3];
    solve_position_tdoa(distances, x->positions, position);

    // Force copy the values to ensure no optimization issues
    t_float x_pos = position[0];
    t_float y_pos = position[1];
    t_float z_pos = position[2];

    // Create and output the list
    t_atom position_list[3];
    SETFLOAT(&position_list[0], x_pos);
    SETFLOAT(&position_list[1], y_pos);
    SETFLOAT(&position_list[2], z_pos);

    outlet_list(x->position_out, &s_list, 3, position_list);

    if (x->debug) {
        post("tetra2pos: relative distances (mm): %.1f %.1f %.1f %.1f", 
             distances[0], distances[1], distances[2], distances[3]);
        post("tetra2pos: time differences (mm): %.1f %.1f %.1f", 
             distances[1] - distances[0], 
             distances[2] - distances[0], 
             distances[3] - distances[0]);
        post("tetra2pos: position (mm): %.1f %.1f %.1f",
             x_pos, y_pos, z_pos);
    }
}

static void tetra2pos_list(t_tetra2pos *x, t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    if (argc != 4) {
        pd_error(x, "tetra2pos: expect 4 absolute distances (mm)");
        return;
    }

    t_float distances[4];
    for (int i = 0; i < 4; i++) {
        distances[i] = atom_getfloat(argv + i);
    }

    t_float position[3];
    solve_position_toa(distances, x->positions, position);

    t_atom position_list[3];
    SETFLOAT(&position_list[0], position[0]);
    SETFLOAT(&position_list[1], position[1]);
    SETFLOAT(&position_list[2], position[2]);
    outlet_list(x->position_out, &s_list, 3, position_list);

    if (x->debug) {
        post("tetra2pos: absolute distances (mm): %.1f %.1f %.1f %.1f", 
             distances[0], distances[1], distances[2], distances[3]);
        post("tetra2pos: position (mm): %.1f %.1f %.1f",
             position[0], position[1], position[2]);
    }
}

static void tetra2pos_positions(t_tetra2pos *x, t_symbol *s, int argc, t_atom *argv) {
    (void)s;  // Unused parameter
    
    if (argc != 12) {
        pd_error(x, "tetra2pos: positions message expects 12 values (x y z for each mic)");
        return;
    }

    // Update positions array
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            x->positions[i][j] = atom_getfloat(argv + (i * 3 + j));
        }
    }

    if (x->debug) {
        post("tetra2pos: mic positions manually set (mm):");
        for (int i = 0; i < 4; i++) {
            post("  %d: %.1f %.1f %.1f", i, x->positions[i][0], x->positions[i][1], x->positions[i][2]);
        }
    }
}

static void tetra2pos_print(t_tetra2pos *x) {
    post("tetra2pos: current mic positions (mm):");
    for (int i = 0; i < 4; i++) {
        post("  %d: %.1f %.1f %.1f", i, x->positions[i][0], x->positions[i][1], x->positions[i][2]);
    }
}

static void *tetra2pos_new(t_floatarg edge) {
    t_tetra2pos *x = (t_tetra2pos *)pd_new(tetra2pos_class);
    
    x->position_out = outlet_new(&x->x_obj, &s_list);
    x->edge_length = edge > 0 ? edge : 1000.0f;
    x->debug = 0;
    
    calculate_positions(x);
    return (void *)x;
}

void tetra2pos_setup(void) {
    tetra2pos_class = class_new(gensym("tetra2pos"),
        (t_newmethod)tetra2pos_new,
        0,
        sizeof(t_tetra2pos),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0);
    
    class_addlist(tetra2pos_class, tetra2pos_list);
    class_addmethod(tetra2pos_class, (t_method)tetra2pos_debug, gensym("debug"), A_FLOAT, 0);
    class_addmethod(tetra2pos_class, (t_method)tetra2pos_edge, gensym("edge"), A_FLOAT, 0);
    class_addmethod(tetra2pos_class, (t_method)tetra2pos_relative, gensym("relative"), A_GIMME, 0);
    class_addmethod(tetra2pos_class, (t_method)tetra2pos_positions, gensym("positions"), A_GIMME, 0);
    class_addmethod(tetra2pos_class, (t_method)tetra2pos_print, gensym("print"), 0);
}
