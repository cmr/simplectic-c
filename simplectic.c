#include "simplectic.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define FIXED(name, type, size) struct name { type elts[size]; }
#define SKEW 0.36602540378 // 0.5 * (sqrt(3.0) - 1.0)
#define UNSKEW 0.36602540378 // (3.0 - sqrt(3.0)) / 6.0
#define SIMPLEX_SIZE 0.70710678119
#define INV_SIMPLEX_SIZE 1.41421356235
#define LAYER_OFFSET_X 0.45534180126
#define LAYER_OFFSET_Y 0.12200846793
#define LAYER_OFFSET_Z 0.35355339059
#define NORM2 8.0
#define NORM3 9.0
#define NORM4 10.0

typedef struct simplectic_point2 {
    int64_t x_cell, y_cell;
    double x_offset, y_offset;
} simplectic_point2;

typedef struct simplectic_point3 {
    int64_t x_cell, y_cell, z_cell;
    double x_offset, y_offset, z_offset;
} simplectic_point3;

typedef struct simplectic_point4 {
    int64_t x_cell, y_cell, z_cell, w_cell;
    double x_offset, y_offset, z_offset, w_offset;
} simplectic_point4;

/*
typedef struct simplectic_pointn {
    const int64_t *cells;
    const double *offsets;
    ssize_t num_dims;
} simplectic_pointn;
*/

struct xorshift {
    uint32_t x, y, z, w;
};

uint32_t xorshift_next(struct xorshift *xs) {
    uint32_t x = xs->x, w = xs->w, t = x ^ (x << 11);
    xs->x = xs->y;
    xs->y = xs->z;
    xs->z = xs->w;
    xs->w = w ^ (w >> 19) ^ (t ^ (t >> 8));
    return xs->w;
}

void simplectic_seed_fill(simplectic_seed *sd, uint32_t seed) {
    for (int i = 0; i < TABLE_SIZE; i++) {
        sd->table[i] = (uint8_t) i;
    }
    // shuffle the table
    struct xorshift xs = {1, seed, seed, seed};
    for (uint32_t i = TABLE_SIZE; i >= 2; i--) {
        uint32_t swap_idx = xorshift_next(&xs) % (i + 1), val = sd->table[swap_idx];
        sd->table[swap_idx] = sd->table[i];
        sd->table[i] = (uint8_t) val;
    }
    memcpy(&sd->table[TABLE_SIZE], &sd->table[0], TABLE_SIZE);
}

FIXED(simplectic2_points, simplectic_point2, 3);

static void simplectic2_points(const double p[static 2], struct simplectic2_points *out) {
    double skew_offset = (p[0] + p[1]) * SKEW,
           x_cell = floor(p[0] + skew_offset),
           y_cell = floor(p[1] + skew_offset),
           unskew_offset = (x_cell + y_cell) * UNSKEW,
           x_origin = x_cell - unskew_offset,
           y_origin = y_cell - unskew_offset,
           dx0 = p[0] - x_origin,
           dy0 = p[1] - y_origin,
           x1_offset = dx0 > dy0 ? 1 : 0,
           y1_offset = dx0 > dy0 ? 0 : 1,
           dx1 = dx0 - x1_offset + UNSKEW,
           dy1 = dy0 - y1_offset + UNSKEW,
           dx2 = dx0 - 1 + 2 * UNSKEW,
           dy2 = dy0 - 1 + 2 * UNSKEW;

    out->elts[0].x_cell = (int64_t) x_cell;
    out->elts[0].y_cell = (int64_t) y_cell;
    out->elts[0].x_offset = dx0;
    out->elts[0].x_offset = dy0;

    out->elts[1].x_cell = (int64_t) (x_cell + x1_offset);
    out->elts[1].y_cell = (int64_t) (y_cell + y1_offset);
    out->elts[1].x_offset = dx1;
    out->elts[1].x_offset = dy1;

    out->elts[2].x_cell = 1 + (int64_t) x_cell;
    out->elts[2].y_cell = 1 + (int64_t) y_cell;
    out->elts[2].x_offset = dx2;
    out->elts[2].x_offset = dy2;
}

static void expand2to3(const simplectic_point2 *p, int64_t cell, double offset, simplectic_point3 *out) {
    out->x_cell = p->x_cell;
    out->x_offset = p->x_offset;

    out->y_cell = p->y_cell;
    out->y_offset = p->y_offset;

    out->z_cell = cell;
    out->z_offset = offset;
}

FIXED(simplectic3_points, simplectic_point3, 6);

static void simplectic3_points(const double point[3], struct simplectic3_points *out) {
    double layer = floor(point[2] * INV_SIMPLEX_SIZE),
           z_offset = point[2] - layer * SIMPLEX_SIZE;
    int64_t layer_int = (int64_t)layer;
    double layer1_point[2], layer2_point[2];

    if (layer_int % 2 == 0) {
        layer1_point[0] = point[0];
        layer1_point[1] = point[1];
        layer2_point[0] = point[0] + LAYER_OFFSET_X;
        layer2_point[1] = point[1] + LAYER_OFFSET_Y;
    } else {
        layer2_point[0] = point[0];
        layer2_point[1] = point[1];
        layer1_point[0] = point[0] + LAYER_OFFSET_X;
        layer1_point[1] = point[1] + LAYER_OFFSET_Y;
    }

    struct simplectic2_points layer1, layer2;

    simplectic2_points(layer1_point, &layer1);
    simplectic2_points(layer2_point, &layer2);

    expand2to3(&layer1.elts[0], layer_int, z_offset, &out->elts[0]);
    expand2to3(&layer1.elts[1], layer_int, z_offset, &out->elts[1]);
    expand2to3(&layer1.elts[2], layer_int, z_offset, &out->elts[2]);

    expand2to3(&layer2.elts[0] + 1, layer_int, z_offset - SIMPLEX_SIZE, &out->elts[3]);
    expand2to3(&layer2.elts[1] + 1, layer_int, z_offset - SIMPLEX_SIZE, &out->elts[4]);
    expand2to3(&layer2.elts[2] + 1, layer_int, z_offset - SIMPLEX_SIZE, &out->elts[5]);
}

static void expand3to4(const simplectic_point3 *p, int64_t cell, double offset, simplectic_point4 *out) {
    out->x_cell = p->x_cell;
    out->x_offset = p->x_offset;

    out->y_cell = p->y_cell;
    out->y_offset = p->y_offset;

    out->z_cell = p->z_cell;
    out->z_offset = p->z_offset;

    out->w_cell = cell;
    out->w_offset = offset;
}

FIXED(simplectic4_points, simplectic_point4, 12);

static void simplectic4_points(const double point[4], struct simplectic4_points *out) {
    double layer = floor(point[3] * INV_SIMPLEX_SIZE),
           w_offset = point[3] - layer * SIMPLEX_SIZE;
    int64_t layer_int = (int64_t)layer;
    double layer1_point[3], layer2_point[3];

    if (layer_int % 2 == 0) {
        layer1_point[0] = point[0];
        layer1_point[1] = point[1];
        layer1_point[2] = point[2];
        layer2_point[0] = point[0] + LAYER_OFFSET_X;
        layer2_point[1] = point[1] + LAYER_OFFSET_Y;
        layer2_point[2] = point[2] + LAYER_OFFSET_Z;
    } else {
        layer2_point[0] = point[0];
        layer2_point[1] = point[1];
        layer2_point[2] = point[2];
        layer1_point[0] = point[0] + LAYER_OFFSET_X;
        layer1_point[1] = point[1] + LAYER_OFFSET_Y;
        layer1_point[2] = point[2] + LAYER_OFFSET_Z;
        layer1_point[2] = point[2] + LAYER_OFFSET_Z;
    }

    struct simplectic3_points layer1, layer2;

    simplectic3_points(layer1_point, &layer1);
    simplectic3_points(layer2_point, &layer2);

    expand3to4(&layer1.elts[0], layer_int, w_offset, &out->elts[0]);
    expand3to4(&layer1.elts[1], layer_int, w_offset, &out->elts[1]);
    expand3to4(&layer1.elts[2], layer_int, w_offset, &out->elts[2]);
    expand3to4(&layer1.elts[3], layer_int, w_offset, &out->elts[3]);
    expand3to4(&layer1.elts[4], layer_int, w_offset, &out->elts[4]);
    expand3to4(&layer1.elts[5], layer_int, w_offset, &out->elts[5]);

    expand3to4(&layer2.elts[0] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[6]);
    expand3to4(&layer2.elts[1] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[7]);
    expand3to4(&layer2.elts[2] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[8]);
    expand3to4(&layer2.elts[3] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[9]);
    expand3to4(&layer2.elts[4] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[10]);
    expand3to4(&layer2.elts[5] + 1, layer_int, w_offset - SIMPLEX_SIZE, &out->elts[11]);
}

#define SMOD(a, b) (a < 0 ? b - (labs(a) % b) : a % b)

static uint8_t seed_get1(const simplectic_seed *sd, uint8_t x) {
    return sd->table[x];
}

static uint8_t seed_get2(const simplectic_seed *sd, uint8_t x, uint8_t y) {
    return sd->table[y + seed_get1(sd, x)];
}

static uint8_t seed_get3(const simplectic_seed *sd, uint8_t x, uint8_t y, uint8_t z) {
    return sd->table[z + seed_get2(sd, x, y)];
}

static uint8_t seed_get4(const simplectic_seed *sd, uint8_t x, uint8_t y, uint8_t z, uint8_t w) {
    return sd->table[w + seed_get3(sd, x, y, z)];
}

#undef SMOD

static void gradient_get2(int64_t idx, double p[const 2]) {
    double diag = 0.70710678118;
    switch (idx % 8) {
        case 0:
            p[0] = diag; p[1] = diag;
            break;
        case 1:
            p[0] = diag; p[1] = -diag;
            break;
        case 2:
            p[0] = -diag; p[1] = diag;
            break;
        case 3:
            p[0] = -diag; p[1] = -diag;
            break;
        case 4:
            p[0] = 1; p[1] = 0;
            break;
        case 5:
            p[0] = -1; p[1] = 0;
            break;
        case 6:
            p[0] = 0; p[1] = 1;
            break;
        case 7:
            p[0] = 0; p[1] = -1;
            break;
    }
}

static void gradient_get3(int64_t idx, double p[const 3]) {
    double diag = 0.70710678118;
    switch (idx % 12) {
        case 0:
            p[0] = diag; p[1] = diag; p[2] = 0;
            break;
        case 1:
            p[0] = diag; p[1] = -diag; p[2] = 0;
            break;
        case 2:
            p[0] = -diag; p[1] = diag; p[2] = 0;
            break;
        case 3:
            p[0] = -diag; p[1] = -diag; p[2] = 0;
            break;
        case 4:
            p[0] = diag; p[1] = 0; p[2] = diag;
            break;
        case 5:
            p[0] = diag; p[1] = 0; p[2] = -diag;
            break;
        case 6:
            p[0] = -diag; p[1] = 0; p[2] = diag;
            break;
        case 7:
            p[0] = -diag; p[1] = 0; p[2] = -diag;
            break;
        case 8:
            p[0] = 0; p[1] = diag; p[2] = diag;
            break;
        case 9:
            p[0] = 0; p[1] = diag; p[2] = -diag;
            break;
        case 10:
            p[0] = 0; p[1] = -diag; p[2] = diag;
            break;
        case 11:
            p[0] = 0; p[1] = -diag; p[2] = -diag;
            break;
    }
}

static void gradient_get4(int64_t idx, double p[const 4]) {
    double diag = 0.57735026919;
    switch (idx % 32) {
        case 0:
            p[0] = diag; p[1] = diag; p[2] = diag; p[3] = 0;
            break;
        case 1:
            p[0] = diag; p[1] = -diag; p[2] = diag; p[3] = 0;
            break;
        case 2:
            p[0] = -diag; p[1] = diag; p[2] = diag; p[3] = 0;
            break;
        case 3:
            p[0] = -diag; p[1] = -diag; p[2] = diag; p[3] = 0;
            break;
        case 4:
            p[0] = diag; p[1] = diag; p[2] = -diag; p[3] = 0;
            break;
        case 5:
            p[0] = diag; p[1] = -diag; p[2] = -diag; p[3] = 0;
            break;
        case 6:
            p[0] = -diag; p[1] = diag; p[2] = -diag; p[3] = 0;
            break;
        case 7:
            p[0] = -diag; p[1] = -diag; p[2] = -diag; p[3] = 0;
            break;
        case 8:
            p[0] = diag; p[1] = diag; p[2] = 0; p[3] = diag;
            break;
        case 9:
            p[0] = diag; p[1] = -diag; p[2] = 0; p[3] = diag;
            break;
        case 10:
            p[0] = -diag; p[1] = diag; p[2] = 0; p[3] = diag;
            break;
        case 11:
            p[0] = -diag; p[1] = -diag; p[2] = 0; p[3] = diag;
            break;
        case 12:
            p[0] = diag; p[1] = diag; p[2] = 0; p[3] = -diag;
            break;
        case 13:
            p[0] = diag; p[1] = -diag; p[2] = 0; p[3] = -diag;
            break;
        case 14:
            p[0] = -diag; p[1] = diag; p[2] = 0; p[3] = -diag;
            break;
        case 15:
            p[0] = -diag; p[1] = -diag; p[2] = 0; p[3] = -diag;
            break;
        case 16:
            p[0] = diag; p[1] = 0; p[2] = diag; p[3] = diag;
            break;
        case 17:
            p[0] = diag; p[1] = 0; p[2] = -diag; p[3] = diag;
            break;
        case 18:
            p[0] = -diag; p[1] = 0; p[2] = diag; p[3] = diag;
            break;
        case 19:
            p[0] = -diag; p[1] = 0; p[2] = -diag; p[3] = diag;
            break;
        case 20:
            p[0] = diag; p[1] = 0; p[2] = diag; p[3] = -diag;
            break;
        case 21:
            p[0] = diag; p[1] = 0; p[2] = -diag; p[3] = -diag;
            break;
        case 22:
            p[0] = -diag; p[1] = 0; p[2] = diag; p[3] = -diag;
            break;
        case 23:
            p[0] = -diag; p[1] = 0; p[2] = -diag; p[3] = -diag;
            break;
        case 24:
            p[0] = 0; p[1] = diag; p[2] = diag; p[3] = diag;
            break;
        case 25:
            p[0] = 0; p[1] = diag; p[2] = -diag; p[3] = diag;
            break;
        case 26:
            p[0] = 0; p[1] = -diag; p[2] = diag; p[3] = diag;
            break;
        case 27:
            p[0] = 0; p[1] = -diag; p[2] = -diag; p[3] = diag;
            break;
        case 28:
            p[0] = 0; p[1] = diag; p[2] = diag; p[3] = -diag;
            break;
        case 29:
            p[0] = 0; p[1] = diag; p[2] = -diag; p[3] = -diag;
            break;
        case 30:
            p[0] = 0; p[1] = -diag; p[2] = diag; p[3] = -diag;
            break;
        case 31:
            p[0] = 0; p[1] = -diag; p[2] = -diag; p[3] = -diag;
            break;
    }
}

double gradient2(const simplectic_seed *seed, const simplectic_point2 *p) {
    double attn = SIMPLEX_SIZE - p->x_offset * p->x_offset
        - p->y_offset * p->y_offset,
        gg[2], attn2 = attn * attn;

    if (attn > 0) {
        gradient_get2(seed_get2(seed, (uint8_t) p->x_cell, (uint8_t)
                    p->y_cell), gg);
        return attn2 * attn2 * (p->x_offset * gg[0] + p->y_offset * gg[1]);
    } else {
        return 0;
    }
}

double simplectic2(const double point[const 2], const simplectic_seed *seed) {
    struct simplectic2_points ps;
    simplectic2_points(point, &ps);

    return
        (gradient2(seed, &ps.elts[0])
         + gradient2(seed, &ps.elts[1])
         + gradient2(seed, &ps.elts[2])) * NORM2;
}

double gradient3(const simplectic_seed *seed, const simplectic_point3 *p) {
    double attn = SIMPLEX_SIZE - p->x_offset * p->x_offset
        - p->y_offset * p->y_offset
        - p->z_offset * p->z_offset,
        gg[3], attn2 = attn * attn;

    if (attn > 0) {
        gradient_get3(seed_get3(seed, (uint8_t) p->x_cell, (uint8_t)
                    p->y_cell, (uint8_t) p->z_cell), gg);
        return attn2 * attn2 * (p->x_offset * gg[0] + p->y_offset * gg[1] + p->z_offset * gg[2]);
    } else {
        return 0;
    }
}

double simplectic3(const double point[const 3], const simplectic_seed *seed) {
    struct simplectic3_points ps;
    simplectic3_points(point, &ps);

    return
        (gradient3(seed, &ps.elts[0])
         + gradient3(seed, &ps.elts[1])
         + gradient3(seed, &ps.elts[2])
         + gradient3(seed, &ps.elts[3])
         + gradient3(seed, &ps.elts[4])
         + gradient3(seed, &ps.elts[5])
         ) * NORM3;
}

double gradient4(const simplectic_seed *seed, const simplectic_point4 *p) {
    double attn = SIMPLEX_SIZE - p->x_offset * p->x_offset
        - p->y_offset * p->y_offset
        - p->z_offset * p->z_offset
        - p->w_offset * p->w_offset,
        gg[4], attn2 = attn * attn;

    if (attn > 0) {
        gradient_get4(seed_get4(seed, (uint8_t) p->x_cell, (uint8_t)
                    p->y_cell, (uint8_t) p->z_cell, (uint8_t) p->w_cell), gg);
        return attn2 * attn2 * (p->x_offset * gg[0] + p->y_offset * gg[1] + p->z_offset * gg[2] + p->w_offset * gg[3]);
    } else {
        return 0;
    }
}

double simplectic4(const double point[const 4], const simplectic_seed *seed) {
    struct simplectic4_points ps;
    simplectic4_points(point, &ps);

    return
        (gradient4(seed, &ps.elts[0])
         + gradient4(seed, &ps.elts[1])
         + gradient4(seed, &ps.elts[2])
         + gradient4(seed, &ps.elts[3])
         + gradient4(seed, &ps.elts[4])
         + gradient4(seed, &ps.elts[5])
         + gradient4(seed, &ps.elts[6])
         + gradient4(seed, &ps.elts[7])
         + gradient4(seed, &ps.elts[8])
         + gradient4(seed, &ps.elts[9])
         + gradient4(seed, &ps.elts[10])
         + gradient4(seed, &ps.elts[11])
         ) * NORM4;
}
