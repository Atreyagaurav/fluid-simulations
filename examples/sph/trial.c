#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>

typedef struct sim_param_t {
  char *fname; /*  File name */
  int nframes; /*  Number of frames */
  int npframe; /*  Steps per frame */
  float h;     /*  Particle size */
  float dt;    /*  Time step */
  float rho0;  /*  Reference density */
  float k;     /*  Bulk modulus */
  float mu;    /*  Viscosity */
  float g;     /*  Gravity strength */
} sim_param_t;
int get_params(int argc, char **argv, sim_param_t *params);

typedef struct sim_state_t {
  int n;
  float mass;
  float *restrict rho;
  float *restrict x;
  float *restrict vh;
  float *restrict v;
  float *restrict a;
} sim_state_t;
sim_state_t *alloc_state(int n){
  sim_state_t *ss;
  ss = malloc(sizeof(sim_state_t));
  ss->n = n;
  ss->rho = malloc(sizeof(float)*n);
  ss->x = malloc(sizeof(float)*2*n);
  ss->v = malloc(sizeof(float)*2*n);
  ss->vh = malloc(sizeof(float)*2*n);
  ss->a = malloc(sizeof(float)*2*n);
  return ss;
};
void free_state(sim_state_t *s){
  free(s->rho);
  free(s->x);
  free(s->v);
  free(s->vh);
  free(s->a);
  free(s);
};

void compute_density(sim_state_t *s, sim_param_t *params) {
  int n = s->n;
  float *restrict rho = s->rho;
  const float *restrict x = s->x;
  float h = params->h;
  float h2 = h * h;
  float h8 = (h2 * h2) * (h2 * h2);
  float C = 4 * s->mass / M_PI / h8;

  memset(rho, 0, n * sizeof(float));
  for (int i = 0; i < n; ++i) {
    rho[i] += 4 * s->mass / M_PI / h2;
    for (int j = i + 1; j < n; ++j) {
      float dx = x[2 * i + 0] - x[2 * j + 0];
      float dy = x[2 * i + 1] - x[2 * j + 1];
      float r2 = dx * dx + dy * dy;
      float z = h2 - r2;
      if (z > 0) {
        float rho_ij = C * z * z * z;
        rho[i] += rho_ij;
        rho[j] += rho_ij;
      }
    }
  }
}

void compute_accel(sim_state_t *state, sim_param_t *params) {
  // Unpack basic parameters
  const float h = params->h;
  const float rho0 = params->rho0;
  const float k = params->k;
  const float mu = params->mu;
  const float g = params->g;
  const float mass = state->mass;
  const float h2 = h * h;
  // Unpack system state
  const float *restrict rho = state->rho;
  const float *restrict x = state->x;
  const float *restrict v = state->v;
  float *restrict a = state->a;
  int n = state->n;

  // Compute density and color
  compute_density(state, params);
  // Start with gravity and surface forces
  for (int i = 0; i < n; ++i) {
    a[2 * i + 0] = 0;
    a[2 * i + 1] = -g;
  }
  // Constants for interaction term
  float C0 = mass / M_PI / ((h2) * (h2));
  float Cp = 15 * k;
  float Cv = -40 * mu;
  // Now compute interaction forces
  for (int i = 0; i < n; ++i) {
    const float rhoi = rho[i];
    for (int j = i + 1; j < n; ++j) {
      float dx = x[2 * i + 0] - x[2 * j + 0];
      float dy = x[2 * i + 1] - x[2 * j + 1];
      float r2 = dx * dx + dy * dy;
      if (r2 < h2) {
        const float rhoj = rho[j];
        float q = sqrt(r2) / h;
        float u = 1 - q;
        float w0 = C0 * u / rhoi / rhoj;
        float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
        float wv = w0 * Cv;
        float dvx = v[2 * i + 0] - v[2 * j + 0];
        float dvy = v[2 * i + 1] - v[2 * j + 1];
        a[2 * i + 0] += (wp * dx + wv * dvx);
        a[2 * i + 1] += (wp * dy + wv * dvy);
        a[2 * j + 0] -= (wp * dx + wv * dvx);
        a[2 * j + 1] -= (wp * dy + wv * dvy);
      }
    }
  }
}
static void reflect_bc(sim_state_t* s);

void leapfrog_step(sim_state_t *s, double dt) {
  const float *restrict a = s->a;
  float *restrict vh = s->vh;
  float *restrict v = s->v;
  float *restrict x = s->x;
  int n = s->n;
  for (int i = 0; i < 2 * n; ++i)
    vh[i] += a[i] * dt;
  for (int i = 0; i < 2 * n; ++i)
    v[i] = vh[i] + a[i] * dt / 2;
  for (int i = 0; i < 2 * n; ++i)
    x[i] += vh[i] * dt;
  reflect_bc(s);
}

void leapfrog_start(sim_state_t *s, double dt) {
  const float *restrict a = s->a;
  float *restrict vh = s->vh;
  float *restrict v = s->v;
  float *restrict x = s->x;
  int n = s->n;
  for (int i = 0; i < 2 * n; ++i)
    vh[i] = v[i] + a[i] * dt / 2;
  for (int i = 0; i < 2 * n; ++i)
    v[i] += a[i] * dt;
  for (int i = 0; i < 2 * n; ++i)
    x[i] += vh[i] * dt;
  reflect_bc(s);
}

static void damp_reflect(int which, float barrier, float *x, float *v,
                         float *vh) {
  // Coefficient of resitiution
  const float DAMP = 0.75;
  // Ignore degenerate cases
  if (v[which] == 0)
    return;
  // Scale back the distance traveled based on time from collision
  float tbounce = (x[which] - barrier) / v[which];
  x[0] -= v[0] * (1 - DAMP) * tbounce;
  x[1] -= v[1] * (1 - DAMP) * tbounce;
  // Reflect the position and velocity
  x[which] = 2 * barrier - x[which];
  v[which] = -v[which];
  vh[which] = -vh[which];
  // Damp the velocities
  v[0] *= DAMP;
  vh[0] *= DAMP;
  v[1] *= DAMP;
  vh[1] *= DAMP;
}

static void reflect_bc(sim_state_t *s) {
  // Boundaries of the computational domain
  const float XMIN = 0.0;
  const float XMAX = 1.0;
  const float YMIN = 0.0;
  const float YMAX = 1.0;
  float *restrict vh = s->vh;
  float *restrict v = s->v;
  float *restrict x = s->x;
  int n = s->n;
  for (int i = 0; i < n; ++i, x += 2, v += 2, vh += 2) {
    if (x[0] < XMIN)
      damp_reflect(0, XMIN, x, v, vh);
    if (x[0] > XMAX)
      damp_reflect(0, XMAX, x, v, vh);
    if (x[1] < YMIN)
      damp_reflect(1, YMIN, x, v, vh);
    if (x[1] > YMAX)
      damp_reflect(1, YMAX, x, v, vh);
  }
}

typedef int (*domain_fun_t)(float, float);
int box_indicator(float x, float y) { return (x < 0.5) && (y < 0.5); }
int circ_indicator(float x, float y) {
  float dx = (x - 0.5);
  float dy = (y - 0.3);
  float r2 = dx * dx + dy * dy;
  return (r2 < 0.25 * 0.25);
}

sim_state_t *place_particles(sim_param_t *param, domain_fun_t indicatef) {
  float h = param->h;
  float hh = h / 1.3;
  // Count mesh points that fall in indicated region.
  int count = 0;
  for (float x = 0; x < 1; x += hh)
    for (float y = 0; y < 1; y += hh)
      count += indicatef(x, y);
  // Populate the particle data structure
  sim_state_t *s = alloc_state(count);
  int p = 0;
  for (float x = 0; x < 1; x += hh) {
    for (float y = 0; y < 1; y += hh) {
      if (indicatef(x, y)) {
        s->x[2 * p + 0] = x;
        s->x[2 * p + 1] = y;
        s->v[2 * p + 0] = 0;
        s->v[2 * p + 1] = 0;
        ++p;
      }
    }
  }
  return s;
}

void normalize_mass(sim_state_t *s, sim_param_t *param) {
  s->mass = 1;
  compute_density(s, param);
  float rho0 = param->rho0;
  float rho2s = 0;
  float rhos = 0;
  for (int i = 0; i < s->n; ++i) {
    rho2s += (s->rho[i]) * (s->rho[i]);
    rhos += s->rho[i];
  }
  s->mass *= (rho0 * rhos / rho2s);
}
sim_state_t *init_particles(sim_param_t *param) {
  sim_state_t *s = place_particles(param, box_indicator);
  normalize_mass(s, param);
  return s;
}

void check_state(sim_state_t *s) {
  for (int i = 0; i < s->n; ++i) {
    float xi = s->x[2 * i + 0];
    float yi = s->x[2 * i + 1];
    assert(xi >= 0 || xi <= 1);
    assert(yi >= 0 || yi <= 1);
  }
}

void write_header(FILE *fp, int n);
void write_frame_data(FILE *fp, int n, float *x, int *c);

int main(int argc, char **argv) {
  sim_param_t params;
  if (get_params(argc, argv, &params) != 0)
    exit(-1);
  sim_state_t *state = init_particles(&params);
  FILE *fp = fopen(params.fname, "wb");
  int nframes = params.nframes;
  int npframe = params.npframe;
  float dt = params.dt;
  int n = state->n;
  write_header(fp, n);
  write_frame_data(fp, n, state->x, NULL);
  compute_accel(state, &params);
  leapfrog_start(state, dt);
  check_state(state);
  for (int frame = 1; frame < nframes; ++frame) {
    for (int i = 0; i < npframe; ++i) {
      compute_accel(state, &params);
      leapfrog_step(state, dt);
      check_state(state);
      printf("%d/%d-%d/%d\r",frame,nframes,i,npframe);
      fflush(stdout);
    }
    write_frame_data(fp, n, state->x, NULL);
  }
  fclose(fp);
  free_state(state);
}

static void default_params(sim_param_t *params) {
  params->fname = "run.out";
  params->nframes = 400;
  params->npframe = 100;
  params->dt = 1e-4;
  params->h = 5e-2;
  params->rho0 = 1000;
  params->k = 1e3;
  params->mu = 0.1;
  params->g = 9.8;
}

static void print_usage() {
  sim_param_t param;
  default_params(&param);
  fprintf(stderr,
          "nbody\n"
          "\t-h: print this message\n"
          "\t-o: output file name (%s)\n"
          "\t-F: number of frames (%d)\n"
          "\t-f: steps per frame (%d)\n"
          "\t-t: time step (%e)\n"
          "\t-s: particle size (%e)\n"
          "\t-d: reference density (%g)\n"
          "\t-k: bulk modulus (%g)\n"
          "\t-v: dynamic viscosity (%g)\n"
          "\t-g: gravitational strength (%g)\n",
          param.fname, param.nframes, param.npframe, param.dt, param.h,
          param.rho0, param.k, param.mu, param.g);
}

int get_params(int argc, char **argv, sim_param_t *params) {
  extern char *optarg;
  const char *optstring = "ho:F:f:t:s:d:k:v:g:";
  int c;
#define get_int_arg(c, field)                                                  \
  case c:                                                                      \
    params->field = atoi(optarg);                                              \
    break
#define get_flt_arg(c, field)                                                  \
  case c:                                                                      \
    params->field = (float)atof(optarg);                                       \
    break
  default_params(params);
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'h':
      print_usage();
      return -1;
    case 'o':
      strcpy(params->fname = malloc(strlen(optarg) + 1), optarg);
      break;
      get_int_arg('F', nframes);
      get_int_arg('f', npframe);
      get_flt_arg('t', dt);
      get_flt_arg('s', h);
      get_flt_arg('d', rho0);
      get_flt_arg('k', k);
      get_flt_arg('v', mu);
      get_flt_arg('g', g);
    default:
      fprintf(stderr, "Unknown option\n");
      return -1;
    }
  }
  return 0;
}


void write_header(FILE *fp, int n) {
  fwrite(&n, sizeof(n), 1, fp);
}

void write_frame_data(FILE *fp, int n, float *x, int *c) {
  for (int i = 0; i < n; ++i) {
    float xi = *(x++);
    float yi = *(x++);
    fwrite(&xi, sizeof(xi), 1, fp);
    fwrite(&yi, sizeof(yi), 1, fp);
    float ci0 = c ? *c++ : 0;
    float ci = ci0;
    fwrite(&ci, sizeof(ci), 1, fp);
  }
}
