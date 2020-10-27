#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>
#include <libconfig.h>

typedef struct sim_param_t {
  char *fname; /*  File name */
  int nframes; /*  Number of frames */
  int npframe; /*  Steps per frame */
  double h;     /*  Particle size */
  double dt;    /*  Time step */
  double rho0;  /*  Reference density */
  double k;     /*  Bulk modulus */
  double mu;    /*  Viscosity */
  double g;     /*  Gravity strength */
} sim_param_t;
int get_params(int argc, char **argv, sim_param_t *params);

typedef struct sim_state_t {
  int n;
  double mass;
  double *restrict rho;
  double *restrict x;
  double *restrict vh;
  double *restrict v;
  double *restrict a;
} sim_state_t;
sim_state_t *alloc_state(int n){
  sim_state_t *ss;
  ss = malloc(sizeof(sim_state_t));
  ss->n = n;
  ss->rho = malloc(sizeof(double)*n);
  ss->x = malloc(sizeof(double)*2*n);
  ss->v = malloc(sizeof(double)*2*n);
  ss->vh = malloc(sizeof(double)*2*n);
  ss->a = malloc(sizeof(double)*2*n);
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
  double *restrict rho = s->rho;
  const double *restrict x = s->x;
  double h = params->h;
  double h2 = h * h;
  double h8 = (h2 * h2) * (h2 * h2);
  double C = 4 * s->mass / M_PI / h8;

  memset(rho, 0, n * sizeof(double));
  for (int i = 0; i < n; ++i) {
    rho[i] += 4 * s->mass / M_PI / h2;
    for (int j = i + 1; j < n; ++j) {
      double dx = x[2 * i + 0] - x[2 * j + 0];
      double dy = x[2 * i + 1] - x[2 * j + 1];
      double r2 = dx * dx + dy * dy;
      double z = h2 - r2;
      if (z > 0) {
        double rho_ij = C * z * z * z;
        rho[i] += rho_ij;
        rho[j] += rho_ij;
      }
    }
  }
}

void compute_accel(sim_state_t *state, sim_param_t *params) {
  // Unpack basic parameters
  const double h = params->h;
  const double rho0 = params->rho0;
  const double k = params->k;
  const double mu = params->mu;
  const double g = params->g;
  const double mass = state->mass;
  const double h2 = h * h;
  // Unpack system state
  const double *restrict rho = state->rho;
  const double *restrict x = state->x;
  const double *restrict v = state->v;
  double *restrict a = state->a;
  int n = state->n;

  // Compute density and color
  compute_density(state, params);
  // Start with gravity and surface forces
  for (int i = 0; i < n; ++i) {
    a[2 * i + 0] = 0;
    a[2 * i + 1] = -g;
  }
  // Constants for interaction term
  double C0 = mass / M_PI / ((h2) * (h2));
  double Cp = 15 * k;
  double Cv = -40 * mu;
  // Now compute interaction forces
  for (int i = 0; i < n; ++i) {
    const double rhoi = rho[i];
    for (int j = i + 1; j < n; ++j) {
      double dx = x[2 * i + 0] - x[2 * j + 0];
      double dy = x[2 * i + 1] - x[2 * j + 1];
      double r2 = dx * dx + dy * dy;
      if (r2 < h2) {
        const double rhoj = rho[j];
        double q = sqrt(r2) / h;
        double u = 1 - q;
        double w0 = C0 * u / rhoi / rhoj;
        double wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
        double wv = w0 * Cv;
        double dvx = v[2 * i + 0] - v[2 * j + 0];
        double dvy = v[2 * i + 1] - v[2 * j + 1];
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
  const double *restrict a = s->a;
  double *restrict vh = s->vh;
  double *restrict v = s->v;
  double *restrict x = s->x;
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
  const double *restrict a = s->a;
  double *restrict vh = s->vh;
  double *restrict v = s->v;
  double *restrict x = s->x;
  int n = s->n;
  for (int i = 0; i < 2 * n; ++i)
    vh[i] = v[i] + a[i] * dt / 2;
  for (int i = 0; i < 2 * n; ++i)
    v[i] += a[i] * dt;
  for (int i = 0; i < 2 * n; ++i)
    x[i] += vh[i] * dt;
  reflect_bc(s);
}

static void damp_reflect(int which, double barrier, double *x, double *v,
                         double *vh) {
  // Coefficient of resitiution
  const double DAMP = 0.75;
  // Ignore degenerate cases
  if (v[which] == 0)
    return;
  // Scale back the distance traveled based on time from collision
  double tbounce = (x[which] - barrier) / v[which];
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
  const double XMIN = 0.0;
  const double XMAX = 1.0;
  const double YMIN = 0.0;
  const double YMAX = 1.0;
  double *restrict vh = s->vh;
  double *restrict v = s->v;
  double *restrict x = s->x;
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

typedef int (*domain_fun_t)(double, double);
int box_indicator(double x, double y) { return (x < 0.5) && (y < 0.5); }
int circ_indicator(double x, double y) {
  double dx = (x - 0.5);
  double dy = (y - 0.3);
  double r2 = dx * dx + dy * dy;
  return (r2 < 0.25 * 0.25);
}

sim_state_t *place_particles(sim_param_t *param, domain_fun_t indicatef) {
  double h = param->h;
  double hh = h / 1.3;
  // Count mesh points that fall in indicated region.
  int count = 0;
  for (double x = 0; x < 1; x += hh)
    for (double y = 0; y < 1; y += hh)
      count += indicatef(x, y);
  // Populate the particle data structure
  sim_state_t *s = alloc_state(count);
  int p = 0;
  for (double x = 0; x < 1; x += hh) {
    for (double y = 0; y < 1; y += hh) {
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
  double rho0 = param->rho0;
  double rho2s = 0;
  double rhos = 0;
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
    double xi = s->x[2 * i + 0];
    double yi = s->x[2 * i + 1];
    assert(xi >= 0 || xi <= 1);
    assert(yi >= 0 || yi <= 1);
  }
}

void write_header(FILE *fp, int n);
void write_frame_data(FILE *fp, int n, double *x, int *c);

int main(int argc, char **argv) {
  sim_param_t params;
  if (get_params(argc, argv, &params) != 0)
    exit(-1);
  sim_state_t *state = init_particles(&params);
  FILE *fp = fopen(params.fname, "wb");
  int nframes = params.nframes;
  int npframe = params.npframe;
  double dt = params.dt;
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
      printf("%d/%d-%d/%d\r",frame+1,nframes,i+1,npframe);
      fflush(stdout);
    }
    write_frame_data(fp, n, state->x, NULL);
  }
  printf("\nOutput written to:%s\n",params.fname);
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

static void read_conf_file(sim_param_t *params, char* confFile){
  config_t conf,*c;
  const char *ofile;
  c = &conf;
  config_init(c);
  config_read_file(c,confFile);
  config_lookup_string(c,"simpara.output",&ofile);
  strcpy(params->fname = malloc(sizeof(char)*(strlen(ofile)+1)),ofile);
  config_lookup_int(c,"simpara.num_frames",&params->nframes);
  config_lookup_int(c,"simpara.num_per_frames",&params->npframe);
  config_lookup_float(c,"simpara.particle.size",&params->h);
  config_lookup_float(c,"simpara.timestep",&params->dt);
  config_lookup_float(c,"simpara.constants.rho",&params->rho0);
  config_lookup_float(c,"simpara.constants.k",&params->k);
  config_lookup_float(c,"simpara.constants.mu",&params->mu);
  config_lookup_float(c,"simpara.constants.g",&params->g);
  config_destroy(c);
}

static void print_usage() {
  sim_param_t param;
  default_params(&param);
  fprintf(stderr,
          "nbody\n"
          "\t-h: print this message\n"
          "\t-c: config File name\n"
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
  const char *optstring = "ho:F:f:t:s:d:k:v:g:c:";
  int c;
#define get_int_arg(c, field)                                                  \
  case c:                                                                      \
    params->field = atoi(optarg);                                              \
    break
#define get_flt_arg(c, field)                                                  \
  case c:                                                                      \
    params->field = (double)atof(optarg);                                       \
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
    case 'c':
      read_conf_file(params, optarg);
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

void write_frame_data(FILE *fp, int n, double *x, int *c) {
  for (int i = 0; i < n; ++i) {
    double xi = *(x++);
    double yi = *(x++);
    fwrite(&xi, sizeof(xi), 1, fp);
    fwrite(&yi, sizeof(yi), 1, fp);
    double ci0 = c ? *c++ : 0;
    double ci = ci0;
    fwrite(&ci, sizeof(ci), 1, fp);
  }
}
