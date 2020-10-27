#include <SDL2/SDL_render.h>
#include <SDL2/SDL_video.h>
#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_timer.h>
#include <stdlib.h>
#include <stdint.h>

#define WIN_WIDTH 500
#define WIN_HEIGHT 500
#define XOFFSET 10
#define YOFFSET 10

typedef struct sim_state {
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

typedef struct {
  SDL_Renderer *renderer;
  SDL_Window *window;
  SDL_Rect container;
  SDL_Rect water;
  double timescale;
  double xscale;
  double yscale;
} App;


void init_everything(FILE *dataFilep, struct sim_state **s, App *a){
  /* initialize the simulation state */
  FILE *fp;
  fp = dataFilep;
  int n;
  fread(&n,sizeof(int),1,fp);
  printf("Number of particles: %d\n",n);
  *s = alloc_state(n);

  /* init the graphics window */
  int rendererFlags, windowFlags;
  rendererFlags = SDL_RENDERER_ACCELERATED;
  windowFlags = 0;
  SDL_Init(SDL_INIT_VIDEO);
  a->window = SDL_CreateWindow("SIMULATION", // creates a window
                                     SDL_WINDOWPOS_CENTERED,
                                     SDL_WINDOWPOS_CENTERED, WIN_WIDTH, WIN_HEIGHT, windowFlags);
  SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
  a->renderer = SDL_CreateRenderer(a->window, -1, rendererFlags);

  a->container.x = XOFFSET;
  a->container.y = YOFFSET;
  a->xscale = (WIN_WIDTH-2*XOFFSET) / 1.0;
  a->yscale = (WIN_HEIGHT-2*YOFFSET) / 1.0;
  a->container.w = (WIN_WIDTH-2*XOFFSET);
  a->container.h = (WIN_WIDTH-2*YOFFSET);
}

void destroy_everything(sim_state_t *s, App *a) {
  free_state(s);
  SDL_DestroyRenderer(a->renderer);
  SDL_DestroyWindow(a->window);
}

void handle_input(void) {
  SDL_Event event;

  while (SDL_PollEvent(&event)) {
    switch (event.type) {
    case SDL_QUIT:
      exit(EXIT_SUCCESS);
      break;
    default:
      break;
    }
  }
}

int next_step(struct sim_state *s, FILE *fp){
  if (feof(fp)){
    return 0;
  }
  int n = s->n;
  float *x;
  x = s->x;
  float xi,yi,ci;
  for (int i = 0; i < n; ++i) {
    fread(&xi, sizeof(xi), 1, fp);
    fread(&yi, sizeof(yi), 1, fp);
    fread(&ci, sizeof(ci), 1, fp);
    *(x++) = xi;
    *(x++) = yi;
  }
  return 1;
}

void draw(struct sim_state *s, App *a){
  SDL_SetRenderDrawColor(a->renderer, 0, 0, 0, 255);
  SDL_RenderClear(a->renderer);
  
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 96, 255);
  SDL_RenderDrawRect(a->renderer, &a->container);
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 255, 255);

  int i, n;
  int x, y;
  n = s->n;
  for (i = 0; i < n; i++) {
    x = XOFFSET + (int)(s->x[2 * i] * a->xscale);
    y = WIN_HEIGHT - YOFFSET - (int)(s->x[2 * i + 1] * a->yscale);
    SDL_RenderDrawPoint(a->renderer, x, y);
  }
  SDL_RenderPresent(a->renderer);
}

int main(int argc, char *argv[])
{
  if (argc<2){
    printf("Usage: %s OUTPUT_FILE\n\n"
	   "OUTPUT_FILE\tOutput file with details for simulation visualization.\n",argv[0]);
    return EXIT_FAILURE;
  }
  struct sim_state *sims;
  App mainapp;
  
  const char *name;
  FILE *fp;
  fp = fopen(argv[1],"r");
  printf("reading %s\n",argv[1]);
  init_everything(fp,&sims,&mainapp);
  while (next_step(sims,fp)){
    draw(sims,&mainapp);
    handle_input();
    SDL_Delay(100);
  }
  destroy_everything(sims,&mainapp);
  return EXIT_SUCCESS;
}
