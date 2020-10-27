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
  int r;
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

typedef struct {
  SDL_Renderer *renderer;
  SDL_Window *window;
  SDL_Rect container;
  SDL_Rect water;
  double timescale;
  double xscale;
  double yscale;
  short int pixels[WIN_HEIGHT][WIN_HEIGHT];
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
  /* need to make program read these from config */
  (*s)->r = (int) (a->xscale * 2e-2);
  printf("radius:%d\n",(*s)->r);
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
  double *x;
  x = s->x;
  double xi,yi,ci;
  for (int i = 0; i < n; ++i) {
    fread(&xi, sizeof(xi), 1, fp);
    fread(&yi, sizeof(yi), 1, fp);
    fread(&ci, sizeof(ci), 1, fp);
    *(x++) = xi;
    *(x++) = yi;
  }
  return 1;
}

void _put_in_pix(App *a, int x, int y, int flag) {
  if (x < (WIN_HEIGHT - XOFFSET) && x > (XOFFSET) &&
      y < (WIN_HEIGHT - YOFFSET) && y > (YOFFSET)) {
    a->pixels[x][y]+=flag;
  }
}


void draw_circ(App *a,int x, int y, int r){
  int r2, xy, p1, p2, rl, rl2;
  r2 = r * r;
  rl = (int)(r*1.25);
  rl2 = rl*rl;
  _put_in_pix(a,x,y,2);
  for (p2 = 1; p2 < r; p2++){
    _put_in_pix(a,x,y+p2,2);
    _put_in_pix(a,x,y-p2,2);
  }
  for (p1 = 1; p1 < r; p1++){
    _put_in_pix(a,x+p1,y,2);
    _put_in_pix(a,x-p1,y,2);
  }
  for (p1 = 1; p1 < rl; p1++) {
    for (p2 = 1; p2 < rl; p2++) {
      xy = (p1 * p1 + p2 * p2);
      if (xy < r2) {
        _put_in_pix(a, x + p1, y + p2,2);
        _put_in_pix(a, x + p1, y - p2,2);
        _put_in_pix(a, x - p1, y - p2,2);
        _put_in_pix(a, x - p1, y + p2,2);
      }else if (xy < rl2) {
        _put_in_pix(a, x + p1, y + p2,1);
        _put_in_pix(a, x + p1, y - p2,1);
        _put_in_pix(a, x - p1, y - p2,1);
        _put_in_pix(a, x - p1, y + p2,1);
      }
    }
  }
}

/* int get_pix(App *a, int x, int y){ */
/*   if (x>=0 && y>=0 && x<=WIN_WIDTH && y<WIN_HEIGHT){ */
/*     return a->pixels[x][y]; */
/*   }else{ */
/*     return 0; */
/*   } */
/* } */

/* int _get_count(App *a, int x, int y, int r) { */
/*   int count = 0; */
/*   int i, j; */
/*   int r2, r2l, rl, xy; */
/*   r2 = r * r; */
/*   rl = (int)(r * 1.25); */
/*   r2l = (int)(r2 * 1.5); */
/*   for (i = 0; i < rl; i++) { */
/*     for (j = 0; j < rl; j++) { */
/*       xy = (i * i + j * j); */
/*       if (xy<r2){ */
/* 	count += get_pix(a, x+i, y+i)*2; */
/* 	count += get_pix(a, x+i, y-i)*2; */
/* 	count += get_pix(a, x-i, y-i)*2; */
/* 	count += get_pix(a, x-i, y+i)*2; */
/*       }else if(xy<r2l){ */
/* 	count += get_pix(a, x+i, y+i); */
/* 	count += get_pix(a, x+i, y-i); */
/* 	count += get_pix(a, x-i, y-i); */
/* 	count += get_pix(a, x-i, y+i); */
/*       } */
/*     } */
/*   } */
/*   return count; */
/* } */

void clean_pixels(App *a){
  int p1, p2;
  for(p1=XOFFSET;p1<WIN_WIDTH-XOFFSET;p1++){
    for(p2=YOFFSET;p2<WIN_HEIGHT-YOFFSET;p2++){
      a->pixels[p1][p2]=0;
    }
  }
}


void draw_raster_pix(App *a, int r){
  int p1, p2;
  for(p1=XOFFSET;p1<WIN_WIDTH-XOFFSET;p1++){
    for(p2=YOFFSET;p2<WIN_HEIGHT-YOFFSET;p2++){
      if (a->pixels[p1][p2]>1){
	SDL_RenderDrawPoint(a->renderer, p1, p2);
      }
    }
  }
}

void draw(struct sim_state *s, App *a){
  SDL_SetRenderDrawColor(a->renderer, 0, 0, 0, 255);
  SDL_RenderClear(a->renderer);
  
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 96, 255);
  SDL_RenderDrawRect(a->renderer, &a->container);
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 255, 255);

  int i, n;
  int x, y, r;
  n = s->n;
  r = s->r;
  clean_pixels(a);
  for (i = 0; i < n; i++) {
    x = XOFFSET + (int)(s->x[2 * i] * a->xscale);
    y = WIN_HEIGHT - YOFFSET - (int)(s->x[2 * i + 1] * a->yscale);
    draw_circ(a, x, y, r);
    }
  draw_raster_pix(a,r);
  SDL_RenderPresent(a->renderer);
}

int main(int argc, char *argv[])
{
  int fps;
  if (argc<2 || argv[1][0]=='-'){
    printf("Usage: %s OUTPUT_FILE [FPS]\n\n"
	   "OUTPUT_FILE\tOutput file generated from simulation"
	   " program to visualize.\n"
	   "FPS (30) \tframes per second for display\n",argv[0]);
    return EXIT_FAILURE;
  }else if (argc<3){
    fps = 30;
  }else{
    fps = atoi(argv[2]);
  }
  struct sim_state *sims;
  App mainapp;
  
  FILE *fp;
  fp = fopen(argv[1],"r");
  printf("reading %s\n",argv[1]);
  init_everything(fp,&sims,&mainapp);
  while (next_step(sims,fp)){
    draw(sims,&mainapp);
    handle_input();
    SDL_Delay(1000/fps);
  }
  destroy_everything(sims,&mainapp);
  return EXIT_SUCCESS;
}
