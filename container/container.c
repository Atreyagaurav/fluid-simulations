#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_timer.h>
#include <libconfig.h>
#include <stdlib.h>

#define CONST_g 9.81

#define WIN_WIDTH 150
#define WIN_HEIGHT 500


struct container{
  double height;
  double fill;
  double diameter;
  double elevation;
  double area;
};

struct pore{
  double diameter;
  double elevation;
  double area;
};

struct sim_state {
  struct container c1;
  struct pore p1;
  int step;
  double delT;
  double timestep;
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

void init_everything(config_t *conf, char *confFile, struct sim_state *s, App *a){
  /* Read the config file for initial conditions */
  config_init(conf);
  config_read_file(conf, confFile);
  config_lookup_float(conf,"container.height",&(s->c1.height));
  config_lookup_float(conf,"container.fill",&(s->c1.fill));
  config_lookup_float(conf,"container.diameter",&(s->c1.diameter));
  config_lookup_float(conf,"container.elevation",&(s->c1.elevation));
  s->c1.area = M_PI_4 * s->c1.diameter * s->c1.diameter;
  config_lookup_float(conf,"pore.diameter",&(s->p1.diameter));
  config_lookup_float(conf,"pore.elevation",&(s->p1.elevation));
  config_lookup_float(conf,"simulation.timestep",&(s->timestep));
  config_lookup_float(conf,"simulation.timescale",&(a->timescale));
  s->p1.area = M_PI_4 * s->p1.diameter * s->p1.diameter;
  s->step = 0;

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

  a->container.x = 10;
  a->container.y = 10;
  a->xscale = (WIN_WIDTH-20) / s->c1.diameter;
  a->yscale = (WIN_HEIGHT-20) / s->c1.height;
  a->container.w = s->c1.diameter * a->xscale;
  a->container.h = s->c1.height * a->yscale;
}

void destroy_everything(config_t *conf) { config_destroy(conf); }

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

int next_step(struct sim_state *s){
  double discharge, height;
  height = s->c1.height * s->c1.fill + s->c1.elevation - s->p1.elevation;
  if (height<=0){
    return 0;
  }
  discharge = s->p1.area * pow(2*CONST_g*height, 0.5);
  s->c1.fill -= discharge * s->delT / s->c1.area / s->c1.height;
  s->step +=1;
  printf("%d:%6.4f\t%6.4f\n",s->step,height,discharge);
  return 1;
}

void draw(struct sim_state *s, App *a){
  a->water = a->container;
  a->water.y += (1 - s->c1.fill)*s->c1.height * a->yscale;
  a->water.h = s->c1.fill * s->c1.height * a->yscale;
  
  SDL_SetRenderDrawColor(a->renderer, 0, 0, 0, 255);
  SDL_RenderClear(a->renderer);
  
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 96, 255);
  SDL_RenderDrawRect(a->renderer, &a->container);
  SDL_SetRenderDrawColor(a->renderer, 96, 96, 255, 255);
  SDL_RenderFillRect(a->renderer, &a->water);
  SDL_RenderPresent(a->renderer);
}



int main(int argc, char *argv[])
{
  config_t simconfig;
  struct sim_state sims;
  App mainapp;
  
  char filename[] = "conf1.conf";
  const char *name;
  init_everything(&simconfig,filename,&sims,&mainapp);
  sims.delT = 0.01;
  config_lookup_string(&simconfig,"name",&name);
  printf("Analysing %s\n",name);

  while (next_step(&sims)){
    draw(&sims,&mainapp);
    handle_input();
    SDL_Delay(sims.delT*1000/mainapp.timescale);
  }
  destroy_everything(&simconfig);
  return 0;
}
