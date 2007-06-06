/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/render.h>
#include <libgtext/graphics.h>
#include <libgtext/element.h>
#include <libgtext/feature_index.h>
#include <libgtext/track.h>
#include <libgtext/line.h>
#include <libgtext/block.h>

typedef struct {
  Range range;
  double dx, dy;
  unsigned int width;
} RenderInfo;

struct Render {
  Diagram* dia;
  Config* cfg;
  Graphics* g;
  RenderInfo info;
};

int print_tracks(void* key, void* value, void* data, Env* env)
{
  printf("track %s: \n", (char*) key);
  print_track((Track*) value);
  return 0;
}

unsigned int render_calculate_height(Render *r, Env* env)
{
  assert(r && env);
  unsigned int height =  diagram_get_total_lines(r->dia, env);
  /* obtain line height and spacer from configuration settings */
  unsigned int line_height = ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_height",
                                                            env)) +
                             ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_vspace",
                                                            env));
  /* get total height of all lines */
  height *= line_height;
  /* add track spacers */
  height += diagram_get_number_of_tracks(r->dia)
              * config_get_num(r->cfg, "format","track_vspace",env);
  /* add header space and footer */
  height += 70 + 40;
  return height;
}

Render* render_new(Diagram *d, Config *cfg, Env *env)
{
  assert(d && cfg && env);
  Render *r = env_ma_malloc(env, sizeof (Render));
  r->dia = d;
  r->cfg = cfg;
  r->info.range = diagram_get_range(d);
  return r;
}

void render_delete(Render *r, Env *env)
{
  assert(r && env);
  env_ma_free(r, env);
}

void render_line(Render *r, Line *line, Env *env)
{
  assert(r && line && env);
  int i;
  Array *blocks = line_get_blocks(line);
  double margins = config_get_num(r->cfg, "format", "margins", env);
  for (i=0; i<array_size(blocks); i++)
  {
    int j;
    Block *block = *(Block**) array_get(blocks, i);
    Array *elems = block_get_elements(block);
    for (j=0;j<array_size(elems); j++)
    {
      Element *elem = *(Element**) array_get(elems, j);
      Range elem_range = element_get_range(elem);
      double elem_start, elem_end, elem_width;

      if (config_get_verbose(r->cfg))
        printf("processing element from %lu to %lu\n",
               elem_range.start,
               elem_range.end);

      /* obtain scaling factor for target image width */
      double factor = (double) r->info.width / (double) (r->info.range.end -
                                                         r->info.range.start);
      if (config_get_verbose(r->cfg))
        printf("scaling factor is %f\n",factor);

      /* calculate coordinates for drawing */
      unsigned long unscaled_offset_start = (elem_range.start
                                             - r->info.range.start);
      unsigned long unscaled_offset_end   = (elem_range.end
                                             - r->info.range.start);

      /* scale coordinates to target image width, clip to margins if needed */
      if (elem_range.start <= r->info.range.start)
        elem_start = margins;
      else
        elem_start = margins + (double) (unscaled_offset_start * factor);
      if (unscaled_offset_end >= r->info.range.end)
        elem_end = (double) (r->info.width - margins);
      else
        elem_end = margins + (double) (unscaled_offset_end * factor);
      elem_width = elem_end - elem_start;

      if (config_get_verbose(r->cfg))
        printf("drawing element from %.2f to %.2f\n",
               (double) elem_start,
               (double) elem_end);

      /* draw each element according to style set in the config */
      const char* type = genome_feature_type_get_cstr(element_get_type(elem));
      const char* style = config_get_cstr(r->cfg, "feature_styles", type, env);
      if (strcmp(style, "box")==0)
      {
        graphics_draw_box(r->g,
                          elem_start,
                          r->info.dy,
                          elem_width,
                          config_get_num(r->cfg, "format", "bar_height", env),
                          config_get_color(r->cfg, type, env),
                          element_get_arrow_status(elem),
                          config_get_num(r->cfg, "format", "arrow_width", env),
                          config_get_num(r->cfg, "format", "stroke_width", env),
                          config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "caret")==0)
      {
        graphics_draw_caret(r->g,
                          elem_start,
                          r->info.dy,
                          elem_width,
                          config_get_num(r->cfg, "format", "bar_height", env),
                          element_get_arrow_status(elem),
                          config_get_num(r->cfg, "format", "arrow_width", env),
                          config_get_num(r->cfg, "format", "stroke_width", env),
                          config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "dashes")==0)
      {
        graphics_draw_dashes(r->g,
                          elem_start,
                          r->info.dy,
                          elem_width,
                          config_get_num(r->cfg, "format", "bar_height", env),
                          element_get_arrow_status(elem),
                          config_get_num(r->cfg, "format", "arrow_width", env),
                          config_get_num(r->cfg, "format", "stroke_width", env),
                          config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "line")==0)
      {
        graphics_draw_horizontal_line(r->g,
                          elem_start,
                          r->info.dy,
                          elem_width);
      }
      else
      {
         graphics_draw_box(r->g,
                          elem_start,
                          r->info.dy,
                          elem_width,
                          config_get_num(r->cfg, "format", "bar_height", env),
                          config_get_color(r->cfg, type, env),
                          element_get_arrow_status(elem),
                          config_get_num(r->cfg, "format", "arrow_width", env),
                          config_get_num(r->cfg, "format", "stroke_width", env),
                          config_get_color(r->cfg, "stroke", env));
      }
    }
  }
  r->info.dy += config_get_num(r->cfg, "format", "bar_height", env) +
                config_get_num(r->cfg, "format", "bar_vspace", env);
}

int render_track(void *key, void* value, void *data, Env *env)
{
  assert(value && key && data && env);
  Render* r = (Render*) data;
  Track* track = (Track*) value;
  Array* lines = track_get_lines(track);
  int i;

  if (config_get_verbose(r->cfg))
    printf("processing track %s\n", (const char*) key);
  
  /* draw track title */
  graphics_draw_text(r->g, r->info.dx, r->info.dy, (const char*) key);
  r->info.dy += graphics_get_text_height(r->g) + 5;
  
  /* render each line */
  for (i=0; i<array_size(lines); i++)
  {
    render_line(r, *(Line**) array_get(lines, i), env);
  }
  
  /* put track spacer after track */
  r->info.dy += config_get_num(r->cfg, "format", "track_vspace", env);
  return 0;
}

void render_to_png(Render *r, char *fn, unsigned int width, Env *env)
{
  assert(r && fn && env && width > 0);
  unsigned int height = render_calculate_height(r, env);
	char str[32];
  
  /* set initial margins, header, target width */
  r->info.dx = config_get_num(r->cfg, "format", "margins", env);
  r->info.dy = 70;
  r->info.width = width;
  
  /* create new Graphics backend */
  r->g = graphics_new_png(fn, width, height, env);
  sprintf(str, "%lu", r->info.range.start);
	graphics_draw_text_left(r->g, config_get_num(r->cfg, "format", "margins", env), 20,str) ;
	sprintf(str, "%lu", r->info.range.end);
	graphics_draw_text_right(r->g, width-config_get_num(r->cfg, "format", "margins", env), 20,str) ;
	
	/* draw scale and location info */
	graphics_draw_scale(r->g,
	                    config_get_num(r->cfg, "format", "margins", env),
											30,
											width-2*config_get_num(r->cfg, "format", "margins", env),
											config_get_color(r->cfg, "stroke", env),
											Both,
                      config_get_num(r->cfg, "format", "stroke_width", env),
											config_get_num(r->cfg, "format", "scale_arrow_height", env),
											config_get_num(r->cfg, "format", "scale_arrow_width", env));
	
  /* process (render) each track */
  Hashtable *tracks = diagram_get_tracks(r->dia);
  hashtable_foreach(tracks, render_track, r, env);
  
  /* write out result file */
  graphics_save(r->g);
  graphics_delete(r->g, env);
}
