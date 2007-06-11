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

enum ClipType
{
  CLIPPED_RIGHT = 1,
	CLIPPED_LEFT = 2,
	CLIPPED_NONE = 3,
	CLIPPED_BOTH = 4
};

typedef struct
{
  double start, end;
	int clip;
} DrawingRange;

struct Render {
  Diagram* dia;
  Config* cfg;
  Graphics* g;
  Range range;
  double dx, dy, margins;
  unsigned int width;
};

int print_tracks(void* key, void* value, void* data, Env* env)
{
  printf("track %s: \n", (char*) key);
  print_track((Track*) value);
  return 0;
}

/*!
Calculates the final height of the image to be created.
\param r Render object
\param env Pointer to Environment object.
\return Height of the image in pixels.
*/
unsigned int render_calculate_height(Render *r, Env* env)
{
  assert(r && env);
  unsigned int lines = diagram_get_total_lines(r->dia, env);
  unsigned int height;
  /* obtain line height and spacer from configuration settings */
  unsigned int line_height = ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_height",
                                                            15,
                                                            env)) +
                             ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_vspace",
                                                            10,
                                                            env));
  /* get total height of all lines */
  height = lines * line_height;
  /* add track spacers */
  height += diagram_get_number_of_tracks(r->dia)
              * (config_get_num(r->cfg, "format","track_vspace", 20, env));
  /* add space for captions above each line */
  height += lines * 10;
  /* add header space and footer */
  height += 70 + 40;
  return height;
}

/*!
Creates a new Render object.
\param d Diagram to render
\param cfg Pointer to Config object.
\param env Pointer to Environment object.
\return Newly creates Render object.
*/
Render* render_new(Diagram *d, Config *cfg, Env *env)
{
  assert(d && cfg && env);
  Render *r = env_ma_malloc(env, sizeof (Render));
  r->dia = d;
  r->cfg = cfg;
  r->range = diagram_get_range(d);
  r->margins = config_get_num(r->cfg, "format", "margins", 10, env);
  return r;
}

/*!
Deletes a Render object
\param r Render object
\param env Pointer to Environment object.
*/
void render_delete(Render *r, Env *env)
{
  assert(r && env);
  env_ma_free(r, env);
}

/*!
Converts an ordinary Range into a Range scaled by the given factor.
If the range exceeds boundaries, it will
\param r Render object
\param env Pointer to Environment object.
\return Height of the image in pixels.
*/
DrawingRange render_convert_coords(Render *r,
                                   Range node_range,
                                   double factor,
                                   bool ensure_caption,
                                   Env *env)
{
  DrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;

  /* subtract viewing offset to get a 0-based Range */
  long unscaled_start = ((long) node_range.start - (long) r->range.start);
  long unscaled_end   = ((long) node_range.end - (long) r->range.start); 

  /* scale coordinates to target image width, clip to margins if needed */
  if (ensure_caption && (long) node_range.start < (long) r->range.start )
  {
    converted_range.start = r->margins;
    converted_range.clip = CLIPPED_LEFT;
  }
  else
    converted_range.start = unscaled_start * factor;
  if (ensure_caption && (long) node_range.end > (long) r->range.end)
  {
    converted_range.end = (double) (r->width - r->margins);
    converted_range.clip = (converted_range.clip == CLIPPED_LEFT ? CLIPPED_BOTH : CLIPPED_RIGHT);
  }
  else
    converted_range.end = unscaled_end * factor;

  return converted_range;
}

void render_line(Render *r, Line *line, Env *env)
{
  assert(r && line && env);
  int i;
  Array *blocks = line_get_blocks(line);
  /* obtain scaling factor for target image width */
  double factor = ((double) r->width
                            -(2*r->margins)) 
                 / (double) (r->range.end 
                            - r->range.start);
  if (config_get_verbose(r->cfg))
     printf("scaling factor is %f\n",factor);
  /* begin drawing block */
  for (i=0; i<array_size(blocks); i++)
  {
    int j;
    Block *block = *(Block**) array_get(blocks, i);
    Array *elems = block_get_elements(block);
    Range block_range = block_get_range(block);
    DrawingRange draw_range;
    const char* caption;
    /* draw block caption */
    draw_range = render_convert_coords(r, block_range, factor, true, env);
    caption = block_get_caption(block);
    if (!caption) caption=""; 
    graphics_draw_text(r->g,
                       draw_range.start,
                       r->dy-graphics_get_text_height(r->g)+3,
                       caption);
    /* draw elements in block */
    for (j=0;j<array_size(elems); j++)
    {
      Element *elem = *(Element**) array_get(elems, j);
      Range elem_range = element_get_range(elem);
      DrawingRange draw_range;
      double elem_start, elem_width;

      if (config_get_verbose(r->cfg))
        printf("processing element from %lu to %lu\n",
               elem_range.start,
               elem_range.end);

      draw_range = render_convert_coords(r, elem_range, factor, false, env);
      elem_start = draw_range.start;
      elem_width = draw_range.end - draw_range.start;

      if (config_get_verbose(r->cfg))
        printf("drawing element from %f to %f\n",
               draw_range.start,
               draw_range.end);

      /* draw each element according to style set in the config */
      const char* type = genome_feature_type_get_cstr(element_get_type(elem));
      const char* style = config_get_cstr(r->cfg,
                                          "feature_styles",
                                          type,
                                          "box",
                                          env);
      if (strcmp(style, "box")==0)
      {
        graphics_draw_box(r->g,
                       elem_start,
                       r->dy,
                       elem_width,
                       config_get_num(r->cfg, "format", "bar_height", 15, env),
                       config_get_color(r->cfg, type, env),
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "caret")==0)
      {
        graphics_draw_caret(r->g,
                       elem_start,
                       r->dy,
                       elem_width,
                       config_get_num(r->cfg, "format", "bar_height", 15, env),
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "dashes")==0)
      {
        graphics_draw_dashes(r->g,
                       elem_start,
                       r->dy,
                       elem_width,
                       config_get_num(r->cfg, "format", "bar_height", 15, env),
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "line")==0)
      {
        graphics_draw_horizontal_line(r->g,
                       elem_start,
                       r->dy,
                       elem_width);
      }
      else
      {
         graphics_draw_box(r->g,
                       elem_start,
                       r->dy,
                       elem_width,
                       config_get_num(r->cfg, "format", "bar_height", 15, env),
                       config_get_color(r->cfg, type, env),
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
    }
  }
  if (i!=array_size(blocks)-1)
  r->dy += config_get_num(r->cfg, "format", "bar_height", 15, env) +
                config_get_num(r->cfg, "format", "bar_vspace", 10, env) +
                graphics_get_text_height(r->g);
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
  graphics_draw_colored_text(r->g,
	                           r->dx,
														 r->dy,
														 config_get_color(r->cfg, "track_title", env),
														 (const char*) key);
  r->dy += graphics_get_text_height(r->g) + 10;

  /* render each line */
  for (i=0; i<array_size(lines); i++)
  {
    render_line(r, *(Line**) array_get(lines, i), env);
  }

  /* put track spacer after track */
  r->dy += config_get_num(r->cfg, "format", "track_vspace", 10, env);
  return 0;
}

void render_to_png(Render *r, char *fn, unsigned int width, Env *env)
{
  assert(r && fn && env && width > 0);
  unsigned int height = render_calculate_height(r, env);
  long i;
  char str[32];

  /* set initial margins, header, target width */
  r->dx = r->margins;
  r->dy = 70;
  r->width = width;

  /* create new Graphics backend */
  r->g = graphics_new_png(fn, width, height, env);
  graphics_set_margins(r->g, r->margins, 0, width, height);
    
  /* draw grid, still to be improved -- ssteinbiss */ 
  for(i = r->range.start - (r->range.start % 1000);
      i<=r->range.end - (r->range.start % 100);
      i+=100)
  {
    long x = i - (long) r->range.start;
    double factor = ((double) r->width
                            -(2*r->margins)) 
                 / (double) (r->range.end 
                            - r->range.start);
    graphics_draw_vertical_line(r->g, x * factor, 0, height);
  }
	
	sprintf(str, "%lu", r->range.start);
	graphics_draw_text_left(r->g, r->margins, 20, str);
	sprintf(str, "%lu", r->range.end);
	graphics_draw_text_right(r->g, width-r->margins, 20, str);

  
	/* draw scale and location info */
	graphics_draw_scale(r->g,
	                    r->margins,
											30,
											width-2*r->margins,
											config_get_color(r->cfg, "stroke", env),
											Right,
                      config_get_num(r->cfg,
                                     "format",
                                     "stroke_width",
                                     1,
                                     env),
											config_get_num(r->cfg,
                                     "format",
                                     "scale_arrow_height",
                                     10,
                                     env),
											config_get_num(r->cfg,
                                     "format",
                                     "scale_arrow_width",
                                     6,
                                     env));
                           
  /* process (render) each track */
  Hashtable *tracks = diagram_get_tracks(r->dia);
  hashtable_foreach(tracks, render_track, r, env);

  /* write out result file */
  graphics_save(r->g);
  graphics_delete(r->g, env);
}
