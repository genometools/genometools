/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <math.h>
#include <libgtext/render.h>
#include <libgtext/graphics.h>
#include <libgtext/element.h>
#include <libgtext/feature_index.h>
#include <libgtext/track.h>
#include <libgtext/line.h>
#include <libgtext/block.h>

enum ClipType
{
  CLIPPED_RIGHT,
	CLIPPED_LEFT,
	CLIPPED_NONE,
	CLIPPED_BOTH
};

typedef struct
{
  double start, end;
	enum ClipType clip;
} DrawingRange;

struct Render 
{
  Diagram* dia;
  Config* cfg;
  Graphics* g;
  Range range;
  double y, margins, factor;
  unsigned int width, height;
};

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

double render_convert_point(Render *r, long pos)
{
  return (double) ((r->factor * MAX(0,(pos-(long) r->range.start))) + r->margins);
}

/*!
Converts an base range into a pixel range.
If the range exceeds visibility boundaries, clipping info is set.
\param r Render object
\param node_range Range to convert.
\param env Pointer to Environment object.
\return Height of the image in pixels.
*/
DrawingRange render_convert_coords(Render *r,
                                   Range node_range,
                                   Env *env)
{
  DrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;

  /* scale coordinates to target image width */
	/* first, check if left side has to be clipped */
  if ((long) node_range.start < (long) r->range.start )
	{
	  converted_range.clip = CLIPPED_LEFT;
		converted_range.start = MAX(0, r->margins-5);
	} 
	else
  {
    converted_range.start = render_convert_point(r, node_range.start);
  }
  
	/* then, check right side. */
  if ((long) node_range.end > (long) r->range.end)
	{
    converted_range.clip = (converted_range.clip == CLIPPED_LEFT ?
                                                      CLIPPED_BOTH : 
                                                      CLIPPED_RIGHT);
		converted_range.end = r->width - r->margins+5;
	}
	else
	{
    converted_range.end = render_convert_point(r, node_range.end);
	}
  return converted_range;
}

/*!
Renders a line.
\param r Render object
\param line Line object to render.
\param env Pointer to Environment object.
*/
void render_line(Render *r, Line *line, Env *env)
{
  assert(r && line && env);
  int i;
  Array *blocks = line_get_blocks(line);

  if (config_get_verbose(r->cfg))
     printf("scaling factor is %f\n",r->factor);
     
  /* begin drawing block */
  for (i=0; i<array_size(blocks); i++)
  {
    Dlistelem *delem;
    Block *block = *(Block**) array_get(blocks, i);
    Dlist *elems = block_get_elements(block);
    Range block_range = block_get_range(block);
    DrawingRange draw_range;
    const char* caption;
    
    /* draw block caption */
    draw_range = render_convert_coords(r, block_range, env);
    caption = block_get_caption(block);
    if (!caption) caption="<unnamed>"; 
    graphics_draw_text(r->g,
                       MAX(r->margins, draw_range.start),
                       r->y-graphics_get_text_height(r->g)+3,
                       caption);
											 
		/* DEBUG */
         graphics_draw_box(r->g,
                       draw_range.start,
                       r->y,
                       draw_range.end - draw_range.start,
                       config_get_num(r->cfg, "format", "bar_height", 15, env),
                       config_get_color(r->cfg, "foo", env),
                       NoArrow,
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
		/* DEBUG */
					
    /* draw elements in block */
    for(delem = dlist_first(elems); delem != NULL; delem = dlistelem_next(delem))
    {
      Element *elem = (Element*) dlistelem_get_data(delem);
      Range elem_range = element_get_range(elem);
      DrawingRange draw_range;
      double elem_start, elem_width, bar_height;
      Color grey;
      
      grey.red=grey.green=grey.blue=.8;
      bar_height = config_get_num(r->cfg, "format", "bar_height", 15, env);

      if (config_get_verbose(r->cfg))
        printf("processing element from %lu to %lu\n",
               elem_range.start,
               elem_range.end);

      draw_range = render_convert_coords(r, elem_range, env);
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
                       r->y,
                       elem_width,
                       bar_height,
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
                       r->y,
                       elem_width,
                       bar_height,
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "dashes")==0)
      {
        graphics_draw_dashes(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      else if (strcmp(style, "line")==0)
      {
        graphics_draw_horizontal_line(r->g,
                       elem_start,
                       r->y,
                       elem_width);
      }
      else
      {
         graphics_draw_box(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       config_get_color(r->cfg, type, env),
                       element_get_arrow_status(elem),
                       config_get_num(r->cfg, "format", "arrow_width", 6, env),
                       config_get_num(r->cfg, "format", "stroke_width", 1, env),
                       config_get_color(r->cfg, "stroke", env));
      }
      
      /* draw arrowheads at clipped margins */
      if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g, r->margins-10, r->y+((bar_height-8)/2), grey, Left);
      if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g, r->width-r->margins+10, r->y+((bar_height-8)/2), grey, Right);
    }
  }
  /* do not add line spacing after the last line of a track */
  if (i!=array_size(blocks)-1)
    r->y += config_get_num(r->cfg, "format", "bar_height", 15, env) +
               config_get_num(r->cfg, "format", "bar_vspace", 10, env) +
               graphics_get_text_height(r->g);
}

/*!
Renders a track.
\param r Render object
\param key CString with track description.
\param value Pointer to Track object.
\param data Render object.
\param env Pointer to Environment object.
*/
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
	                           r->margins,
														 r->y,
														 config_get_color(r->cfg, "track_title", env),
														 (const char*) key);
  r->y += graphics_get_text_height(r->g) + 10;

  /* render each line */
  for (i=0; i<array_size(lines); i++)
  {
    render_line(r, *(Line**) array_get(lines, i), env);
  }

  /* put track spacer after track */
  r->y += config_get_num(r->cfg, "format", "track_vspace", 10, env);
  return 0;
}

/*!
Formats a given position number for short display in the ruler.
\param txt String to write result into.
\param pos Position.
*/
void format_ruler_label(char* txt, long pos)
{
/*  if (pos >= 1000)
	  sprintf(txt, "%.1lfk", pos/1000.0);
  else */
    sprintf(txt, "%li", pos);
}

void render_ruler(Render *r, Env* env)
{
  double step, minorstep, vmajor, vminor;
  long base_length, tick;
	Color rulercol, gridcol;
	char str[BUFSIZ];
	
	rulercol.red = rulercol.green = rulercol.blue = .2;
	gridcol.red = gridcol.green = gridcol.blue = .9;
	
	/* determine range and step of the scale */
  base_length = range_length(r->range);

  /* determine tick steps */
  step = pow(10,ceil(log10(base_length))-1);
  minorstep = step/4.0;

  /* calculate starting positions */
  vminor = (double) (floor(r->range.start / minorstep))*minorstep;
  vmajor = (double) (floor(r->range.start / step))*step;

  /* draw major ticks */
  for(tick = vmajor; tick <= r->range.end; tick += step)
  {
	  if (tick < r->range.start) continue;
    graphics_draw_vertical_line(r->g, render_convert_point(r, tick), 30, rulercol, 10);
    format_ruler_label(str, tick);
    graphics_draw_text_centered(r->g, render_convert_point(r, tick), 20, str);
  }
	
	/* draw minor ticks */
  if (minorstep >= 1)
  {
    for(tick = vminor; tick <= r->range.end; tick += minorstep)
    {
		  if (tick < r->range.start) continue;
			if(strcmp(config_get_cstr(r->cfg, "format","show_grid", "no", env), "yes") == 0)
			  graphics_draw_vertical_line(r->g, render_convert_point(r, tick), 40, gridcol, r->height);
			graphics_draw_vertical_line(r->g, render_convert_point(r, tick), 35, rulercol, 5);
    }
  }
}

void render_to_png(Render *r, char *fn, unsigned int width, Env *env)
{
  assert(r && fn && env && width > 0);
  unsigned int height;  

  /* set initial margins, header, target width */
  r->y = 70;
  r->width = width;
  height = render_calculate_height(r, env);
  r->height = height;
	
  /* calculate scaling factor */
  r->factor = ((double) r->width
                 -(2*r->margins)) 
               / (r->range.end 
                 - r->range.start);

  /* create new Graphics backend */
  r->g = graphics_new_png(fn, width, height, env);
  graphics_set_margins(r->g, r->margins, 0, width, height);
  
  render_ruler(r, env);

  /* process (render) each track */
  Hashtable *tracks = diagram_get_tracks(r->dia);
  hashtable_foreach(tracks, render_track, r, env);

  /* write out result file */
  graphics_save(r->g);
  graphics_delete(r->g, env);
}
