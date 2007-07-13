-- This is the GenomeTools annotation viewer config file.
-- All options must be set inside the 'config' table.
-- All items in this table must be tables, called 'sections'.

config =
{
  -- Defines a color for a certain feature type.
  -- The items in the 'colors' section must be given in
  -- the following form:
  --      <name> = {red=<val>,green=<val>,blue=<val>}
  -- where <val> is a decimal value between 0 and 1.
  colors =
  {
    stroke          ={red=0.0,green=0.0,blue=0.0},
    track_title     ={red=0.6,green=0.6,blue=0.7},
    exon            ={red=0.6,green=0.6,blue=0.9},
    CDS             ={red=0.9,green=0.9,blue=0.2},
    mRNA            ={red=0.4,green=0.5,blue=0.6},
    TF_binding_site ={red=0.8,green=0.6,blue=0.6},
    gene            ={red=0.9,green=0.9,blue=1.0},
    intron          ={red=0.2,green=0.2,blue=0.6},
    repeat_region   ={red=0.8,green=0.3,blue=0.3},
    long_terminal_repeat ={red=0.9,green=0.9,blue=0.4},
    LTR_retrotransposon  ={red=0.8,green=0.5,blue=0.5},
  },
  -- Defines how a feature is drawn. 
  -- Possible choices: "line", "box", "caret", "dashes"
  feature_styles =
  {
    exon            = "box",
    CDS             = "box",
    TF_binding_site = "box",
    mRNA            = "box",
    gene            = "box",
    intron          = "caret",
  },
  -- Defines which feature types are displayed in another 
  -- feature's track.
  collapse =
  {
    to_parent = {"exon","intron","CDS"},
  },
  -- Defines precedence of feature types when overlapping
  -- in a collapsed parent track.
  -- read "=" as ">" or "dominates"
  dominate =
  {
    CDS = {"exon","intron","mRNA","gene"},
    exon = {"mRNA","gene"},
    intron = {"mRNA", "gene"},
  },
  -- Defines various format options for drawing.
  format =
  {
    margins = 30,      -- space left and right of diagram, in pixels
    bar_height = 15,   -- height of a feature bar, in pixels
    bar_vspace = 10,   -- space between feature bars, in pixels
    track_vspace = 10, -- space between tracks, in pixels
    scale_arrow_width = 6,     -- width of scale arrowheads, in pixels
    scale_arrow_height = 10,   -- height of scale arrowheads, in pixels
    arrow_width = 6,   -- width of feature arrowheads, in pixels
    stroke_width = .5,  -- width of outlines, in pixels
    show_grid = "yes", -- shows light vertical lines for orientation
  },
}
