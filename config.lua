-- This is the GenomeTools annotation viewer config file.
-- All options must be set inside the 'config' table.
-- All items in this table must be tables, called 'sections'.
-- The items in the 'colors' section must be given in
-- the following form:
--      <name> = {red=<val>,green=<val>,blue=<val>}
-- where <val> is a decimal value between 0 and 1.

config = 
{
  -- Defines which tracks are shown, independent of view range
  tracks =
  {
    hide_always = {},
  },
  threshold =
  {
    -- TODO: let here be thresholds for view ranges
    -- and the types still displayed in this view
    -- e.g. exon = 10000
    -- means exon elements are nor displayed when
    -- a range > 10000 is viewed
  },
  -- Defines a color for a certain feature type.
  colors = 
  {
    stroke          ={red=0.0,green=0.0,blue=0.0},
    track_title     ={red=0.6,green=0.6,blue=0.7},
    exon            ={red=0.7,green=0.7,blue=0.9},
    CDS             ={red=0.9,green=0.9,blue=0.2},
    mRNA            ={red=0.4,green=0.5,blue=0.6},
    TF_binding_site ={red=0.8,green=0.6,blue=0.6},
    gene            ={red=0.8,green=0.2,blue=0.2},
    intron          ={red=0.2,green=0.2,blue=0.6},
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
    --to_parent = {"exon","intron"},
  },
  dominate = 
  {
    CDS = {"exon"},
    exon = {"mRNA","gene"},
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
    stroke_width = 1,  -- width of outlines, in pixels
  },
} 
