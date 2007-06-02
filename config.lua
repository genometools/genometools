-- This is the GenomeTools annotation viewer config file.
-- All options must be set inside the 'config' table.
-- All items in this table must be tables, called 'sections'.
-- The items in the 'colors' section must be given in
-- the following form:
--      <name> = {red=<val>,green=<val>,blue=<val>}
-- where <val> is a decimal value between 0 and 1.

config = 
{
  -- Defines a color for a certain feature type.
  colors = 
  {
    stroke  ={red=0.0,green=0.0,blue=0.0},
    exon    ={red=0.1,green=0.2,blue=0.3},
    cds     ={red=0.4,green=0.5,blue=0.6},
  },
  -- Defines how a feature is drawn. 
  -- Possible choices: "line", "box", "caret", "dashes"
  feature_styles = 
  {
    exon    = "box",
    CDS     = "box",
    tf_binding_factor = "box",
    mRNA    = "box",
    gene    = "line",
  },
  -- Defines how the free space in a collapsing feature track is drawn. 
  -- Possible choices: "line", "box", "caret", "dashes"
  space_styles = 
  {
    exon    = "caret",
    CDS     = "dashes",
  },
  -- Defines which feature types are displayed in another 
  -- feature's track.
  collapse = 
  {
    to_parent = {"exon", "CDS"},
  },
  -- Defines various format options for drawing.
  format =
  {
    margins = 10,
    bar_height = 8,
    bar_vspace = 5,
    arrow_width = 6,
    stroke_width = 1,
  },
} 
