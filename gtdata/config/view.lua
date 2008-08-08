-- This is the GenomeTools annotation viewer config file.
-- All options must be set inside the 'config' table.
-- All items in this table must be tables, called 'sections'.

config =
{
  gene = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.9, green=0.9, blue=1.0},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = false,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  mRNA = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.4, green=0.5, blue=0.6},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = false,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  exon = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.6, green=0.6, blue=0.9},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  CDS = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.9, green=0.9, blue=0.2},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  TF_binding_site = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.8, green=0.6, blue=0.6},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = false,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  intron = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=1.0, green=1.0, blue=1.0},
    style              = "caret",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  repeat_region = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.8, green=0.3, blue=0.3},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = false,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  LTR_retrotransposon = {
    -- Color definitions
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.8, green=0.5, blue=0.5},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  long_terminal_repeat = {
    -- Color definitions
    -- RGB triplets {red=<val>,green=<val>,blue=<val>}
    -- where <val> is a decimal value between 0 and 1.
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.9, green=0.9, blue=0.4},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
    protein_match = {
    -- Color definitions
    -- RGB triplets {red=<val>,green=<val>,blue=<val>}
    -- where <val> is a decimal value between 0 and 1.
    stroke             = {red=0.0, green=0.0, blue=0.0},
    stroke_marked      = {red=1.0, green=0.0, blue=0.0},
    fill               = {red=0.1, green=0.1, blue=0.5},
    style              = "box",
    -- Collapsing options
    collapse_to_parent = true,
    split_lines        = true,
    -- Caption options
    max_capt_show_width= nil,
    -- Display this track only if the viewport is not wider than this
    -- number of nucleotides. Set to 0 to disable type track.
    max_show_width     = nil,
    -- Limit the number of tracks
    max_num_lines      = 10,
  },
--------------------------------------
  five_prime_splice_site = {
    collapse_to_parent = true,
  },
--------------------------------------
  three_prime_splice_site = {
    collapse_to_parent = true,
  },
--------------------------------------
  expressed_sequence_match = {
    fill               = {red=0.2, green=0.2, blue=0.8},
    max_show_width     = 10000,
    max_num_lines      = 10,
    max_capt_show_width= 5000,
  },
--------------------------------------
  binding_site = {
    fill               = {red=0.7, green=0.2, blue=0.8},
    max_show_width     = 10000,
    max_num_lines      = 10,
    max_capt_show_width= 10000,
  },
--------------------------------------
  SNP = {
    fill               = {red=0.1, green=0.8, blue=0.8},
    max_show_width     = 10000,
    max_num_lines      = 10,
    max_capt_show_width= 10000,
  },
--------------------------------------
  chromosome = {
    fill               = {red=0.1, green=0.8, blue=0.8},
    max_show_width     = nil,
    max_num_lines      = 10,
  },
--------------------------------------
  substitution = {
    fill               = {red=1.0, green=0.1, blue=0.05},
    max_show_width     = nil,
    max_num_lines      = 10,
    max_capt_show_width= 1000,
  },
  -- Defines various format options for drawing.
  format =
  {
    margins = 30,      -- space left and right of diagram, in pixels
    bar_height = 15,   -- height of a feature bar, in pixels
    bar_vspace = 10,   -- space between feature bars, in pixels
    track_vspace = 20, -- space between tracks, in pixels
    scale_arrow_width = 6,     -- width of scale arrowheads, in pixels
    scale_arrow_height = 10,   -- height of scale arrowheads, in pixels
    arrow_width = 6,   -- width of feature arrowheads, in pixels
    stroke_width = .5, -- width of outlines, in pixels
    stroke_marked_width = 1.5, -- width of outlines for marked elements, in pixels
    show_grid = true, -- shows light vertical lines for orientation
    min_len_block = 5, -- minimum length of a block in which single elements are shown
    track_title_color     = {red=0.7, green=0.7, blue=0.7},
    default_stroke_color  = {red=0.1, green=0.1, blue=0.1},
  },
}
