derives_from = {}

describe.feature("gene", function(gene)
  it("contains a transcript", function()
    expect(gene:has_child_of_type("mRNA")
             or gene:has_child_of_type("tRNA")
             or gene:has_child_of_type("snRNA")
             or gene:has_child_of_type("snoRNA")
             or gene:has_child_of_type("rRNA")
             or gene:has_child_of_type("SLRNA")
             or gene:has_child_of_type("ncRNA")).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in gene:get_children() do
      local gene_rng = gene:get_range()
      local child_rng = child:get_range()
      expect(child_rng:get_start() >= gene_rng:get_start()
               and child_rng:get_end() <= gene_rng:get_end()).should_be(true)
    end
  end)

  it("does not overlap an unrelated feature", function()
    fs = feature_index:get_features_for_range(gene:get_seqid(), gene:get_range())
    -- expect the gene and the polypeptide in this region
    expect(#fs).should_be_smaller_than(3)
  end)

  it("adheres to a ID naming scheme, if coding", function()
    if gene:has_child_of_type("mRNA") then
      expect(gene:get_attribute("ID")).should_match("^"..gene:get_seqid()
                                                       ..".%d%d%d%d")
    else
      expect(gene:get_attribute("ID")).should_match("RNA")
    end
  end)

  it("is not suspiciously short (>30nt)", function()
    local rng = gene:get_range()
    expect(rng:get_end() - rng:get_start() + 1).should_be_larger_than(30)
  end)

  it("does not span a contig separator sequence (100 Ns)", function()
    local seq = string.lower(gene:extract_sequence("gene", false, region_mapping))
    expect(seq).should_not_match(string.rep("n", 100))
  end)
end)

describe.feature("mRNA", function(mrna)
  it("has a coding sequence", function()
    expect(mrna:has_child_of_type("CDS")).should_be(true)
  end)

  it("has CDS with no internal stop codons", function()
    local seq = mrna:extract_sequence("CDS", true, region_mapping)
    expect(gt.translate_dna(string.sub(seq, 1, -3))).should_not_match("\*")
  end)

  it("has CDS ending on a stop codon", function()
    local seq = mrna:extract_sequence("CDS", true, region_mapping)
    expect(gt.translate_dna(string.sub(seq, -3, -1))).should_match("\*")
  end)
end)

describe.feature("polypeptide", function(pp)
  it("is unique in its Derives_from", function()
    local dfrom = pp:get_attribute("Derives_from")
    expect(dfrom).should_not_be(nil)
    expect(derives_from).should_not_have_key(dfrom)
    derives_from[dfrom] = true
  end)

  it("has a product name associated with it", function()
    expect(pp:get_attribute("product")).should_not_be(nil)
    if (pp:get_attribute("product") ~= nil) then
      expect(string.len(pp:get_attribute("product"))).should_be_larger_than(0)
    end
  end)
end)

describe.meta(function(meta)
  it("only uses valid URLs in the feature-ontology directive", function()
    if meta:get_directive() == "feature-ontology" then
      url = meta:get_data()
      expect(string.find(url, "https?://[%w-_%.%?%.:/%+=&]+") or
             string.find(url, "ftp://[%w-_%.%?%.:/%+=&]+")).should_be_truthy()
    end
  end)
end)

describe.region(function(region)
  it("has sequence IDs starting with organism shorthand", function()
    expect(region:get_seqid()).should_match("^LmjF.%d%d$")
  end)

  it("starts at coordinate 1", function()
    expect(region:get_range():get_start()).should_be(1)
  end)
end)

describe.comment(function(comment)
  it("is not longer than 80 characters", function()
    expect(string.len(comment:get_comment())).should_be_smaller_than(80)
  end)
end)

describe.sequence(function(seq)
  it("has descriptions starting with organism shorthand", function()
    expect(seq:get_description()).should_match("^LmjF.")
  end)

  it("is not empty", function()
    expect(seq:get_sequence_length()).should_be_larger_than(0)
  end)
end)
