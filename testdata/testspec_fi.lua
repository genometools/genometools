rm = gt.region_mapping_new_seqfile("LmjF_v6.1_all_20131105.fa")

derives_from = {}

describe.feature("gene", function(gene)
  it("does not overlap an unrelated feature", function()
    fs = featureidx:get_features_for_range(gene:get_seqid(), gene:get_range())
    -- expect the gene and the polypeptide in this region
    expect(#fs).should_be_smaller_than(3)
  end)
end)
