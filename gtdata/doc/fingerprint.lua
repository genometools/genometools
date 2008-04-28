print ([[

If neither option -check nor option -duplicates is used, the fingerprints for
all sequences are shown on stdout.

Examples:
---------

Make sure a sequence file contains no duplicates (not the case here):

$ gt fingerprint -duplicates testdata/U89959_ests.fas
97350772e29b7881b149afbedf09290d 2
gt fingerprint: error: duplicates found: 1 out of 199 (0.503%)

Return values:
--------------

0  everything went fine (-check: the comparison was successful;
                         -duplicates: no duplicates found)
1  an error occured     (-check: the comparison was not successful;
                         -duplicates: duplicates found)]])
