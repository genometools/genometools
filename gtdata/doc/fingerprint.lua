print ([[

If neither option -check nor option -duplicates is used, the fingerprints for
all sequences are shown on stdout.

Examples:
---------

Make sure a sequence file contains no duplicates (not the case here):

$ gt fingerprint -duplicates U89959_ests.fas
97350772e29b7881b149afbedf09290d 2
gt fingerprint: error: duplicates found: 1 out of 199 (0.503%)

Extract sequence with given fingerprint:

$ gt fingerprint -extract 6d3b4b9db4531cda588528f2c69c0a57 U89959_ests.fas
>SQ;8720010
TTTTTTTTTTTTTTTTTCCTGACAAAACCCCAAGACTCAATTTAATCAATCCTCAAATTTACATGATACCAACGTAATGGGAGCTTAAAAATA

Return values:
--------------

0  everything went fine (-check: the comparison was successful;
                         -duplicates: no duplicates found)
1  an error occured     (-check: the comparison was not successful;
                         -duplicates: duplicates found)]])
