print ([[

Two feature node graphs are considered to be repeated if they have the same
depth-first traversal and each feature node pair is ``similar''.

Two feature nodes are ``similar'', if they have the same sequence ID, feature
type, range, strand, and phase.

For such a repeated feature node graph the one with the higher score (of the
top-level feature) is kept. If only one of the feature node graphs has a defined
score, this one is kept.]])
