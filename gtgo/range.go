package gt

// Range represents a genomic ranges. The start must always be smaller or
// equal than the end.
type Range struct {
	Start uint
	End   uint
}
