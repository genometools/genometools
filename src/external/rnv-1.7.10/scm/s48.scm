; loader to test rx with scheme48

(define (program-vicinity) #f)
(define (in-vicinity vicinity string) string)

(define char->integer char->ascii)
(define integer->char ascii->char)

(load "u.scm")
(load "xml-ranges.scm")
(load "rx-ranges.scm")

(load "rx.scm")
