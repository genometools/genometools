; $Id: M.scm 345 2004-01-25 09:44:04Z dvd $

(define (result v) (lambda (inp) (list (cons v inp))))
(define (zero inp) '())
(define (unit inp)
  (if (pair? inp)
    (list (cons (car inp) (force (cdr inp))))
    '()))

(define bind
  (lambda (p f)
    (lambda (inp)
      (let ((rl (p inp)))
	(apply
	  append
	  (map
	    (lambda (l) (if (pair? l) ((f (car l)) (cdr l)) '()))
	    rl))))))

(define (plus p q) (lambda (inp) (append (p inp) (q inp))))

; the rest is experiments

(define (sat p)
  (bind unit (lambda (x) (if (p x) (result x) zero))))

(define (char x)
  (sat (lambda (y) (eqv? x (integer->char y)))))

(define digit
  (sat (lambda (y) (char-numeric? (integer->char y)))))

(define upper
  (sat (lambda (y) (char-upper-case? (integer->char y)))))

(define lower
  (sat (lambda (y) (char-lower-case? (integer->char y)))))

(define letter (plus lower upper))
(define alnum (plus letter digit))

(define (many o)
  (let (
      (some
	(bind o (lambda (x)
	(bind (many o) (lambda (xs)
	  (result (cons x xs))))))))
    (plus some (result '()))))
	


