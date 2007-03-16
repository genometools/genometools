; $Id: u.scm,v 1.5 2004/01/25 09:44:04 dvd Exp $
; unicode

(define (utf8->lazy-list s) ; no checks, everything we get is unicode
  (letrec (
      (ux
	(lambda (c sl n) 
	  (if (= n 0) (cons c (delay (left sl)))
	    (ux 
	      (+ (* c 64) (remainder (char->integer (car sl)) 64)) 
	      (cdr sl) (- n 1)))))
      (left
	(lambda (sl)
	  (if (null? sl) '()
	    (let ((c (char->integer (car sl))) (sl (cdr sl)))
	      (cond 
		((< c #x80) (ux c sl 0))
		((< c #xE0) (ux (remainder c #x20) sl 1))
		((< c #xF0) (ux (remainder c #x10) sl 2))
		((< c #xF8) (ux (remainder c #x08) sl 3))
		((< c #xFC) (ux (remainder c #x04) sl 4))
		((< c #xFE) (ux (remainder c #x02) sl 5))))))))
    (left (string->list s))))      

; binary search on list of ranges
(define (u-in-ranges u ranges)
  (let between ((n 0) (m (- (vector-length ranges) 1)))
    (if (> n m) #f
      (let* ((i (quotient (+ n m) 2)) (r (vector-ref ranges i)))
	(cond
	  ((< u (car r)) (between n (- i 1)))
	  ((> u (cdr r)) (between (+ i 1) m))
	  (else r))))))
