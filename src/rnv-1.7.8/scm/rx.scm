; $Id: rx.scm,v 1.15 2004/01/29 18:14:36 dvd Exp $

; XML Schema Datatypes Regular Expressions
;   This is a full conformant implementation of the regular expressions
;   syntax in Scheme;   the only implementation-dependant feature I am aware of
;   is integer<->char conversions. Normal people (including SCM) use byte value;
;   and everything works. scheme48 adds 1000 to the byte value to make the code
;   general. To make it work with scheme48, redefine char->integer,integer->char
;   to char->ascii,ascii->char.

; Synopsis:
;   (rx-compile "string")
;     if the computed value is true, it is a pattern suitable for passing to
;   (rx-match pattern "string")
;     the returned value is true if the string matches the pattern

; Example
;   (rx-match (rx-compile "\\c\\i+") "fo:block")

; Performance:
;   it is approximately three times slower than 'regex. It memoizes patterns
;   to speed up processing, however, all non-core algorithms are as simple as
;   possible.

(load  (in-vicinity (program-vicinity) "u.scm"))
(load  (in-vicinity (program-vicinity) "xml-ranges.scm"))
(load  (in-vicinity (program-vicinity) "rx-ranges.scm"))

(define rx-patterns '(none empty any choice group more except range class char))

(define rx-regex-cache-size 16)
(define rx-pattern-cache-size 256)
(define rx-memo-cache-size 1024)

(define (rx-null? p)
  (case (car p)
    ((empty) #t)
    ((choice) (or (rx-null? (cadr p)) (rx-null? (cddr p))))
    ((group) (and (rx-null? (cadr p)) (rx-null? (cddr p))))
    ((more) (rx-null? (cdr p)))
    (else #f)))

(define rx-none #f)
(define rx-empty #f)
(define rx-any #f)

(define rx-regex-cache (make-vector rx-regex-cache-size #f))
(define (rx-regex-hash regex)
  (let digit ((v 0) (l (string->list regex)))
    (if (null? l) v
      (digit (+ (* v 31) (char->integer (car l))) (cdr l)))))

(define rx-memo-cache (make-vector rx-memo-cache-size #f))

(define rx-newpat ; #f flushes
  (letrec (
      (cache '())
      (pat=?
	(lambda (p1 p2)
	  (and (equal? (car p1) (car p2))
	    (case (car p1)
	      ((none empty any) #t)
	      ((choice group except range)
		(and (eqv? (cadr p1) (cadr p2)) (eqv? (cddr p1) (cddr p2))))
	      ((more class char)
		(eqv? (cdr p1) (cdr p2)))))))
      (old
	(lambda (p que)
	  (and (pair? que)
	    (or (and (pat=? p (car que)) (car que))
	      (old p (cdr que)))))))
    (lambda (p)
      (if p
        (let ((que (assv (car p) cache)))
	  (if que
	    (or (old p (cdr que))
	      (begin (set-cdr! que (cons p (cdr que))) p))
	    (begin (set! cache (cons (list (car p) p) cache)) p)))
	(begin
	  (if (> (apply + (map length cache)) rx-pattern-cache-size) (set! cache '()))
	  (set! rx-none (rx-newpat '(none)))
	  (set! rx-empty (rx-newpat '(empty)))
	  (set! rx-any (rx-newpat '(any)))
	  #f)))))	
	
(define (rx-more p)
  (case (car p)
    ((none empty more) p)
    (else (rx-newpat `(more . ,p)))))

(define rx-choice
  (letrec (
      (same
	(lambda (p1 p2)
	  (or
	    (and (eqv? 'choice (car p1))
	      (or (eqv? (cddr p1) p2) (same (cadr p1) p2)))
	    (eqv? p1 p2)))))
    (lambda (p1 p2)
      (or
	(and (eq? rx-none p1) p2)
	(and (eqv? rx-none p2) p1)
	(and (eqv? 'choice (car p2))
	  (rx-choice (rx-choice p1 (cadr p2)) (cddr p2)))
	(and (same p1 p2) p1)
	(and (rx-null? p1) (eq? rx-empty p2) p1)
	(and (rx-null? p2) (eq? rx-empty p1) p2)
	(rx-newpat `(choice ,p1 . ,p2))))))

(define (rx-group p1 p2)
  (cond
    ((or (eq? rx-none p1) (eq? rx-none p2)) rx-none)
    ((eq? rx-empty p1) p2)
    ((eqv? rx-empty p2) p1)
    (else (rx-newpat `(group ,p1 . ,p2)))))

; parser symbols: chr esc cls ncl end

(define (rx-compile regex)
  (letrec (
      (rxll (utf8->lazy-list regex))
      (nextc
	(lambda ()
	  (and (pair? rxll)
	    (let ((c (car rxll))) (set! rxll (force (cdr rxll))) c))))
      (newpat
	(let ((cache '()))
	  (lambda (pattern)
	    (car
	      (or (member pattern cache)
		(begin (set! cache (cons pattern cache)) cache))))))
      (errors #f)
      (error!
	(lambda msg
	  (or errors
	    (begin
	      (display (string-append "bad regex '" regex "'"))
	      (if (pair? msg)
		(begin (display ":")
		  (for-each (lambda (x) (display " ") (display x)) msg)))
	      (newline)
	      (set! errors #t)))))
      (chclass
	(lambda ()
	  (if (not (= (nextc) 123)) (begin (error! "{ expected") #f)
	    (let i ((c (nextc)) (l '()))
	      (cond
		((not c) (error! "class name expected") #f)
		((>= c 128) (error! "illegal class name") #f)
		((= c 125)
		  (let (
		     (clsym
		       (string->symbol
			 (list->string
			  `(#\u #\-
			 ,@(map (lambda (i)
			     (char-downcase (integer->char i)))
			     (reverse l)))))))
		   (or (assv clsym rx-ranges) (error! "unknown class" clsym))
		   clsym))
	       (else (i (nextc) (cons c l))))))))
      (esc
	(lambda ()
	  (let ((c (nextc)))
	    (if (not c) (begin (error! "escape expected") '(end))
	      (if (>= c 128) (begin (error! "illegal escape") `(esc ,c))
		(case (integer->char c)
		  ((#\p) `(cls . ,(chclass)))
		   ((#\P) `(ncl . ,(chclass)))
		  ((#\s) '(cls . S))
		   ((#\S) '(ncl . S))
		  ((#\i) '(cls . I))
		   ((#\I) '(ncl . I))
		  ((#\c) '(cls . C))
		   ((#\C) '(ncl . C))
		  ((#\d) '(cls . U-Nd))
		   ((#\D) '(ncl . U-Nd))
		  ((#\w) '(cls . W))
		   ((#\W) '(ncl . W))
		  ((#\n) '(esc . 10))
		  ((#\r) '(esc . 13))
		  ((#\t) '(esc . 9))
		  ((#\\ #\| #\. #\- #\^ #\? #\* #\+
		      #\{ #\} #\[ #\] #\( #\))
		   `(esc . ,c))
		  (else (error! "unknown escape") `(esc . ,c))))))))
      (sym #f)
      (getsym
	(lambda ()
	  (set! sym
	    (let ((c (nextc)))
	      (or (and (not c) '(end))
		(or (and (= c 92) (esc))
		  (and (= c 46) '(ncl . NL)))
	       `(chr . ,c))))))
      (chgroup
	(let (
	    (check-range
	      (lambda ()
		(and (eqv? (cdr sym) 'chr) (memv (cdr sym) '(45 91 93))
		  (error! "illegal range") #t))))
	  (lambda()
	    (let range ((p rx-none))
	      (if
		 (and (not (eqv? (car p) 'none))
		   (and (or (equal? sym '(chr . 45)) (equal? sym '(chr . 93)))))
		 p
	      (case (car sym)
		((chr esc)
		  (check-range)
		  (let ((c (cdr sym))) (getsym)
		    (if (equal? sym '(chr . 45))
		      (if (and (pair? rxll) (eqv? (car rxll) 91))
			(rx-choice p (rx-newpat `(char . ,c)))
			(begin (getsym)
			  (case (car sym)
			    ((chr esc) (check-range)
			      (let ((p (rx-choice p (rx-newpat `(range ,c . ,(cdr sym))))))
				(getsym) (range p)))
			    (else (error! "illegal range") (getsym) (range p)))))
		      (range (rx-choice p (rx-newpat `(char . ,c)))))))
	        ((cls)
	          (let ((c (cdr sym)))
		    (getsym) (range (rx-choice p (rx-newpat `(class . ,c))))))
	        ((ncl)
	          (let ((c (cdr sym)))
		    (getsym)
		    (range (rx-choice p
		      (rx-newpat `(except ,rx-any . ,(rx-newpat `(class . ,c))))))))
	        (else (error! "missing ]") (getsym) p)))))))
      (chexpr
	(lambda()
	  (let (
	      (p
		(if (equal? sym '(chr . 94))
		  (begin (getsym) (rx-newpat `(except ,rx-any . ,(chgroup))))
		  (chgroup))))
	    (if (equal? sym '(chr . 45))
	      (begin (getsym)
		(or (equal? sym '(chr . 91)) (error! "[ expected"))
		(getsym)
		(let ((p (rx-newpat `(except ,p . ,(chgroup)))))
		  (getsym)
		  (or (equal? sym '(chr . 93)) (error! "] expected"))
		  p))
	      p))))
      (atom
	(lambda ()
	  (case (car sym)
	    ((chr)
	      (case (cdr sym)
		((91) ; [
		  (getsym)
		  (let ((p (chexpr)))
		    (or (equal? sym '(chr . 93)) (error! "missing ]"))
		    (getsym)
		    p))
		((40) ; (
		  (getsym)
		  (let ((p (expr)))
		    (or (equal? sym '(chr . 41)) (error! "missing )"))
		    (getsym)
		    p))
		((41 42 43 63 93 123 124 125) ; must be escaped: () * + ? ] { | }
		  (error! "unescaped " (integer->char (cdr sym)))
		  (getsym) rx-none)
		(else
		  (let ((p (rx-newpat `(char . ,(cdr sym))))) (getsym) p))))
	    ((esc)
	      (let ((p (rx-newpat `(char . ,(cdr sym))))) (getsym) p))
	    ((cls)
	      (let ((p (rx-newpat `(class . ,(cdr sym))))) (getsym) p))
	    ((ncl)
	      (let ((p (rx-newpat `(except ,rx-any .
		      ,(rx-newpat `(class . ,(cdr sym))))))) (getsym) p))
	    (else (error! sym) (getsym) rx-none))))
      (number
	(lambda ()
	  (let digit ((n 0))
	    (or
	      (and (eqv? (car sym) 'chr)
		(let ((d (cdr sym)))
		  (and (<= d 57) (>= d 48)
		    (begin (getsym) (digit (+ (* n 10) (- d 48)))))))
	       n))))
      (quantifier
	(lambda (p0)
	  (let ((n0 (number)))
	    (let from ((p rx-empty) (n n0))
	      (or (and (> n 0) (from (rx-group p p0) (- n 1)))
		(and (eqv? (car sym) 'chr)
		  (or (and (eqv? (cdr sym) 125) p)
		    (and (eqv? (cdr sym) 44)
		      (begin
			(getsym)
			(if (and (eqv? (car sym) 'chr) (eqv? (cdr sym) 125))
			  (rx-group p
			    (rx-choice rx-empty (rx-more p0)))
			  (let till ((p p) (n (- (number) n0)))
			    (if (> n 0)
			      (till
				(rx-group p
				  (rx-choice rx-empty p0))
				(- n 1))
				p)))))))
		 (begin (error! "bad quantifier") p))))))
      (piece
	(lambda ()
	  (let ((p (atom)))
	    (if (eqv? (car sym) 'chr)
	      (case (cdr sym)
		((63) ; ?
		  (getsym) (rx-choice p rx-empty))
		((42) ; *
		  (getsym) (rx-choice (rx-more p) rx-empty))
		((43) ; +
		  (getsym) (rx-more p))
		((123) ; {
		  (getsym)
		    (let ((p (quantifier p)))
		      (or (equal? sym '(chr . 125)) (error! "missing }"))
		      (getsym)
		      p))
		(else p))
	      p))))
      (branch
	(lambda ()
	  (let loop ((p rx-empty))
	    (if(or (eqv? (car sym) 'end)
	      (and (eqv? (car sym) 'chr)
		(or (eqv? (cdr sym) 124) (eqv? (cdr sym) 41))))
	      p
	      (loop (rx-group p (piece)))))))
      (expr
	(lambda ()
	  (let loop ((p (branch)))
	    (if (equal? sym '(chr . 124))
	      (begin (getsym)
		(loop (rx-choice p (branch))))
	      p)))))
    (let* (
        (n (remainder (rx-regex-hash regex) (vector-length rx-regex-cache)))
	(cell (vector-ref rx-regex-cache n)))	
      (if (and cell (equal? (car cell) regex)) (cdr cell)
        (begin
          (rx-newpat #f) (getsym) 	
          (let ((p (expr)))
            (or (equal? sym '(end)) (error! "junk after end"))
	    (vector-set! rx-regex-cache n (cons regex p))
            (and (not errors) p)))))))

(define (rx-deriv p c)
  (letrec (
      (in-class?
	(letrec (
	    (or*
	      (lambda list
		(and (pair? list)
		  (or (car list) (apply or* (cdr list))))))
	    (in-xml-ranges (lambda (i) (u-in-ranges c (cdr (assv i xml-ranges)))))
	    (in-rx-ranges (lambda (i) (u-in-ranges c (cdr (assv i rx-ranges)))))
	    (join-rx-ranges (lambda (ranges) (apply or* (map in-rx-ranges ranges)))))
	  (lambda (c id)
	    (case id
	      ((NL) (or (= c 13) (= c 10)))
	      ((S) (or (= c 13) (= c 10) (= c 9) (= c 32)))
	      ((I) (or (= c 58) (= c 95)
		     (apply or* (map in-xml-ranges '(base-char ideographic)))))
	      ((C) (or (in-class? c 'I) (= c 45) (= c 46)
		     (apply or* (map in-xml-ranges '(digit combining-char extender)))))
	      ((W) (not
		     (or (in-class? c 'U-P) (in-class? c 'U-Z) (in-class? c 'U-C))))
	      ((U-C) (join-rx-ranges '(U-Cc U-Cf U-Co)))
	      ((U-L) (join-rx-ranges '(U-Lu U-Ll U-Lt U-Lm U-Lo)))
	      ((U-M) (join-rx-ranges '(U-Mn U-Mc U-Me)))
	      ((U-N) (join-rx-ranges '(U-Nd U-Nl U-No)))
	      ((U-P) (join-rx-ranges '(U-Pc U-Pd U-Ps U-Pe)))
	      ((U-S) (join-rx-ranges '(U-Sm U-Sc U-Sk U-So)))
	      ((U-Z) (join-rx-ranges '(U-Zl U-Zs U-Zp)))
	      (else (in-rx-ranges id)))))))
    (let* (
        (n
	  (remainder (+ (* c (length rx-patterns)) (length (memv (car p) rx-patterns)))
	    (vector-length rx-memo-cache)))
	(cell (vector-ref rx-memo-cache n)))
      (if (and cell (eqv? (car cell) p) (eqv? (cadr cell) c)) (cddr cell)
        (let ((q
            (case (car p)
              ((none empty) rx-none)
              ((any) rx-empty)
              ((choice) (rx-choice (rx-deriv (cadr p) c) (rx-deriv (cddr p) c)))
              ((group)
                (let ((q (rx-group (rx-deriv (cadr p) c) (cddr p))))
                  (if (rx-null? (cadr p)) (rx-choice q (rx-deriv (cddr p) c)) q)))
              ((more) (rx-group (rx-deriv (cdr p) c) (rx-choice rx-empty p)))
              ((except)
                (if (and (rx-null? (rx-deriv (cadr p) c))
                    (not (rx-null? (rx-deriv (cddr p) c))))
                  rx-empty
                  rx-none))
              ((range) (if (and (<= (cadr p) c) (<= c (cddr p))) rx-empty rx-none))
              ((class) (if (in-class? c (cdr p)) rx-empty rx-none))
              ((char) (if (= (cdr p) c) rx-empty rx-none)))))
	  (vector-set! rx-memo-cache n (cons p (cons c q)))
	  q)))))

(define (rx-match r s)
  (and r
    (let ((sll (utf8->lazy-list s)))
      (let delta ((p r) (sll sll))
	(if (null? sll) (rx-null? p)
	  (let* (
	      (c (car sll))
	      (p (rx-deriv p c)))
	    (and (not (eq? p rx-none))
	      (delta p (force (cdr sll))))))))))
