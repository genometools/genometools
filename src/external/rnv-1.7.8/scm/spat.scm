; $Id: spat.scm,v 1.5 2004/01/31 15:45:40 dvd Exp $

(define (spat->regex spat)
  (letrec (
      (l (string->list spat))
      (errors #f)
      (error!
	(lambda msg
	  (set! errors #t)
	  (for-each display `("bad s-pattern '" ,spat "': " ,@msg))
	  (newline)))
      (nextc
	(lambda ()
	  (and (pair? l) (let ((c (car l))) (set! l (cdr l)) c))))
      (nextsym
	(let (
	    (literal
	      (lambda ()
		(let ch ((ll '()) (c (nextc)))
		  (case c ((#f) (error!) "")
		    ((#\") (list->string (reverse ll)))
		    ((#\\) (let ((d (nextc)))
			(ch (if (eqv? d #\")
			    (cons d ll)
			    (cons d (cons c ll)))
			  (nextc))))
		    (else (ch (cons c ll) (nextc)))))))
	    (ident
	      (lambda (c)
		(string->symbol (list->string (reverse
		  (let ch ((il (list c)) (c (nextc)))
		    (cond
		      ((or (not c) (char-whitespace? c)) il)
		      ((char=? #\= c) (set! l (cons c l)) il)
		      (else (ch (cons c il) (nextc)))))))))))
	  (lambda ()
	    (let ((c (nextc)))
	      (and c
		(cond
		  ((char=? #\" c) `(lit . ,(literal)))
		  ((char=? #\= c) '(=))
		  ((char-whitespace? c) (nextsym))
		  (else `(id . ,(ident c)))))))))
      (code
	(lambda ()
	  (let tok ((prog '()) (ding '(=)) (sym (nextsym)))
	    (if sym
	      (case (car sym)
		((lit) (tok prog (cons (cdr sym) ding) (nextsym)))
		((id)
		  (let ((sym2 (nextsym)))
		    (if (and sym2 (eqv? (car sym2) '=))
		      (tok (cons (reverse ding) prog) (list (cdr sym)) (nextsym))
		      (tok prog (cons (cdr sym) ding) sym2))))
		(else (error! sym " unexpected") #f))
	      (cons (reverse ding) prog)))))
       (splice
	 (lambda (code)
	   (letrec (
	       (resolve
		 (lambda (piece parents)
		   (cond
		     ((string? piece) piece)
		     ((symbol? piece)
		       (if (memv piece parents)
			 (begin (error! "recursion through " piece) "")
			 (let ((entry (assv piece code)))
			   (if entry
			     (apply string-append
			       (map (lambda (p)
				   (let ((parents (cons piece parents)))
				     (resolve p parents)))
				 (cdr entry)))
			     (begin (error! piece " unresolved") "")))))))))
	     (resolve 'start '())))))
    (let ((regex (splice (code))))
      (and (not errors) regex))))
