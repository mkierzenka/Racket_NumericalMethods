#lang racket

(require "Continued_Fractions.rkt")

;;-- Using the function real->decimal-string we can get fractions
;;--   accurate to a number of decimal places, but not expressions


;;Exploring real->decimal-string function, checked with wolfram alpha

;(real->decimal-string 2/3 30)
;"0.666666666666666666666666666667" -- correct (rounded last)
;(real->decimal-string 7/17 30)
;"0.411764705882352941176470588235" -- correct
;(real->decimal-string 13/17 30)
;"0.764705882352941176470588235294" -- correct 
;(real->decimal-string (sqrt 2) 30)
;"1.414213562373095145474621858739" -- wrong, only 16 decimals accurate


;;-- But we can use Newton's Method, if it is a nice function
;;-- Newton's Method approx doubles precision each iteration
;;-- For example, to get sqrt(2) to 30 decimal places, we only need a small
;;--   number of iterations (one if you start from 16 decimal accuracy)

(define sqrt-2 (get-rational (get-frac (sqrt 2) (expt 10 -16))))
(define two-30
  (real->decimal-string (- sqrt-2 (/ (- (sqr sqrt-2) 2) (* 2 sqrt-2))) 30))
;;two-30 is (sqrt 2) accurate to 30 decimals (checked with wolfram alpha)

;;-- So we have some decimal arbitrary precision

;;-- But what if we were to look at something like e^1, how could we get
;;--   that to 30 decimal accuracy.

(define newton (λ (x) (- x (/ (- (sqr x) 2) (* 2 x)))))

(define newton1 (λ (x) (- x (/ (exp x) (exp x)))))


;; Number Number -> [Listof Number]
;; Returns a list of 15 successive newton approximations
(define (my-rec a n)
  (cond
    [(> n 15) '()]
    [else (cons (real->decimal-string a 30) (my-rec (newton a) (add1 n)))]))
(my-rec 1 0)