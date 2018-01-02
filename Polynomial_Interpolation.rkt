#lang racket
(module+ test (require rackunit))
(require test-engine/racket-tests)

(require "Matrix_List_Helpers.rkt"
         "Linear_Systems.rkt"
         "Polynomial_Ops.rkt")


;; A Polynomial is represented as a list of coefficients, in ascending degree
;; Eg. (list 1 2 3) corresponds to 1+2x+3x^2

;; A Matrix is:
;; (list (list x0 x1 x2 ...xn)
;;       (list x0 x1 x2 ...xn)
;;       (list x0 x1 x2 ...xn) ...)
;; Interpretation: A Matrix represents a matrix where each sub-list is a row



;; [Listof (list Number Number)] Number -> Number
;; Uses Polynomial Interpolation to evaluate p(x) given a number of points of p
(module+ test
  (check-equal? (poly-interp (list '(0 0) '(2 2)) 1) 1) ;y=x
  (check-equal? (poly-interp (list '(-3 9) '(0 0) '(3 9)) 2) 4) ;y=x^2
  (check-equal? (poly-interp (list '(-2 4) '(1 1) '(8 64)) 3) 9);y=x^2
  (check-equal? (poly-interp (list '(-1 0) '(1 2) '(2 9) '(3 28)) 4) 65);y=1+x^3
  (check-equal? (poly-interp (list '(-0.5 10.75) '(-2 31) '(-2.5 58.25) '(2.5 112) '(3 226)) 5.5) 2292.25)) ;y=7-8x+3x^3+2x^4

(define (poly-interp points x)
  (local ((define f (calc-poly points)))
    (eval-poly f x)))


;; [Listof (list Number Number)] -> [Listof Number]
;; Fits a polynomial to the list of points, degree = length-1
(module+ test
  (check-equal? (calc-poly (list '(0 0) '(2 2))) '(0 1)) ;y=x
  (check-equal? (calc-poly (list '(-3 9) '(0 0) '(3 9))) '(0 0 1)) ;y=x^2
  (check-equal? (calc-poly (list '(-2 4) '(1 1) '(8 64))) '(0 0 1));y=x^2
  (check-equal? (calc-poly (list '(-1 0) '(1 2) '(2 9) '(3 28))) '(1 0 0 1));y=1+x^3
  (check-equal? (calc-poly (list '(-0.5 10.75) '(-2 31) '(-2.5 58.25) '(2.5 112) '(3 226)))
                '(7.0 -8.0 0.0 3.0 2.0))) ;y=7-8x+3x^3+2x^4
(define (calc-poly lop)
  (local ((define syst (create-system lop))
          (define A (first syst))
          (define rhs (second syst)))
    (linsys-solver A rhs)))


;; [Listof (list Number Number)] -> (list Matrix [Listof Number])
;; Turns a list of points into a system of linear equations and rhs pair
(module+ test
  (check-equal? (create-system (list '(0 0) '(2 2)))
                (list '((1 0) (1 2)) '(0 2)))
  (check-equal? (create-system (list '(-3 9) '(0 0) '(3 9)))
                (list '((1 -3 9) (1 0 0) (1 3 9)) '(9 0 9)))
  (check-equal? (create-system (list '(-1 0) '(1 2) '(2 9) '(3 28)))
                (list '((1 -1 1 -1) (1 1 1 1) (1 2 4 8) (1 3 9 27)) '(0 2 9 28))))
(define (create-system lop)
  (local ((define template (make-list (length lop) 1))
          (define list-of-degs (build-list (length lop) values))
          (define make-row (λ (p) (map (λ (d) (expt (first p) d)) list-of-degs)))
          (define A (map make-row lop))
          (define rhs (map (λ (p) (second p)) lop)))
    (list A rhs)))