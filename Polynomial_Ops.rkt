#lang racket
(module+ test (require rackunit))

(provide eval-poly)


;;--Polynomials are represented as a list of coefficients
;;--Least significant is first
;;(list 1 2 3) is 1+2x+3x^2

;; [Listof Number] Number -> Number
;; Evaluate given polynomial at x
(module+ test
  (check-equal? (eval-poly '(1) 4) 1)
  (check-equal? (eval-poly '(1 1 1) 2) (+ 1 2 4))
  (check-equal? (eval-poly '(6 0 2) 2) (+ 6 (* 2 (sqr 2))))
  (check-equal? (eval-poly '(1 0 0) 4) 1))
(define (eval-poly coeffs0 x)
  (local (; [Listof Number] -> Number
          ; Evaluates polynomial at x, recursively
          ; Input is coeffs not yet evaluated
          (define (eval-runner coeffs)
            (cond
              [(empty? (rest coeffs)) (first coeffs)]
              [else (+ (first coeffs) (* x (eval-runner (rest coeffs))))])))
    (eval-runner coeffs0)))
          

;; [Listof Number] -> [Listof Number]
;; Return derivative of a polynomial, derivative of constant is empty list
(module+ test
  (check-equal? (deriv '(3)) '())
  (check-equal? (deriv '(1 2)) '(2))
  (check-equal? (deriv '(1 2 3)) '(2 6)))
(define (deriv coeffs0)
  (local ((define length0 (length coeffs0))
          ; [Listof Number] Number -> [Listof Number]
          ; Differentiates polynomial, recursively
          ; Input is coeffs not yet diffed, p is power of 1st coeff
          (define (deriv-runner coeffs p)
            (cond
              [(empty? coeffs) '()]
              [else (cons (deriv-term (first coeffs) p)
                          (deriv-runner (rest coeffs) (add1 p)))]))
          ; Number Number -> Number
          ; Return coefficient of derivative of a*x^p
          (define (deriv-term a p)
            (* a p)))
    (deriv-runner (rest coeffs0) 1)))


;; [Listof Number] -> [Listof Number]
;; Returns indefinite integral as list of polynomial with input coeffs
;;   Constant "+C" is hardcoded as 0, the first coeff in output
(module+ test
  (check-equal? (poly-integral '(1)) '(0 1))
  (check-equal? (poly-integral '(1 2)) '(0 1 1))
  (check-equal? (poly-integral '(2 3 4)) '(0 2 3/2 4/3)))
(define (poly-integral coeffs0)
  (local ((define degree (length coeffs0))
          (define C 0)
          (define list-of-degrees (build-list degree values))
          ; Number Number -> Number
          ; Return coefficient of integral of a*x^p
          (define (integrate-term a p)
            (cond
              [(equal? p 0) a]
              [else (/ a p)])))
    (cons C
          (map (λ (x) (integrate-term (list-ref coeffs0 x) (add1 x)))
               list-of-degrees))))


;; [Listof Number] Number Number -> Number
;; Integrates polynomial from a to b
(module+ test
  (check-equal? (poly-def-integral '(2) 0 2) 4)
  (check-equal? (poly-def-integral '(0 1) 0 1) 1/2))
(define (poly-def-integral poly a b)
  (local ((define anti-deriv (poly-integral poly))
          (define anti-deriv@a (eval-poly anti-deriv a))
          (define anti-deriv@b (eval-poly anti-deriv b)))
    (- anti-deriv@b anti-deriv@a)))