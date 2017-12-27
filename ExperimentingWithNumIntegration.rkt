#lang racket
(require "Numerical_Integration.rkt")
(module+ test (require rackunit))
(require test-engine/racket-tests)

(define EPSILON (expt 10 -15))
  

;;--Exploring the difference between computing using numerical integration
;;--  and closed form solutions. Looking for example of error propagation

;;--Considerthe integral from 0 to 1 of (x^j)*(e^(x-1)), with j from 1 to 25
;;--Compare the results-integrated and (calc-iterations 25), do you see any
;;--  evidence of errors?

;; Numerical Integration
(define results-integrated
  (map (λ (j) (fnInt (λ (x) (* (exp (- x 1)) (expt x j))) 0 1
                     '(-1 0 1)
                     '(1/6 4/6 1/6)))
       (build-list 25 add1)))


;; Number -> [Listof Number]
;; Calculates the recursion (numerically) from v_0 to v_max-j, as a list of nums
;; v_0 = 1 - (1/e)
;; v_j = 1 - (j * v_(j-1))
(module+ test
  (check-equal? (calc-iterations 1) (list (/ 1 (exp 1))))
  (check-equal? (calc-iterations 2) (list (/ 1 (exp 1))
                                          (- 1 (* 2 (/ 1 (exp 1)))))))
(define (calc-iterations max-j)
  (local ((define v0 (- 1 (/ 1 (exp 1))))
          ; Number Number -> Number
          ; Calculates the list of v's from j and v_(j-1) up to v_max-j
          (define (calc-i-runner j prev-v)
            (cond
              [(> j max-j) '()]
              [else (local ((define current-val (- 1 (* j prev-v))))
                      (cons current-val
                            (calc-i-runner (add1 j) current-val)))])))
    (calc-i-runner 1 v0)))