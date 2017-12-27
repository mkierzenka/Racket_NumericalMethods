#lang racket
(require plot
         "ODE_Vector_Solvers.rkt")
;;--This file contains a few Restricted Three Body Problems, and solutions
;;--   using adaptive and non-adaptive techniques


(plot-new-window? #t)

;;Problem 1
;(define mu (/ 1 82.45))
;(define 3body-x0 (list 1.2 0 0 -1.04935750983))
;(define 3body-period 6.192169331396)

;;Problem 2
;(define mu 0.012277471)
;(define 3body-x0 (list 0.994 0 0 -2.00158510637908252))
;(define 3body-period 17.0652165601579625)

;;Problem 3
(define mu 0.012277471)
(define 3body-x0 (list 0.994 0 0 -2.031732629557))
(define 3body-period 11.124340337266)



(define lam (- 1 mu))

#;(define numSteps 100000)
#;(define 3body-steps (build-list (add1 numSteps)
                                  (λ (x) (* x (/ 3body-period numSteps)))))

;; Number [Listof Number] -> [Listof Number]
;; Derivative function for 3-body problem
(define (deriv3body t X)
  (local ((define x1 (first X))
          (define x2 (second X))
          (define y1 (third X))
          (define y2 (fourth X))
          (define r1 (sqrt (+ (sqr (+ x1 mu)) (sqr y1))))
          (define r2 (sqrt (+ (sqr (- x1 lam)) (sqr y1))))
          ;(define r1 (sqrt (+ (sqr (+ x1 mu)) (sqr y1))))
          ;(define r2 (sqrt (+ (sqr (- x1 lam)) (sqr y1))))
          (define x1-prime x2)
          (define x2-prime (- (+ (* 2 y2) x1)
                              (/ (* lam (+ x1 mu)) (expt r1 3))
                              (/ (* mu (- x1 lam)) (expt r2 3))))
          (define y1-prime y2)
          (define y2-prime (- (+ (* -2 x2) y1)
                              (/ (* lam y1) (expt r1 3))
                              (/ (* mu y1) (expt r2 3)))))
    (list x1-prime x2-prime y1-prime y2-prime)))



(define RESULT
  (adaptive-stepper deriv3body 3body-x0 0 3body-period 0.01 (expt 10 -3)))
(define pts (map (λ (TX) (second TX)) RESULT))
(length pts) ;;Number of steps
(define 3body-xs (map (λ (X) (first X)) pts))
(define 3body-ys (map (λ (X) (third X)) pts))

(plot (points (map vector 3body-xs 3body-ys)
              #:x-min -1.5
              #:x-max 1.5
              #:y-min -1.5
              #:y-max 1.5
              #:size 1))


;Non-Adaptive
#;(define RES (steppr-vec rk4-vec deriv3body 3body-x0 3body-steps))
#;(define 3body-x1s (map (λ (X) (first X)) RES))
#;(define 3body-y1s (map (λ (X) (third X)) RES))