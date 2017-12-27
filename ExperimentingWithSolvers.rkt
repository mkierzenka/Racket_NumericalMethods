#lang racket
#;(require plot)
(require "Solvers.rkt")


;;-- This file has a few examples to compare equation solvers (finding roots)
;;--   and their accuracy. Change the particular solvers used to compare results.
;;-- Solvers to experiment with: my-solver, my-solver-newton, my-solver-sec


;; A particular example for the equations below
(define KEPLER (λ (x e r) (- (- x (* e (sin x))) r)))
(define KEPLER-prime (λ (x e) (- 1 (* e (cos x)))))

#;(plot (list (axes)
              (function (λ (r) (my-solver-newton
                                (λ (x) (KEPLER x 0.04 r))
                                (λ (x) (KEPLER-prime x 0.04)) -100 0.000001))
                        0
                        (* 3 pi)
                        #:color "red"
                        #:label "e = 0.04")
              (function (λ (r) (my-solver-newton
                                (λ (x) (KEPLER x 0.5 r))
                                (λ (x) (KEPLER-prime x 0.5)) -100 0.000001))
                        0
                        (* 3 pi)
                        #:color "blue"
                        #:label "e = 0.5")
              (function (λ (r) (my-solver-newton
                                (λ (x) (KEPLER x 0.96 r))
                                (λ (x) (KEPLER-prime x 0.96)) -10 0.000001))
                        0
                        (* 3 pi)
                        #:color "green"
                        #:label "e = 0.96")))


;(define list-of-r (map (λ (x) (* (* 3/100 pi) x)) (build-list 100 add1)))
;(define outputs (map (λ (r) (my-solver (λ (x) (KEPLER x 0.5 r)) -100 100 0.000001)) list-of-r))

;;-- Solving for well-known functions from their derivatives
;;-- Notice that the derivative does not need to be perfect
;(my-solver-sec sin 3.5 4 (expt 10 -13))
;(my-solver-newton sin (λ (x) (cos x)) 3.5 (expt 10 -13))
;(my-solver-newton sin (λ (x) 1) 3.5 (expt 10 -13))
;(my-solver-newton atan (λ (x) pi) 10 (expt 10 -13))