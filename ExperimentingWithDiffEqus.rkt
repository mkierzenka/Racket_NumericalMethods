#lang racket
(require plot
         "ODE_Solvers.rkt")

;;--A few examples of ordinary differential equation solvers in action
;;--Each is plotted against the derivative's true solution

#;(define ex-ts (build-list 51 (λ (a) (* a (/ pi 50)))))
#;(define ex-xs (stepper heun-step (λ (t x) (cos t)) 0 ex-ts))

#;(define ex2-ts (build-list 6 (λ (a) (* a (/ pi 5)))))
#;(define ex2-xs (stepper heun-step (λ (t x) (exp t)) 1 ex2-ts))
;;compare to forward euler step and other step methods

#;(define ex3-ts (build-list 21 (λ (a) (* a (/ 1 20)))))
#;(define ex3-xs (stepper forward-euler-step (λ (t x) (* -100.0 x)) 1 ex3-ts))


#;(plot-new-window? #t)

#;(plot (list (axes)
            (points (map vector ex-ts ex-xs))
            (function sin 0 pi #:color "red" #:label "sin")))

#;(plot (list (axes)
            (points (map vector ex2-ts ex2-xs))
            (function exp 0 pi #:color "red" #:label "sin")))

#;(plot (list (axes)
            (points (map vector ex3-ts ex3-xs))
            (function exp 0 pi #:color "red" #:label "e^x")))