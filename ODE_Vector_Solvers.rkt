#lang racket
(module+ test (require rackunit))
(require "Matrix_List_Helpers.rkt")

(provide adaptive-stepper
         adaptive-step
         rkf45
         stepper-vec
         rk4-vec)
;;--Functions for numerically solving Ordinary Differential Equations
;;--Initial Value Problems, only
;;--As opposed to ODE_Solvers.rkt, this file considers Vector Functions
;;--Adaptive Solvers of this form included at the end
;;--Derivative Function is: [Number [Listof Number]] -> [Listof Number]




;; StepFunction Derivative [Listof Number] [Listof Number]
;;                                               -> [Listof [Listof Number]]

;; StepFunction is [Derivative Number [Listof Number] Number -> [Listof Number]]
;; Computes list of x's visited while traveling along list of t's
;;   starting from initial, using single-step function
(module+ test
  (check-equal? (stepper-vec forward-euler-step-vec (λ (t x) '(0)) '(0) '(0 1))
                '((0) (0)))
  (check-equal? (stepper-vec forward-euler-step-vec (λ (t x) '(3)) '(1) '(1 2 3))
                '((1) (4) (7))))
(define (stepper-vec single-step deriv x0 lot-all)
  (local (;; [Listof Number] [Listof Number] -> [Listof [Listof Number]]
          ;; Uses current list of t's and previous x to make a list of x's
          (define (stepper-r x-prev lot)
            (cond
              [(empty? (rest lot)) '()]
              [else
               (local ((define new-x (single-step deriv (first lot)
                                                  x-prev (second lot))))
                 (cons new-x (stepper-r new-x (rest lot))))])))
    (cons x0 (stepper-r x0 lot-all))))


;;---------------------------Single Step Functions------------------------------
;; Derivative Number [Listof Number] Number -> [Listof Number]
;; Takes one step along derivative starting at t0,x0 to t1 returns x1
;;   x0 and x1 are [Listof Number]
(module+ test
  (check-equal? (forward-euler-step-vec (λ (t x) '(1 0)) 0 '(0 0) 1) '(1 0))
  (check-equal? (forward-euler-step-vec (λ (t x) '(3)) 0 (list 2) 1) '(5))
  (check-equal? (forward-euler-step-vec (λ (t x) '(-5)) 1 (list -5) 3)
                '(-15)))
(define (forward-euler-step-vec deriv t0 x0 t1)
  (axpy (deriv t0 x0) (- t1 t0) x0))

;; Derivative Number [Listof Number] Number -> [Listof Number]
;; Takes half step, then uses that derivative value to get x1
(module+ test
  (check-equal? (modified-euler-vec (λ (t x) '(0)) 0 '(0) 1) '(0))
  (check-equal? (modified-euler-vec (λ (t x) '(3)) 0 '(2) 1) '(5))
  (check-equal? (modified-euler-vec (λ (t x) '(-5)) 1 '(-5) 3) '(-15)))
(define (modified-euler-vec deriv t0 x0 t1)
  (local ((define mid-t (/ (+ t0 t1) 2))
          (define half-step (forward-euler-step-vec deriv t0 x0 mid-t))
          (define new-deriv-val (deriv mid-t half-step)))
    (axpy new-deriv-val (- t1 t0) x0)))


;; Derivative Number [Listof Number] Number -> [Listof Number]
;; Heun's Method, average derivative (using forward euler vec) at t0,x0 and t1
;;   then use to calculate better x1
(module+ test
  (check-equal? (heun-step-vec (λ (t x) '(0)) 0 '(0) 1) '(0))
  (check-equal? (heun-step-vec (λ (t x) '(3)) 0 '(2) 1) '(5))
  (check-equal? (heun-step-vec (λ (t x) '(-5)) 1 '(-5) 3) '(-15)))
(define (heun-step-vec deriv t0 x0 t1)
  (local ((define x1-try (forward-euler-step-vec deriv t0 x0 t1))
          (define deriv-start (deriv t0 x0))
          (define deriv-end (deriv t1 x1-try))
          (define avg-deriv (map (λ (s) (/ s 2))
                                 (axpy deriv-start 1 deriv-end))))
    (axpy avg-deriv (- t1 t0) x0)))


;; Derivative Number [Listof Number] Number -> [Listof Number]
;; Uses 4th-order Runge-Kutta method to find new x
(module+ test
  (check-equal? (rk4-vec (λ (t x) '(1 0)) 0 '(0 0) 1) '(1 0))
  (check-equal? (rk4-vec (λ (t x) '(3)) 0 (list 2) 1) '(5))
  (check-equal? (rk4-vec (λ (t x) '(-5)) 1 (list -5) 3) '(-15)))
(define (rk4-vec deriv t0 x0 t1)
  (local ((define h (- t1 t0))
          (define k1 (deriv t0 x0))     
          (define t-mid (/ (+ t0 t1) 2))
          (define k2 (deriv t-mid (axpy k1 (/ h 2) x0)))
          (define k3 (deriv t-mid (axpy k2 (/ h 2) x0)))
          (define k4 (deriv t1 (axpy k3 h x0)))
          (define x1-1 (axpy k2 2 k1))
          (define x1-2 (axpy k3 2 x1-1))
          (define x1-3 (axpy k4 1 x1-2)))
    (axpy x1-3 (/ h 6) x0)))



;;----------------------------------Examples------------------------------------
#;(define ex-ts (build-list 51 (λ (a) (* a (/ pi 50)))))
#;(define ex-xs (map (λ (x) (first x))
                     (stepper heun-step-vec (λ (t x) (list (cos t)))
                              '(0)
                              ex-ts)))

#;(define ex2-ts (build-list 6 (λ (a) (* a (/ pi 5)))))
#;(define ex2-xs (map (λ (x) (first x))
                      (stepper heun-step-vec (λ (t x) (list (exp t)))
                               '(1)
                               ex2-ts)))
;;compare to forward euler step and other step methods

#;(plot (list (axes)
              (points (map vector ex-ts ex-xs))
              (function sin 0 pi #:color "red" #:label "sin")))

#;(plot (list (axes)
              (points (map vector ex2-ts ex2-xs))
              (function exp 0 pi #:color "red" #:label "e^x")))


;;--------------------------------Adaptive Stepper------------------------------
;; Stepper Derivative [Listof Number] Number Number Number Number
;;                                  -> [Listof (list Number [Listof Number])]
;; Computes list of t's and matching x's visited traveling from t0 to t-final,
;;   using adaptive step sizes, with initial step length h
(define (adaptive-stepper a-step deriv x0 t0 t-final h tol)
  (local ((define t-btwn (+ t0 h))
          (define x&error (a-step deriv t0 x0 t-btwn))
          (define x1 (first x&error))
          (define err (argmax (λ (x) x) (second x&error))))
    (cond
      [(> t0 t-final) '()]
      [(> err tol) (adaptive-stepper a-step deriv x0 t0 t-final (/ h 2) tol)]
      [(< err (/ tol 4))
       (cons (list t-btwn x1)
             (adaptive-stepper a-step deriv x1 t-btwn t-final (* 2 h) tol))]
      [else (cons (list t-btwn x1)
                  (adaptive-stepper a-step deriv x1 t-btwn t-final h tol))])))


;; Derivative Number [Listof Number] Number -> [Listof [Listof Number] Number]
;; Heuns Method, average derivative (using euler 1step) at t0,x0 and t1,x1
;;   then use to calculate better x1.
;; Also returns error estimate (compared to forward euler)
(module+ test
  (check-equal? (adaptive-step (λ (t x) '(0)) 0 '(0) 1) '((0) (0)))
  (check-equal? (adaptive-step (λ (t x) '(3)) 0 '(2) 1) '((5) (0)))
  (check-equal? (adaptive-step (λ (t x) '(-5)) 1 '(-5) 3) '((-15) (0))))
(define (adaptive-step deriv t0 x0 t1)
  (local ((define x1-try (forward-euler-step-vec deriv t0 x0 t1))
          (define deriv-start (deriv t0 x0))
          (define deriv-end (deriv t1 x1-try))
          (define avg-deriv (map (λ (s) (/ s 2))
                                 (axpy deriv-start 1 deriv-end)))
          (define x1 (axpy avg-deriv (- t1 t0) x0))
          (define error-estimate (map (λ (x) (abs x)) (axpy x1 -1 x1-try))))
    (list x1 error-estimate)))

;; Derivative Number [Listof Number] Number -> [Listof [Listof Number] Number]
;; Runge-Kutta-Fehlberg 4 and 5, adaptive stepper (rk45)
;; Source:
;;http://mathfaculty.fullerton.edu/mathews//n2003/rungekuttafehlberg/RungeKuttaFehlbergProof.pdf
;;Note, there is a mistake in this text. Correct coeff is 2197/4104 not 2197/4101
(module+ test
  (check-equal? (rkf45 (λ (t x) '(0)) 0 '(0) 1) '((0) (0)))
  (check-equal? (rkf45 (λ (t x) (list (+ (sqr t) (* 2 t) t 5))) 0 '(2) 1) '((53/6) (0)))
  (check-equal? (rkf45 (λ (t x) (list (+ (* 4 (expt t 3)) (- (sqr t)) (* 2 t) 2))) 1 '(-5) 3) '((235/3) (0))))
(define (rkf45 deriv t0 x0 t1)
  (local ((define h (- t1 t0))
          (define k1 (scale-vector h (deriv t0 x0)))
          (define k2 (scale-vector h (deriv (+ t0 (/ h 4)) (axpy k1 1/4 x0))))
          (define k3 (scale-vector h (deriv (+ t0 (* 3/8 h))
                                            (axpy k2 9/32 (axpy k1 3/32 x0)))))
          (define k4 (scale-vector h (deriv (+ t0 (* 12/13 h))
                                            (axpy k3 7296/2197 (axpy k2 -7200/2197 (axpy k1 1932/2197 x0))))))
          (define k5 (scale-vector h (deriv (+ t0 h)
                                            (axpy k4 -845/4104 (axpy k3 3680/513 (axpy k2 -8 (axpy k1 439/216 x0)))))))
          (define k6 (scale-vector h (deriv (+ t0 (* 1/2 h))
                                            (axpy k5 -11/40 (axpy k4 1859/4104 (axpy k3 -3544/2565 (axpy k2 2 (axpy k1 -8/27 x0))))))))
          (define rk4-guess (axpy k5 -1/5 (axpy k4 2197/4104 (axpy k3 1408/2565 (axpy k1 25/216 x0)))))
          (define rk5-guess (axpy k6 2/55 (axpy k5 -9/50 (axpy k4 28561/56430 (axpy k3 6656/12825 (axpy k1 16/135 x0))))))
          (define error (map (λ (x) (abs x)) (axpy rk4-guess -1 rk5-guess))))
    (list rk5-guess error)))


;;---------------------------------Examples-------------------------------------
#;(define res (adaptive-stepper adaptive-step (λ (t x) (list (cos t))) '(0) 0 (* 2 pi) 1
                                (expt 10 -2)))
#;(define ts (map (λ (X) (first X)) res))
#;(define xs (map (λ (X) (first (second X))) res))

#;(plot (list (axes)
              (function sin 0 (* 2 pi))
              (points (map vector ts xs))))
;; The points should closely follow the actual function
