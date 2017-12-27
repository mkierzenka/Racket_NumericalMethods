#lang racket
(module+ test (require rackunit))
(require "Matrix_List_Helpers.rkt")
(provide forward-euler-step
         modified-euler
         heun-step
         stepper)

;;--Contains a variety of numerical solvers for Ordinary Differential Equations
;;--Note, these are for initial value problems
;;--Both single-step methods and a "stepper" function are included and provided

;;--Derivative Function is: [Number Number -> Number]




;; StepFunction Derivative Number [Listof Number] -> [Listof Number]

;; StepFunction is [Derivative Number Number Number -> Number]
;; Computes list of x's visited while traveling along list of t's
;;   starting from initial, using single-step function
(module+ test
  (check-equal? (stepper forward-euler-step (λ (t x) 0) 0 '(0 1)) '(0 0))
  (check-equal? (stepper forward-euler-step (λ (t x) 3) 1 '(1 2 3)) '(1 4 7)))
(define (stepper single-step deriv x0 lot-all)
  (local (;; Number [Listof Number] -> [Listof Number]
          ;; Uses current list of t's and previous x to make a list of x's
          (define (stepper-r x-prev lot)
            (cond
              [(empty? (rest lot)) '()]
              [else (local ((define new-x (single-step deriv
                                                       (first lot)
                                                       x-prev
                                                       (second lot))))
                      (cons new-x (stepper-r new-x (rest lot))))])))
    (cons x0 (stepper-r x0 lot-all))))

;;-----------------------Single Step Functions----------------------------------
;; Derivative Number Number Number -> Number
;; Takes one step along derivative starting at t0,x0 to t1 returns x1
(module+ test
  (check-equal? (forward-euler-step (λ (t x) 0) 0 0 1) 0)
  (check-equal? (forward-euler-step (λ (t x) 3) 0 2 1) 5)
  (check-equal? (forward-euler-step (λ (t x) -5) 1 -5 3) -15))
(define (forward-euler-step deriv t0 x0 t1)
  (+ x0 (* (- t1 t0) (deriv t0 x0))))


;; Derivative Number Number Number -> Number
;; Takes half step, then uses that derivative value to get x1
(module+ test
  (check-equal? (modified-euler (λ (t x) 0) 0 0 1) 0)
  (check-equal? (modified-euler (λ (t x) 3) 0 2 1) 5)
  (check-equal? (modified-euler (λ (t x) -5) 1 -5 3) -15))
(define (modified-euler deriv t0 x0 t1)
  (local ((define mid-t (/ (+ t0 t1) 2))
          (define half-step (forward-euler-step deriv t0 x0 mid-t))
          (define new-deriv-val (deriv mid-t half-step)))
    (+ x0 (* (- t1 t0) new-deriv-val))))


;; Derivative Number Number Number -> Number
;; Heun's Method, average derivative (using forward euler) at t0,x0 and t1
;;   then use to calculate better x1
(module+ test
  (check-equal? (heun-step (λ (t x) 0) 0 0 1) 0)
  (check-equal? (heun-step (λ (t x) 3) 0 2 1) 5)
  (check-equal? (heun-step (λ (t x) -5) 1 -5 3) -15))
(define (heun-step deriv t0 x0 t1)
  (local ((define x1-try (forward-euler-step deriv t0 x0 t1))
          (define deriv-start (deriv t0 x0))
          (define deriv-end (deriv t1 x1-try))
          (define better-deriv (/ (+ deriv-start deriv-end) 2)))
    (+ x0 (* (- t1 t0) better-deriv))))