#lang racket
(module+ test (require rackunit))
(require test-engine/racket-tests)
(provide my-solver
         my-solver-sec
         my-solver-newton)

;;--Includes a variety of solvers to find the root of a function, note that each
;;--  requires particular inputs
;;--Bisection, Newton's Method, and Secant Method

;; [Number -> Number] Number Number Number -> Number
;; Finds roots of input function f via bisection on [a0,b0] with tol
;; Error if (f a0) (f b0) have same sign
(define (my-solver f a0 b0 tol)
  (local (; Number Number -> Number
          ; If range btwn a and b is within tol, return a
          ; else recur on reduced range
          (define (bisect-runner a b)
            (local ((define f@a (f a))
                    (define f@b (f b))
                    (define midpt (/ (+ a b) 2))
                    (define f@mid (f midpt)))
              (cond
                [(equal? (my-sgn f@a) (my-sgn f@b)) (error "Same sign")]
                [(<= (abs (- a b)) tol) a]
                [(equal? (my-sgn f@a) (my-sgn f@mid)) (bisect-runner midpt b)]
                [(equal? (my-sgn f@b) (my-sgn f@mid)) (bisect-runner a midpt)]))))
    (bisect-runner (min a0 b0) (max a0 b0))))


;; Number -> Number
;; Returns 1 if x is positive or 0, else returns -1
(define (my-sgn x)
  (cond
    [(negative? x) -1]
    [else 1]))


;; [Number -> Number] Number Number Number -> Number
;; Uses Secant method to find root of f with tolerance tol, starting with a0,b0
#;(module+ test
    (check-equal? (my-solver-sec (λ (x) x) -2 3 0.001) 0)
    (check-within (my-solver-sec sin 2 4 0.00001) pi 0.00001)
    (check-error (my-solver-sec (λ (x) (+ (sqr x) 1)) -2 2 0.0001)))
(define (my-solver-sec f a0 b0 tol)
  (local (; Number Number -> Number
          ; Find secant through A = a,(f a) and B = b,(f b)
          ; If distance btwn A and B is within tol, return a
          ; else recur on new point and closest of A or B
          ; Accumulator its: number of iterations so far
          (define (sec-runner a b its)
            (cond
              [(> its 50) (error "Maximum number of iterations (50) exceeded")]
              [(<= (abs (- a b)) tol) a]
              [else  (local ((define f@a (f a))
                             (define f@b (f b))
                             (define inv-slope (/ (- b a) (- f@b f@a)))
                             (define closer-root (- a (* f@a inv-slope))))
                       (sec-runner closer-root (closest a b closer-root) (add1 its)))])))
    (sec-runner a0 b0 0)))


;; Number Number Number -> Number
;; Returns the number (a or b) that is closest to target, or a if equidistant
#;(module+ test
    (check-equal? (closest 1 2 3) 2)
    (check-equal? (closest -2 10 0) -2)
    (check-equal? (closest 6 6 7) 6)
    (check-equal? (closest -2 3 0) -2)
    (check-equal? (closest 5 8 5) 5))
(define (closest a b target)
  (local ((define a-dist (abs (- a target)))
          (define b-dist (abs (- b target))))
    (cond
      [(< b-dist a-dist) b]
      [else a])))


;; [Number -> Number] [Number -> Number] Number Number -> Number
;; Uses Newton's method with f, f-prime, and initial guess to get root with tol
;; Limited to 100 iterations, since a solution is not guaranteed
(module+ test
  (check-equal? (my-solver-newton (λ (x) x) (λ (x) 1) 2 0.001) 0)
  (check-within (my-solver-newton sin cos 2 0.00001) pi 0.00001)
  (check-error (my-solver-newton (λ (x) (+ (sqr x) 1)) (λ (x) (* 2 x)) 2 0.0001)))
(define (my-solver-newton f f-prime r0 tol)
  (local (; Number -> Number
          ; Find derivative at r
          ; If f@r/f-prime@r is within tol, return r
          ; else return next approx using Newton's method
          ; Accumulator its: number of iterations so far, stop if too many
          (define (newt-runner r its)
            (local ((define f@r (f r))
                    (define f-prime@r (f-prime r)))
              (cond
                [(> its 100) (error "Number of iterations " its)]
                [(equal? f-prime@r 0) (error "Derivative at " r " is 0")]
                [(<= (abs (/ f@r f-prime@r)) tol) r]
                [else (local ((define closer-root (- r (/ f@r f-prime@r))))
                        (newt-runner closer-root (add1 its)))]))))
    (newt-runner r0 0)))


