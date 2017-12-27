#lang racket
(module+ test (require rackunit))
(require test-engine/racket-tests)

;; [Number -> Number] Number Number Number -> Number
;; Finds roots of input function f via bisection on [a0,b0] to within tolerance
;; Error if (f a0) and (f b0) have same sign
(module+ test
  (check-within (bisect (λ (x) x) -3 3 0.001) 0 0.001)
  (check-within (bisect (λ (x) (sin x)) 2 4 0.001) pi 0.001))

(define (bisect f a0 b0 tol)
  (local (; Number Number -> Number
          ; If range between a and b is within tol, return a
          ; else recur on reduced range.
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
(check-expect (my-sgn 3) 1)
(check-expect (my-sgn 0) 1)
(check-expect (my-sgn -3) -1)

(define (my-sgn x)
  (cond
    [(negative? x) -1]
    [else 1]))