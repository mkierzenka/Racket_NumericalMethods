#lang racket
(module+ test (require rackunit))

(provide get-frac
         get-rational)

;;-- Reference: Numerical Methods in Scientific Computing Vol. 1
;;--        By: Germund Dahlquist and Åke Björck

;;-- Continued fractions can be used to represent real numbers & have the form: 
;;                      1
;;  x =   a0 +  -----------------
;;                         1
;;               a1 + -----------
;;                            1
;;                     a2 + -----
;;                           ...

;;-- These can be infinite or finite, depending on your choice of x
;;-- Specifically, we will focus on positive numbers which can be represented by
;;-- x = b0 + 1    1    1   
;;           ---  ---  ---   ...
;;           b1+  b2+  b3+
;;-- The position of the addition symbol is shorthand for the layering shown
;;--   in the image

;;-- We will represent a sequence of such additions as a list of the b values
;;-- ContinuedFraction is thus a [Listof Number]


;;-- Sidenote, the built-in function (rationalize x tol) does the same thing as
;;--   my function and the outputs seem to almost always agree
;;-- But when you push the functions to double precision, they sometimes are
;;--   a little different. But Racket prioritizes smallest denominator, it seems
;; Try this: (check-rationalize (sin 765/738))


;; Number -> Number
;; Takes a number, calls my function (tol 10^-16) and built-in func to compare
;;   results
(define (check-rationalize x)
  (local ((define mine (get-rational (get-frac (inexact->exact x) (expt 10 -16))))
          (define rack (rationalize (inexact->exact x) (expt 10 -16))))
    (println mine)
    (println rack)
    (- mine rack)))

;; Number -> ContinuedFraction
;; Takes a number and returns a continued fraction that matches it to some tol
(module+ test
  (check-equal? (get-frac pi 0.01) '(3 7))
  (check-equal? (get-frac pi 0.0001) '(3 7 15))
  (check-equal? (get-frac 2 0.001) '(2)))
(define (get-frac x0 tol)
  (local ((define nmax 100) ;;x=rational => finite expansion, but could be big
          (define b0 (inexact->exact (floor x0)))
          (define p-1 1)
          (define p0 b0)
          (define q-1 0)
          (define q0 1)
          ;; (list Number Number Number Number) Number Number -> ContinuedFraction
          ;; Takes (pn-1,pn-2,qn-1,qn-2),current xn, n to calculate next level or stop
          ;; pn = pn-2 + bn*pn-1, qn = qn-2 + bn*qn-2, bn=(floor xn)
          (define (get-frac-runner old-stuff xn n)
            (local ((define pn-1 (first old-stuff))
                    (define pn-2 (second old-stuff))
                    (define qn-1 (third old-stuff))
                    (define qn-2 (fourth old-stuff))
                    (define bn (inexact->exact (floor xn)))
                    (define dif (- xn bn));These 2 lines are to remove 1/0 error
                    (define next-x (if (equal? dif 0) xn (/ 1 (- xn bn))))
                    (define pn (+ pn-2 (* bn pn-1)))
                    (define qn (+ qn-2 (* bn qn-1)))
                    (define new-stuff (list pn pn-1 qn qn-1)))
              (cond
                [(> n nmax) (error "Maximum Depth Reached")]
                [(< (abs (- x0 (/ pn qn))) tol) (cons bn '())]
                [else (cons bn (get-frac-runner new-stuff next-x (add1 n)))]))))
    (cond
      [(zero? (- x0 b0)) (list b0)]
      [else (local ((define x1 (/ 1 (- x0 b0))))
              (cons b0 (get-frac-runner (list p0 p-1 q0 q-1) x1 0)))])))


;; ContinuedFraction -> Rational
;; Evaluates the continued fraction to a fraction
(module+ test
  (check-equal? (get-rational '(3 7)) (+ 3 (/ 1 7)))
  (check-equal? (get-rational '(1 2 3)) (+ 1 (/ 1 (+ 2 (/ 1 3)))))
  (check-equal? (get-rational '(1 2 3 4)) 43/30)
  (check-equal? (get-rational '(3 7 15 1 292))
                (+ 3 (/ 1 (+ 7 (/ 1 (+ 15 (/ 1 (+ 1 (/ 1 292))))))))))
(define (get-rational lob)
  (cond
    [(empty? (rest lob)) (first lob)]
    [else (+ (first lob) (/ 1 (get-rational (rest lob))))]))