#lang racket

(require "Continued_Fractions.rkt")


;;--This file includes a couple ways to calculate pi to some precision

;;--First is the Chudnovsky Algorithm, used as recently as 2016 to calculate
;;--  pi to world record lengths
;;--https://en.wikipedia.org/wiki/Chudnovsky_algorithm


;; Number -> Rational
;; Uses k-final terms of Chudnovsky Algorithm to calculate pi
(define (calc-pi k-final)
  (local (;Number -> Rational
          ;Calculates pi from k to k-final
          (define (calc-pi-runner k)
            (cond
              [(> k k-final) 0]
              [else
               (local ((define next (/ (* (factorial (* 6 k)) (+ (* 545140134 k) 13591409))
                                       (* (factorial (* 3 k)) (expt (factorial k) 3) (expt -262537412640768000 k)))))
                 (+ next (calc-pi-runner (add1 k))))])))
    (/ (* 426880 (sqrt-n-decimals 10005 50 100)) (calc-pi-runner 0))))


;; Number Number Number -> Rational
;; Returns the squareroot of a number as a rational, but the decimals expansion
;;   is accurate to n places
(define (sqrt-n-decimals a prev-approx n)
  (local ((define newton (λ (x) (- x (/ (- (sqr x) a) (* 2 x)))))
          (define new-approx (newton prev-approx))
          (define new-str (real->decimal-string new-approx n))
          (define prev-str (real->decimal-string prev-approx n)))
    (cond
      [(equal? new-str prev-str) new-approx]
      [else (sqrt-n-decimals a new-approx n)])))


;; Number -> Number
;; Returns n!
(define (factorial n)
  (cond
    [(equal? n 0) 1]
    [else (* n (factorial (sub1 n)))]))




;; Approximation Using continued fractions
;; Note, to truely have pi requires an infinite fraction
;; Thanks to the On-line Encyclopedia of Integer Sequences
;; http://oeis.org/A001203
(define pi-approx-cont-fraction
  (list 3 7 15 1 292 1 1 1 2 1 3 1 14 2 1 1 2 2 2 2 1 84 2 1 1 15 3 13 1 4 2 6 6
        99 1 2 2 6 3 5 1 1 6 8 1 7 1 2 3 7 1 2 1 1 12 1 1 1 3 1 1 8 1 1 2 1 6 1        
        1 5 2 2 3 1 2 4 4 16 1 161 45 1 22 1 2 2 1 4 1 2 24 1 2 1 3 1 2 1))

#;(get-rational pi-approx-cont-fraction)
#;(real->decimal-string (get-rational pi-approx-cont-fraction) 68)


#;(real->decimal-string (calc-pi 10) 300)
;;If this is not accurate enough for you, increase digits printed and increase
;;  accuracy of calculated square-root in the function








