#lang racket

(require plot
         "Matrix_List_Helpers.rkt"
         "Linear_Systems.rkt"
         "Polynomial_Ops.rkt"
         "Polynomial_Interpolation.rkt")


;;--This file experiments with Polynomial Interpolation through an example
;; Look how the polynomial fitted through the points changes depending on which
;;   points are used in the approximation

(define POINTS (list '(5 1.226) '(10 1.498) '(15 1.822) '(20 2.138) '(30 2.662)
                     '(40 2.840) '(50 2.807) '(60 2.542) '(70 2.210) '(80 1.877)
                     '(90 1.539) '(100 1.201)))
(define 6-POINTS (map (λ (i) (list-ref POINTS i)) (list 1 3 5 7 9 11)))
(define 6PT-x (map (λ (pt) (first pt)) 6-POINTS))
(define 6PT-y (map (λ (pt) (second pt)) 6-POINTS))
(define Others (filter (λ (pt) (not (member pt 6-POINTS))) POINTS))
(define Others-x (map (λ (pt) (first pt)) Others))
(define Others-y (map (λ (pt) (second pt)) Others))


(define All-x (map (λ (x) (first x)) POINTS))
(define All-y (map (λ (x) (second x)) POINTS))

(define P6 (calc-poly 6-POINTS))
(define P6-func (λ (x) (eval-poly P6 x)))

(define P12 (calc-poly POINTS))
(define P12-func (λ (x) (eval-poly P12 x)))

(define Others-approx (map P6-func Others-x))
(define error (axpy Others-y -1 Others-approx))

(plot-new-window? #true)

(plot (list (axes)
            (function P6-func 0 100)
            (points (map vector Others-x Others-approx)
                    #:size 2)
            (points (map vector All-x All-y))
            (function P12-func 0 100
                      #:color "red")))
;;The Plot:
;;- The polynomial fitted to the 6 points
;;- The calculated values for the other 6 points
;;- The polynomial fitted to 12 points (red)

;; Look at the accuracy of of the points we calculated vs. the real values
;; Also, look how the polynomial changes shape from 6-pts to 12-pts
;; Note that it become "curvier", that is, has more peaks as the number of
;;   points considered in the approximation increases.