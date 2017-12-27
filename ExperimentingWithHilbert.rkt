#lang racket
(module+ test (require rackunit))

(require "Matrix_List_Helpers.rkt"
         "Linear_Systems.rkt")
;;--The Hilbert Matrix is constructed using a simple rule, but is often used to
;;--  observe the accuracy of linear system solvers.
;;-- Run check-hilbert to see how numerical methods introduce error
;;--   A perfect world would mean the result of this function is a list of zeros

;; Number -> Matrix
;; Creates Hilbert Matrix for N dimension (square)
(module+ test
  (check-equal? (make-hilbert 2) (list (list 1 1/2)
                                       (list 1/2 1/3))))
(define (make-hilbert n)
  (local (;; Number -> [Listof Number]
          ;; Creates the kth row of hilbert matrix, length n
          (define (make-hilbert-row k)
            (map (λ (i) (/ 1 (+ k i))) (build-list n add1))))
    (map make-hilbert-row (build-list n values))))

;; Number -> Matrix
;; Creates Hilbert Matrix for N dimension (square) with decimals
(module+ test
  (check-equal? (make-hilbert-d 2) (list (list 1.0 0.5)
                                         (list 0.5 (exact->inexact 1/3)))))
(define (make-hilbert-d n)
  (local (;; Number -> [Listof Number]
          ;; Creates the kth row of hilbert matrix, length n, decimals
          (define (make-hilbert-row k)
            (map (λ (i) (exact->inexact (/ 1 (+ k i)))) (build-list n add1))))
    (map make-hilbert-row (build-list n values))))


;; Number -> [Listof Matrix [Listof Number]]
;; Makes Hilbert matrix and solution of given nxn size
(define (make-hilbert-pair n)
  (local ((define h (make-hilbert n))
          (define list-of-1 (map (λ (x) 1)
                                 (my-make-list 1 n)))
          (define rhs (m-dot-v h list-of-1)))
    (list h rhs)))

;; Number -> [Listof Matrix [Listof Number]]
;; Makes Hilbert matrix and solution of given nxn size, inexact
(define (make-hilbert-pair-d n)
  (local ((define h (make-hilbert-d n))
          (define list-of-1 (map (λ (x) 1)
                                 (my-make-list 1 n)))
          (define rhs (m-dot-v h list-of-1)))
    (list h rhs)))

;; Number -> [Listof Number]
;; Displays solution to Hilbert Matrix and rhs, decimals
(define (check-hilbert n)
  (local ((define h-pair (make-hilbert-pair-d n))
          (define h (first h-pair))
          (define rhs (second h-pair)))
    (axpy rhs -1 (m-dot-v h (linsys-solver h rhs)))))