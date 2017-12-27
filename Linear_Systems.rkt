#lang racket
(module+ test (require rackunit))
(require "Matrix_List_Helpers.rkt")

(provide linsys-solver) ;Solves a system of linear equations

;; A Matrix is:
;; (list (list x0 x1 x2 ...xn)
;;       (list x0 x1 x2 ...xn)
;;       (list x0 x1 x2 ...xn) ...)
;; Interpretation: A Matrix represents a matrix where each sub-list is a row

;; An UTMatrix is a Matrix that is upper-triangular


(define m-ut-ex1 '((1 1)
                   (0 1)))
(define rhs-ut-ex1 (list 5 3))
(define sol-ex1 (list 2 3))

(define m-ex2 '((1 1)
                (2 1)))
(define m-ut-ex2 '((1  1)
                   (0 -1)))
(define rhs-ex2 (list 2 3))
(define rhs-ut-ex2 (list 2 -1))
(define sol-ex2 (list 1 1))

(define m-ex3 '((1 2 3)
                (4 5 6)
                (7 8 9)))
(define m-ut-ex3 '((1  2  3)
                   (0 -3 -6)
                   (0  0  0)))
(define rhs-ex3 (list 2 2 2))
(define rhs-ut-ex3 (list 2 -6 0))

(define m-ex4 '(( 1  2 1)
                ( 2  2 3)
                (-1 -3 0)))
(define rhs-ex4 (list 0 3 2))
(define sol-ex4 (list 1 -1 1))

(define filips-example (list (list 6 2)
                             (list 12 192)))
(define filips-rhs (list 1.5 500.1))


;; Matrix [Listof Number] -> [Listof Matrix [Listof Number]]
;; Solves a general system of linear equations A*<x1,x2...xn> = <c1,c1...cn>
(module+ test
  (check-equal? (linsys-solver m-ex2 rhs-ex2) sol-ex2)
  (check-equal? (linsys-solver m-ex4 rhs-ex4) sol-ex4))
(define (linsys-solver A rhs)
  (local ((define a&rhs (list A rhs))
          (define new-a&rhs (elim-ab a&rhs))
          (define ut-A (first new-a&rhs))
          (define ut-rhs (second new-a&rhs)))
    (linsys-solver-ut ut-A ut-rhs)))

  
;; Matrix [Listof Number] -> [Listof Number]
;; Solves the system of linear equations using backward substitution:
;; A*<x1,x2...xn> = <c1,c2...cn>   -> A is coefficient matrix, <c1...cn> is rhs
;; ASSUME Matrix is upper-triangular and all diagonal entries are non-zero
(module+ test
  (check-equal? (linsys-solver-ut m-ut-ex1 rhs-ut-ex1) ;x+y=5, y=3
                sol-ex1)
  (check-equal? (linsys-solver-ut (list (list 1 1) (list 0 3)) (list 5 1))
                (list 14/3 1/3))
  (check-equal? (linsys-solver-ut (list (list 2 1) (list 0 -1)) (list 3 -1))
                (list 1 1)))

(define (linsys-solver-ut A rhs0)
  (local (; Number [Listof Number] -> [Listof Number]
          ; Uses current righthand side to solve the kth element of row
          (define (linsys-solver-r k rhs)
            (replace (linsys-solver-sub (get-ith A k) rhs k) k rhs)))
    (foldr (λ (row-index current-rhs) (linsys-solver-r row-index current-rhs))
           rhs0
           (my-make-list 1 (length A)))))


;; [Listof Number] [Listof Number] Number -> Number
;; Uses row of A and col of constants/solution to get kth solution
;; Assume row[k] isn't 0
(module+ test
  (check-equal? (linsys-solver-sub (list 1 3) (list 1 3) 1) -8)
  (check-equal? (linsys-solver-sub (list 1 1) (list 5 1) 1) 4))
(define (linsys-solver-sub coefs sols k)
  (/ (- (get-ith sols k) (dot-prod-n (add1 k) coefs sols))
     (get-ith coefs k)))


;-------------------------------Elimination-------------------------------------
;; [Listof Matrix [Listof Number]] -> [Listof UTMatrix [Listof Number]]
;; Make A and rhs into upper triangular form
(module+ test
  (check-equal? (elim-ab (list m-ex2 rhs-ex2))
                (list m-ut-ex2 rhs-ut-ex2))
  (check-equal? (elim-ab (list m-ex3 rhs-ex3))
                (list m-ut-ex3 rhs-ut-ex3)))
(define (elim-ab a&rhs)
  (foldl elim-col&rhs
         a&rhs
         (my-make-list 1 (sub1 (length (first a&rhs))))))


;; Number [Listof Matrix [Listof Number]] -> [Listof Matrix [Listof Number]]
;; Makes A and rhs into uppertriangular
(module+ test
  (check-equal? (elim-col&rhs 1 (list (list (list 1 1) (list 2 1)) (list 2 3)))
                (list (list (list 1 1) (list 0 -1)) (list 2 -1))))
(define (elim-col&rhs col m&rhs)
  (local ((define A (first m&rhs))
          (define N (length A))
          (define rhs (second m&rhs))
          (define Akk (get-ith (get-ith A col) col))
          (define elim-A
            (map (λ (x)
                   (local ((define coef (/ (get-ith (get-ith A x) col) Akk)))
                     (axpy (get-ith A col) (- coef) (get-ith A x))))
                 (my-make-list (add1 col) N)))
          (define elim-rhs
            (map (λ (x)
                   (local ((define coef (/ (get-ith (get-ith A x) col) Akk)))
                     (+ (get-ith rhs x) (* (- coef) (get-ith rhs col)))))
                 (my-make-list (add1 col) N)))
          (define new-A (append (sub-list A 1 col) elim-A))
          (define new-B (append (sub-list rhs 1 col) elim-rhs)))
    (list new-A new-B)))