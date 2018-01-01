#lang racket
(module+ test (require rackunit))
(require test-engine/racket-tests)

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
(define m-ut-ex2 '((2  1)
                   (0 1/2)))
(define rhs-ex2 (list 2 3))
(define rhs-ut-ex2 (list 3 1/2))
(define sol-ex2 (list 1 1))

(define m-ex3 '((1 2 3)
                (4 5 6)
                (7 8 9)))
(define m-ut-ex3 '((7 8 9)
                   (0 6/7 12/7)
                   (0  0  0)))
(define rhs-ex3 (list 2 2 2))
(define rhs-ut-ex3 (list 2 12/7 0))

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
;; Uses Partial pivoting and Gaussian Elimination
(module+ test
  (check-equal? (linsys-solver m-ex2 rhs-ex2) sol-ex2)
  (check-equal? (linsys-solver m-ex4 rhs-ex4) sol-ex4)
  (check-equal? (linsys-solver '((0 2 0) (2 0 0) (0 0 2)) '(4 4 4)) (list 2 2 2))
  (check-equal? (linsys-solver '((0 0 3) (0 2 0) (2 0 0)) '(4 4 4)) (list 2 2 4/3))
  (check-error (linsys-solver '((0 0 3) (0 0 0) (2 0 0)) '(4 4 4)))) ;Singular Matrix
(define (linsys-solver A rhs)
  (local ((define a&rhs (combine A rhs))
          (define new-a&rhs (elim-ab a&rhs))
          (define ut-A (get-M new-a&rhs))
          (define ut-rhs (get-rhs new-a&rhs)))
    (backward-sub ut-A ut-rhs)))

  
;; Matrix [Listof Number] -> [Listof Number]
;; Solves the system of linear equations using backward substitution:
;; A*<x1,x2...xn> = <c1,c2...cn>   -> A is coefficient matrix, <c1...cn> is rhs
;; ASSUME Matrix is upper-triangular and all diagonal entries are non-zero
(module+ test
  (check-equal? (backward-sub m-ut-ex1 rhs-ut-ex1) ;x+y=5, y=3
                sol-ex1)
  (check-equal? (backward-sub (list (list 1 1) (list 0 3)) (list 5 1))
                (list 14/3 1/3))
  (check-equal? (backward-sub (list (list 2 1) (list 0 -1)) (list 3 -1))
                (list 1 1)))

(define (backward-sub A rhs0)
  (local (; Number [Listof Number] -> [Listof Number]
          ; Uses current righthand side to solve the kth element of row
          (define (b-sub-r k rhs)
            (replace (b-sub-row (get-ith A k) rhs k) k rhs)))
    (foldr (λ (row-index current-rhs) (b-sub-r row-index current-rhs))
           rhs0
           (my-make-list 1 (length A)))))


;; [Listof Number] [Listof Number] Number -> Number
;; Uses row of A and col of constants/solution to get kth solution
;; Assume row[k] isn't 0
(module+ test
  (check-equal? (b-sub-row (list 1 3) (list 1 3) 1) -8)
  (check-equal? (b-sub-row (list 1 1) (list 5 1) 1) 4))
(define (b-sub-row coefs sols k)
  (/ (- (get-ith sols k) (dot-prod-n (add1 k) coefs sols))
     (get-ith coefs k)))


;-------------------------------Elimination-------------------------------------
;; Ext-Matrix -> Ext-Matrix
;; Make A and rhs (in one) into upper triangular form, uses partial pivoting
(module+ test
  (check-equal? (elim-ab (combine m-ex2 rhs-ex2))
                (combine m-ut-ex2 rhs-ut-ex2))
  (check-equal? (elim-ab (combine m-ex3 rhs-ex3))
                (combine m-ut-ex3 rhs-ut-ex3)))
(define (elim-ab a&rhs0)
  (local (;; Number Ext-Matrix -> Ext-Matrix
          ;; Eliminates col column and remaining cols in Ext-Matrix (not rhs)
          ;; Uses partial pivoting (through split func)
          #;(module+ test
              (check-equal? (elim-col&rhs 1 (list (list 1 1 2) (list 2 1 3)))
                            (list (list 2 1 3) (list 0 1/2 1/2))))
          (define (elim-col&rhs col a&rhs)
            (cond
              [(empty? a&rhs) '()]
              [else (local ((define split (my-split a&rhs (λ (row) (abs (get-ith row col)))))
                            (define r (first split))
                            (define new-a&rhs (new-elim col r (second split))))
                      (cons r (elim-col&rhs (add1 col) new-a&rhs)))])))
    (elim-col&rhs 1 a&rhs0)))


;; Number [Listof Number] Ext-Matrix -> Ext-Matrix
;; Eliminate kth column from rest of Matrix using the given row
;; Throws errors if Matrix is singular
(module+ test
  (check-equal? (new-elim 1 (list 1 1 2) (list (list 2 1 3)))
                (list (list 0 -1 -1)))
  (check-equal? (new-elim 2 (list 1 1 2) (list (list 2 1 3)))
                (list (list 1 0 1)))
  (check-equal? (new-elim 1 (list 1 1 2) (list (list 2 4 3) (list 4 1 1)))
                (list (list 0 2 -1) (list 0 -3 -7)))
  (check-error (new-elim 1 '(0 0 3) '((1 0 0) (2 0 0))))
  (check-error (new-elim 1 '(0 0 3) '((0 2 0) (2 0 0)))))
(define (new-elim col row rest-a&rhs)
  (local ((define A (get-M rest-a&rhs))
          (define N (length A))
          (define Akk (get-ith row col)))
    (map (λ (x)
           (local ((define coef (safe-coef (get-ith (get-ith rest-a&rhs x) col) Akk)))
             (axpy row (- coef) (get-ith rest-a&rhs x))))
         (my-make-list 1 N))))


;; [Listof Any] [Any -> Number] -> (list Any [Listof Any])
;; Returns the Any which maximizes the score function and the rest of the list
(module+ test
  (check-equal? (my-split '(1 2 3 2 1) (λ (x) (abs x))) (list 3 '(1 2 2 1)))
  (check-equal? (my-split '(1 2 3 3 2 3 1) (λ (x) (abs x))) (list 3 '(1 2 3 2 3 1)))
  (check-equal? (my-split '(1 2 -3 -2 1) (λ (x) (abs x))) (list -3 '(1 2 -2 1))))
(define (my-split loa scorer)
  (local ((define best-val (argmax scorer loa))
          (define others (remove best-val loa)))
    (list best-val others)))

;; Matrix [Listof Number] -> Ext-Matrix
;; Combines a Matrix and rhs into one extended matrix
;; ASSUME Number of rows in Matrix = Number of elements in list
(module+ test
  (check-equal? (combine (list (list 1 2) (list 4 5)) (list 3 6))
                (list (list 1 2 3) (list 4 5 6)))
  (check-equal? (combine (list (list 3 5 3) (list 1 1 1) (list 2 2 2)) (list 8 8 8))
                (list (list 3 5 3 8) (list 1 1 1 8) (list 2 2 2 8))))
(define (combine A rhs)
  (local ((define N (length A))
          (define indices (my-make-list 1 N)))
    (map (λ (k) (append (get-ith A k) (list (get-ith rhs k)))) indices)))

;; Ext-Matrix -> Matrix
;; Extracts Matrix from the combined Matrix and rhs
(module+ test
  (check-equal? (get-M (combine (list (list 1 2) (list 4 5)) (list 3 6)))
                (list (list 1 2) (list 4 5)))
  (check-equal? (get-M (combine (list (list 3 5 3) (list 1 1 1) (list 2 2 2)) (list 8 8 8)))
                (list (list 3 5 3) (list 1 1 1) (list 2 2 2))))
(define (get-M m&rhs)
  (local ((define m-length (length m&rhs)))
    (map (λ (row) (take row m-length)) m&rhs)))

;; Ext-Matrix -> [Listof Number]
;; Extracts rhs from the combined Matrix and rhs
(module+ test
  (check-equal? (get-rhs (combine (list (list 1 2) (list 4 5)) (list 3 6))) (list 3 6))
  (check-equal? (get-rhs (combine (list (list 3 5 3) (list 1 1 1)) (list 8 8))) (list 8 8)))
(define (get-rhs m&rhs)
  (map (λ (row) (last row)) m&rhs))

;; Number Number -> Number
;; Performs division to calculate coefficient or throws proper error
(module+ test
  (check-equal? (safe-coef 0 5) 0)
  (check-equal? (safe-coef 1 2) 1/2)
  (check-error (safe-coef 2 0))
  (check-error (safe-coef 0 0)))
(define (safe-coef elem Akk)
  (cond
    [(and (zero? elem) (zero? Akk)) (error "Singular Matrix, Infinitely many solutions")]
    [(zero? Akk) (error "Singular Matrix, No solutions")]
    [else (/ elem Akk)]))