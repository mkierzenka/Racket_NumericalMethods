#lang racket
(module+ test (require rackunit))

(provide replace  ;Replaces an element of a [Listof Number] with another number
         get-ith  ;Get's ith element of a [Listof Any]
         dot-prod-n ;Calculates dot product of 2 lists starting at an index
         dot-prod ;Calculates full dot product of 2 lists
         sub-list ;Returns inclusive sub-list of [Listof Any] from a to b
         m-dot-v  ;Dot product of Matrix and Vector
         my-make-list ;Make (list a,a+1,a+2...b)
         axpy         ;A*X+Y = [Listof Number]
         scale-vector);Multiplies each element in list by a constant



;; Matrix [Listof Number] -> [Listof Number]
;; Dotproduct of Matrix and a Vector
#;(module+ test
    (check-equal? (m-dot-v (list (list 1 2) (list 2 3)) (list 2 2))
                  (list 6 10))
    (check-equal? (m-dot-v (list (list 1 2 3) (list 4 5 6) (list 7 8 9))
                           (list 1 1 1))
                  (list 6 15 24)))
(define (m-dot-v m v)
  (map (λ (current-row) (dot-prod-n 1 current-row v))
       m))

;; [Listof Number] [Listof Number] -> Number
;; Returns sum of previos products btwn 2 lists in step
#;(module+ test
  (check-equal? (dot-prod (list 1 1) (list 5 1)) 6)
  (check-equal? (dot-prod (list 1 2 3) (list 1 2 3)) (+ 1 4 9))
  (check-equal? (dot-prod (list 4 5 6) (list 1 3 10)) 79))
(define (dot-prod c s)
  (dot-prod-n 1 c s))

;; Number [Listof Number] [Listof Number] -> Number
;; Returns sum of previous products btwn 2 lists in step, start at n (1-based)
#;(module+ test
    (check-equal? (dot-prod-n 2 (list 1 1) (list 5 1)) 1)
    (check-equal? (dot-prod-n 1 (list 1 2 3) (list 1 2 3)) (+ 1 4 9))
    (check-equal? (dot-prod-n 2 (list 4 5 6) (list 1 3 10)) 75))
(define (dot-prod-n n c s)
  (cond
    [(> n (length c)) 0]
    [else (+ (* (get-ith c n) (get-ith s n)) (dot-prod-n (add1 n) c s))]))


;; [Listof Number] Number Number -> [Listof Number]
;; Returns inclusive sublist from a to b
#;(module+ test
    (check-equal? (sub-list (list 1 2 3 4) 1 2) (list 1 2))
    (check-equal? (sub-list (list 1 2 3 4 5) 1 3) (list 1 2 3))
    (check-equal? (sub-list (list 1 2 3 4) 2 4) (list 2 3 4))
    (check-equal? (sub-list (list 2 2 2 2 3) 1 4) (list 2 2 2 2))
    (check-equal? (sub-list (list 4 6 5 3 1 2 3) 2 5) (list 6 5 3 1))
    (check-equal? (sub-list (list 1 2 3 4 5) -1 6) (list 1 2 3 4 5))
    (check-equal? (sub-list (list 1 2 3) 1 1) (list 1)))
(define (sub-list lon a b)
  (local ((define length-interval (add1 (- b a)))
          (define keep-from-a (last-k a lon))
          (define flip-k-f-a (reverse keep-from-a))
          (define new-list-length (length keep-from-a)))
    (reverse (last-k (add1 (- new-list-length length-interval)) flip-k-f-a))))


; Number [Listof Number] -> [Listof Number]
; Returns k+1 -> end elements of the list
; ASSUME 1 based
#;(module+ test
    (check-equal? (last-k 1 (list 1 2 3)) (list 1 2 3))
    (check-equal? (last-k 1 (list 1 2 3 4 5)) (list 1 2 3 4 5))
    (check-equal? (last-k 3 (list 5 4 3 2 1)) (list 3 2 1))
    (check-equal? (last-k 2 (list 1 2 3 4)) (list 2 3 4))
    (check-equal? (last-k 2 (list 1 2 3)) (list 2 3))
    (check-equal? (last-k 4 (list 1 2 3 4 5 6)) (list 4 5 6)))
(define (last-k k lon)
  (cond
    [(or (<= k 1) (empty? lon)) lon]
    [else (last-k (sub1 k) (rest lon))]))


;; Number Number [Listof Number] -> [Listof Number]
;; Places a into lon at index(1-based) then return, keeping rest of lon the same
#;(module+ test
    (check-equal? (replace 2 1 (list 1 2 3)) (list 2 2 3))
    (check-equal? (replace 4 2 (list 3 3 3 3 3)) (list 3 4 3 3 3))
    (check-equal? (replace 7 3 (list 8 8 8)) (list 8 8 7)))
(define (replace a index lon)
  (append (sub-list lon 1 (sub1 index))
          (list a)
          (sub-list lon (add1 index) (length lon))))


;; Number Number -> [Listof Number]
;; Returns a list of integers btwn a and b, inclusive
#;(module+ test
    (check-equal? (my-make-list 2 4) (list 2 3 4)))
(define (my-make-list a b)
  (map (λ (x) (+ (sub1 a) x))
       (build-list (add1 (- b a)) add1)))


;; [Listof Any] Number -> Any
;; Returns index-th element of list, 1 based index
#;(module+ test
    (check-equal? (get-ith (list 1 2 3 4) 2) 2)
    (check-equal? (get-ith (list 2 3 4 5) 1) 2))
(define (get-ith lon index)
  (first (sub-list lon index index)))


;; [Listof Number] Number [Listof Number] -> [Listof Number]
;; ax + y for each element y in second list
;; if one list is shorter than the other, fill with either ax or y
(module+ test
  (check-equal? (axpy (list 1 2 3) 5 (list 1 1 0)) (list 6 11 15))
  (check-equal? (axpy (list 1 2) 2 (list 2 2 3)) (list 4 6 3))
  (check-equal? (axpy (list 1 1 1) 5 (list 2 2)) (list 7 7 5))
  (check-equal? (axpy (list 1 1 1 1 1) 5 (list 2 2)) (list 7 7 5 5 5))
  (check-equal? (axpy (list 2 4 6 8) 2 '()) (list 4 8 12 16)))
(define (axpy lox0 a loy0)
  (local (; [Listof Number] [Listof Number] -> [Listof Number]
          ; Axpy of lox-c and loy-c
          (define (axpy-run lox-current loy-current)
            (cond
              [(and (empty? loy-current) (empty? lox-current)) '()]
              [(empty? loy-current)
               (cons (* a (first lox-current))
                     (axpy-run (rest lox-current) loy-current))]
              [(empty? lox-current)
               (cons (first loy-current)
                     (axpy-run lox-current (rest loy-current)))]
              [else (cons (+ (* (first lox-current) a) (first loy-current))
                          (axpy-run (rest lox-current) (rest loy-current)))])))
    (axpy-run lox0 loy0)))



;; Number [Listof Number] -> [Listof Number]
;; Multiplies each element in list by factor
(module+ test
  (check-equal? (scale-vector 1 '(1 1 1)) '(1 1 1))
  (check-equal? (scale-vector 2 '(1 2 3)) '(2 4 6))
  (check-equal? (scale-vector -1 '(2 3 4)) '(-2 -3 -4)))
(define (scale-vector factor coeffs)
  (axpy coeffs factor '()))