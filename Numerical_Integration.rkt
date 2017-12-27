#lang racket
(module+ test (require rackunit))
(require test-engine/racket-tests
         "Matrix_List_Helpers.rkt")

(provide fnInt
         split-integral
         adaptive-integrate)

;; [Number -> Number] Number Number [Listof Number] [Listof Number] -> Number
;; Calculates the numerical integral of f over [a,b]
;;   1st list is list of x's to sample on the normalized interval [-1,1]
;;   2nd list is corresponding weights for each x above
;;   fnInt(f, a->b) = (1/2)(b-a)*sum_over_all_samples(weight_i*f(sample_i))
;; ASSUME: Lists have equal length, a<b
(module+ test
  (check-equal? (fnInt (λ (x) 2) 0 2 (list -1) (list 1)) 4)
  (check-equal? (fnInt (λ (x) 2) 0 2 (list 0) (list 1)) 4)
  (check-equal? (fnInt (λ (x) 2) 0 2 (list 1) (list 1)) 4)
  (check-equal? (fnInt (λ (x) 2) -2 4 (list -1) (list 1)) 12)
  (check-equal? (fnInt (λ (x) x) 1 2 (list -1) (list 1)) 1)
  (check-equal? (fnInt (λ (x) x) 0 2 (list 1) (list 1)) 4)
  (check-equal? (fnInt (λ (x) x) 0 2 (list -1 1) (list 1/2 1/2)) 2)
  (check-equal? (fnInt (λ (x) (+ x 2)) -1 1 (list -1 0 1) (list 1/3 1/3 1/3))
                4))
  
(define (fnInt f a b norm_samples weights)
  (local ((define int-length (- b a))
          (define scale (/ int-length 2))
          (define mid-pt (/ (+ a b) 2))
          (define samples (map (λ (x) (+ (* x scale) mid-pt))
                               norm_samples))
          (define values (map f samples))
          (define intermed-res (dot-prod weights values)))
    (* int-length intermed-res)))

;;-- Samples and Weights, along with highest degree polynomial
;;--   that will yield an accurate numerical solution

;;           Samples                     Weights         Highest Degree
;;           -1 0 1                    1/6 4/6 1/6            2nd
;Gaussian
;;    (sqrt 1/3) (- (sqrt 1/3))          1/2 1/2              3rd
;;   (sqrt 3/5) 0 (- (sqrt 3/5))     5/18 8/18 5/18           5th    


;; [Number -> Number] Number Number Number -> Number
;; Calculates integral of f on [a,b] by first splitting into n sub-intervals
;;   and summing fnInt over each
(module+ test
  (check-equal? (split-integral (λ (x) 2) -1 1 3)
                (fnInt (λ (x) 2) -1 1
                       (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       '(5/18 8/18 5/18)))
  (check-equal? (split-integral (λ (x) (+ (* 3 x) 4)) 0 6 4)
                (fnInt (λ (x) (+ (* 3 x) 4)) 0 6
                       (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       '(5/18 8/18 5/18)))
  (check-within (split-integral
                 (λ (x) (+ 3 (* 4 x) (* 5 (sqr x)) (* 6 (expt x 3)))) -4 4 8)
                (fnInt (λ (x) (+ 3 (* 4 x) (* 5 (sqr x)) (* 6 (expt x 3)))) -4 4
                       (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       (list 5/18 8/18 5/18))
                (expt 10 -14)))
(define (split-integral f a b n)
  (local ((define part (/ (- b a) n))
          (define new-partition (build-list (add1 n) (λ (x) (+ a (* x part)))))
          ; [Listof Number] -> Number
          ; Integrates function over each interval pair in list, sums together
          ;   Uses 3 point Gaussian weights and samples
          (define (get-sols lon)
            (cond
              [(empty? (rest lon)) 0]
              [else (+ (fnInt f (first lon) (second lon)
                              (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                              (list 5/18 8/18 5/18)) (get-sols (rest lon)))])))
    (get-sols new-partition)))


;; [Number -> Number] Number Number Number
;; Return integral of f on [a,b] approx to within tolerance
;; Algorithm:
;;   Calculates fnInt on [a,b] using 2-point Gaussian and 3-point Gaussian
;;   if the two are within tolerance, return 3-point result
;;   else recur on each half of [a,b] and sum result
;; ASSUME a<b
(module+ test
  (check-equal? (adaptive-integrate (λ (x) 2) -1 1 (expt 10 -13))
                (fnInt (λ (x) 2) -1 1 (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       '(5/18 8/18 5/18)))
  (check-equal? (adaptive-integrate (λ (x) (+ (* 3 x) 4)) 0 6 (expt 10 -13))
                (fnInt (λ (x) (+ (* 3 x) 4)) 0 6
                       (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       '(5/18 8/18 5/18)))
  (check-equal? (adaptive-integrate (λ (x) (+ 3 (* 4 x) (* 5 (sqr x)) (* 6 (expt x 3)))) -4 4 (expt 10 -13))
                (fnInt (λ (x) (+ 3 (* 4 x) (* 5 (sqr x)) (* 6 (expt x 3)))) -4 4
                       (list (sqrt 3/5) 0 (- (sqrt 3/5)))
                       (list 5/18 8/18 5/18))))
(define (adaptive-integrate f a b tol)
  (local ((define try1 (fnInt f a b (list (sqrt 1/3) (- (sqrt 1/3))) '(1/2 1/2)))
          (define try2 (fnInt f a b (list (sqrt 3/5) 0 (- (sqrt 3/5))) '(5/18 8/18 5/18)))
          (define midpt (/ (+ a b) 2))
          (define new-tol (/ tol 2)))
    (cond
      [(<= (abs (- try1 try2)) tol) try2]
      [else (+ (adaptive-integrate f a midpt new-tol) (adaptive-integrate f midpt b new-tol))])))