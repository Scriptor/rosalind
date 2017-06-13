(ns rosalind.core
  (:require [clojure.java.io :as io]
            [clojure.string :as str]))

(def getters
  "A collection of getters for each base."
  (mapv (fn [c] #(get % c)) [\A \C \G \T]))

(def complements
  "Map of nucleotides and their complements."
  {\A \T
   \T \A
   \C \G
   \G \C})

(defn base-count
  "Returns a vector of the count for each base in the given string."
  [s]
  (->> s
       frequencies
       ((apply juxt getters))))

(defn transcribe
  "Given a DNA string return the transcribed RNA sequence."
  [s]
  (str/replace s #"T" "U"))

(defn rev-complement
  "Returns the reverse complement of the given DNA sequence."
  [s]
  (->> s
       (replace complements)
       reverse
       (apply str)))

(defn point-mut-count
  "Given two DNA sequences count the number of point mutation differences."
  [s t]
  (->> (map not= s t)
       (filter identity)
       count))

;; To answer this question we need to consider all possible pairings:
;; (note that each of the below needs to be doubled)
;; - Two hodos: (k / total) + (k-1 / dtotal)
;; - hodo and het: (k / total) + (m / dtotal)
;; - two hets: (m / total) + (m - 1) / dtotal
;; - het and hore: (m / total) + (n / dtotal)
;; - two hores: (n / total) + (n - 1 / dtotal)
;; - hodo and hore: (k / total) + (n / dtotal)
(defn dominant-chance
  "Given population numbers for:
     k - homozygous dominant
     m - heterozygous
     n - homozygous recessive
   Return the probability that any random pair's offspring will have at least
   one copy of the dominant allele."
  [k m n]
  (let [total (+ k m n)
        kk (*   (/ k total) (/ (dec k) (dec total)))
        km (* 2 (/ k total) (/ m (dec total)))
        mm (*   (/ m total) (/ (dec m) (dec total)))
        mn (* 2 (/ m total) (/ n (dec total)))
        nn (*   (/ n total) (/ (dec n) (dec total)))
        kn (* 2 (/ k total) (/ n (dec total)))]
    (+ kk km (* 0.75 mm) (* 0.5 mn) kn)))

(defn run
  "Takes a function and executes it against the dataset resource."
  [f]
  (let [s (slurp (io/resource "dna.txt"))]
    (f s)))
