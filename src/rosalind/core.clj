(ns rosalind.core
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [rosalind.fasta :as fasta]
            [clojure.pprint :as pp]))

(def getters
  "A collection of getters for each base."
  (mapv (fn [c] #(get % c)) [\A \C \G \T]))

(def complements
  "Map of nucleotides and their complements."
  {\A \T
   \T \A
   \C \G
   \G \C})

(def codon-map
  "Map of RNA codons to their amino acid codes (or to the stop code)."
  {"AUA" "I", "CAU" "H", "GAG" "E", "GUG" "V", "UGA" "Stop", "CUU" "L", "CGU" "R",
   "CAA" "Q", "AGC" "S", "GAA" "E", "GUC" "V", "GUA" "V", "UAC" "Y", "GAU" "D",
   "ACA" "T", "CUC" "L", "UCC" "S", "GCC" "A", "UAA" "Stop", "CUA" "L", "GGC" "G",
   "CGC" "R", "ACU" "T", "ACC" "T", "AGU" "S", "UUC" "F", "GGG" "G", "CGG" "R",
   "CAG" "Q", "AAC" "N", "ACG" "T", "UGC" "C", "CUG" "L", "AAG" "K", "GGU" "G",
   "AAU" "N", "UGG" "W", "CGA" "R", "CCA" "P", "CCU" "P", "AUU" "I", "CCC" "P",
   "CAC" "H", "AGA" "R", "UCU" "S", "UCA" "S", "UUU" "F", "GCG" "A", "UAG" "Stop",
   "GAC" "D", "AUG" "M", "UAU" "Y", "UGU" "C", "CCG" "P", "GGA" "G", "AGG" "R",
   "UCG" "S", "GUU" "V", "AAA" "K", "AUC" "I", "GCU" "A", "UUG" "L", "UUA" "L",
   "GCA" "A"})

(defn stop?
  "Returns whether the given string corresponds to a stop codon as used in
   `codon-map'"
  [s]
  (= s "Stop"))

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
       str/join))

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

(defn wascally
  "Fibonacci except each rabbit pair produces k more pairs instead of just one.
   n is # of months/generations."
  [n k]
  (loop [n       (dec n)
         adults  0
         bunnies 1]
    (if (zero? n)
      (+ adults bunnies)
      (recur (dec n) (+ adults bunnies) (* adults k)))))

(defn subs-locs
  "Return the positions of the given DNA substring t in the larger string s.
   Note that the positions must be 1-based."
  [s t]
  (let [tlen (count t)]
    (->> s
         (iterate #(subs % 1))
         (take-while #(>= (count %) tlen))
         (keep-indexed (fn [i s']
                         (when (= t (subs s' 0 tlen))
                           ;; Increment because result must be 1-based
                           (inc i)))))))

(defn gc-content
  "Return the gc content of a given DNA string"
  [s]
  (-> (/ (count (filter #{\C \G} s))
         (count s))
      float
      (* 100)))

(defn compute-fasta-gc
  "Given a fasta file find the string with the highest GC content and
   return its label and GC percentage."
  [f]
  (let [fasta (fasta/import-fasta f)]
    (->> fasta
         (map #(vector (first %) (float (gc-content (second %)))))
         (sort-by second >)
         first)))

(defn translation
  "Given the RNA string `s' return a sequence of the amino acid codes it
   will translate to.`"
  [s]
  (->> s
       (partition 3)
       (map (comp codon-map str/join))
       (take-while (complement stop?))
       str/join))

(defn run
  "Takes a function and executes it against the dataset resource."
  [f]
  (let [s (slurp (io/resource "dna.txt"))]
    (f s)))
