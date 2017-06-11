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

(defn run
  "Takes a function and executes it against the dataset resource."
  [f]
  (let [s (slurp (io/resource "dna.txt"))]
    (f s)))
