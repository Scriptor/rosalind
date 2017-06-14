(ns rosalind.fasta
  (:require [clojure.java.io :as io]))

(defn parsed-fasta
  "Given a sequence of FASTA-formatted lines return a map of
   labels to their corresponding strings."
  [lines]
  (->> lines
       (partition-by #(= (nth % 0) \>))
       (map (partial apply str))
       (apply hash-map)))

(defn import-fasta
  [f]
  (parsed-fasta (line-seq (io/reader f))))
