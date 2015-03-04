(in-package :bioinformatics)

(defun compute-frequency-array-from-file (filename)
  (let* ((lines (read-file-lines filename))
	 (genome (first lines))
	 (kmer-length (parse-integer (second lines))))
    (compute-frequency-array genome kmer-length)))


(defun find-vibrio-pattern-occurences (pattern)
  (let ((vibrio (first (read-file-lines "c:/Users/serge.demarre/Downloads/Vibrio_cholerae.txt"))))
    (find-pattern-occurences vibrio pattern)))

(defun find-clump-in-genome (genome kmer-length sliding-window-length min-occurences)
  (iter (for sliding-window-position from 0 to (- (length genome) sliding-window-length))
	(multiple-value-bind (kmers occurences)
	    (most-frequent-kmer (subseq genome sliding-window-position sliding-window-length) kmer-length)
	  (when (>= occurences min-occurences)
	    (collect kmers)))))
