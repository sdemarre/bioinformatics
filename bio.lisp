(in-package :bio)

(defun most-frequent-kmer (genome kmer-length)
  (let* ((kmer-freq (make-hash-table :test #'equal))
	 (max-kmer-freq
	  (iter (for position from 0 to (- (length genome) kmer-length))
		(let ((kmer (subseq-of-length genome position kmer-length)))
		  (maximize (incf (gethash kmer kmer-freq 0)))))))
    (values
     (iter (for (kmer freq) in-hashtable kmer-freq)
	   (when (= freq max-kmer-freq)
	     (collect kmer)))
     max-kmer-freq)))

(defun problem-2 (dataset-filename)
  (let* ((lines (read-file-lines dataset-filename))
	 (genome (first lines))
	 (kmer-length (parse-integer (second lines))))
    (most-frequent-kmer genome kmer-length)))

(defun compute-frequency-array (genome kmer-length)
  (let ((frequency-array (make-array (expt 4 kmer-length) :initial-element 0))
	(genome-length (length genome)))
    (iter (for position from 0 to (- genome-length kmer-length))
	  (incf (elt frequency-array (pattern-to-number (subseq-of-length genome position kmer-length)))))
    frequency-array))

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
