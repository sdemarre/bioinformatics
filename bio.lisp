(in-package :bio)


(defun maybe-trim-eol (line)
  (if (char= #\Return (elt line (1- (length line))))
      (subseq line 0 (1- (length line)))
      line))

(defun read-file-lines (filename)
  (iter (for line in-file filename using #'read-line)
	(collect (maybe-trim-eol line))))

(defun subseq-of-length (seq start length)
  (subseq seq start (+ start length)))
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

(defparameter *base-values* '((#\A . 0) (#\C . 1) (#\G . 2) (#\T . 3)))
(defparameter *base-complements* '((#\A . #\T) (#\T . #\A) (#\C . #\G) (#\G . #\C)))
(defun base-to-value (base)
  (cdr (assoc base *base-values*  :test #'char=)))
(defun value-to-base (value)
  (car (rassoc value *base-values*)))
(defun pattern-to-number (pattern)
  "ACGT->{00 01 10 11}"
  (let ((result 0))
    (iter (for letter in-string pattern)
	  (setf result (* result 4))
	  (incf result (base-to-value letter)))
    result))

(defun number-to-pattern (number length)
  (let ((letters (iter (for i from 1 to length)
		       (multiple-value-bind (div rem) (floor number 4)
					     (collect (value-to-base rem))
					     (setf number div)))))
    (coerce (reverse letters) 'string)))

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

(defun base-complement (base)
  (cdr (assoc base *base-complements* :test #'char=)))
(defun reverse-complement (pattern)
  (reverse (map 'string #'base-complement pattern)))
(defun reverse-complement-file (filename)
  (let ((lines (read-file-lines filename)))
    (reverse-complement (first lines))))

(defun find-pattern-occurences (genome kmer)
  (let ((kmer-length (length kmer)))
   (iter (for position from 0 to (- (length genome) kmer-length))
	 (when (string= genome kmer :start1 position :end1 (+ position kmer-length))
	   (collect position)))))
(defun find-pattern-occurences-from-file (filename)
  (let* ((lines (read-file-lines filename))
	 (kmer (first lines))
	 (genome (second lines)))
    (find-pattern-occurences genome kmer)))
(defun find-vibrio-pattern-occurences (pattern)
  (let ((vibrio (first (read-file-lines "c:/Users/serge.demarre/Downloads/Vibrio_cholerae.txt"))))
    (find-pattern-occurences vibrio pattern)))

(defun find-clump-in-genome (genome kmer-length sliding-window-length min-occurences)
  (iter (for sliding-window-position from 0 to (- (length genome) sliding-window-length))
	(multiple-value-bind (kmers occurences)
	    (most-frequent-kmer (subseq genome sliding-window-position sliding-window-length) kmer-length)
	  (when (>= occurences min-occurences)
	    (collect kmers)))))