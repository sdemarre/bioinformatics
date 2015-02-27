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

(defun pattern-to-number (pattern)
  "ACGT->{00 01 10 11}"
  (let ((result 0))
    (iter (for letter in-string pattern)
	  (setf result (* result 4))
	  (incf result (cond
			 ((char= letter #\A) 0)
			 ((char= letter #\C) 1)
			 ((char= letter #\G) 2)
			 ((char= letter #\T) 3))))
    result))

(defun number-to-pattern (number length)
  (let ((letters (iter (for i from 1 to length)
		       (multiple-value-bind (div rem) (floor number 4)
					     (collect (case rem
							(0 #\A)
							(1 #\C)
							(2 #\G)
							(3 #\T)))
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
