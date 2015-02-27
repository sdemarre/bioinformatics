(in-package :bioinformatics)

(defun read-file-lines (filename)
  (iter (for line in-file filename using #'read-line)
	(collect line)))

(defun subseq-of-length (seq start length)
  (subseq seq start (+ start length)))
(defun most-frequent-kmer (genome k)
  (let ((kmer-freq (make-hash-table :test #'equal)))
    (let ((max-freq (iter (for position from 0 to (- (length genome) k))
			  (let ((kmer (subseq-of-length genome position k)))
			    (maximizing (incf (gethash kmer kmer-freq 0)))))))
      (list
       (iter (for (kmer freq) in-hashtable kmer-freq)
	     (when (= freq max-freq)
	       (collect kmer)))
       max-freq))))


(defun count-kmer (genome kmer)
  (iter (for position from 0 to (- (length genome) (length kmer)))
	(let ((sub (subseq-of-length genome position (length kmer))))
	  (counting (string= kmer sub)))))

(defun problem-1 (filename)
  (let* ((lines (read-file-lines filename))
	 (genome (first lines))
	 (kmer (second lines)))
    (count-kmer genome kmer)))

(defun problem-2 (filename)
  (let* ((lines (read-file-lines filename))
	 (genome (first lines))
	 (k (parse-integer (second lines))))
    (most-frequent-kmer genome k)))

(defun probability-of-pattern (string-length number-symbols pattern n)
  "Let Pr(N, A, Pattern, t) denote the probability that a string Pattern appears t or more times in a random string of length N formed from an alphabet of A letters"
 )

