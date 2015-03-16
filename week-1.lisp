(in-package :bioinformatics)

(defun most-frequent-kmers (genome kmer-length)
  (let* ((kmer-freq-hash (make-hash-table :test #'equal))
	 (max-kmer-freq 
	  (iter (for position from 0 to (- (length genome) kmer-length))
		(let ((kmer (subseq-of-length genome position kmer-length)))
		  (maximizing (incf (gethash kmer kmer-freq-hash 0)))))))
    (values
     (iter (for (kmer freq) in-hashtable kmer-freq-hash)
	   (when (= freq max-kmer-freq)
	     (collect kmer))))))

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

