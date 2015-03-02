(in-package :bioinformatics)

(defun maybe-trim-eol (line)
  "depending on the platform and the git config, we might or might not have ^M at end of line, so check and remove if needed."
  (if (char= #\Return (elt line (1- (length line))))
      (subseq line 0 (1- (length line)))
      line))

(defun read-file-lines (filename)
  "list of lines from file"
  (iter (for line in-file filename using #'read-line)
	(collect line)))

(defun subseq-of-length (seq start length)
  (subseq seq start (+ start length)))

(defparameter *base-values* '((#\A . 0) (#\C . 1) (#\G . 2) (#\T . 3)) "used when converting a dna string to a number")
(defun base-to-value (base)
  "ATCG->2 bits."
  (cdr (assoc base *base-values*  :test #'char=)))
(defun value-to-base (value)
  "2 bits->ATCG"
  (car (rassoc value *base-values*)))
(defun pattern-to-number (pattern)
  "ATCG->{00 01 10 11}"
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


(defparameter *base-complements* '((#\A . #\T) (#\T . #\A) (#\C . #\G) (#\G . #\C)))
(defun base-complement (base)
  (cdr (assoc base *base-complements* :test #'char=)))
(defun reverse-complement (pattern)
  (reverse (map 'string #'base-complement pattern)))
(defun reverse-complement-file (filename)
  (let ((lines (read-file-lines filename)))
    (reverse-complement (first lines))))

(defun find-pattern-occurences (genome kmer)
  "returns a list of positions in the genome where kmer was found"
  (let ((kmer-length (length kmer)))
   (iter (for position from 0 to (- (length genome) kmer-length))
	 (when (string= genome kmer :start1 position :end1 (+ position kmer-length))
	   (collect position)))))
(defun find-pattern-occurences-from-file (filename)
  "filename names a file with 2 lines: the genome, and the kmer to be found.
returns a list of positions in the genome where kmer was found."
  (let* ((lines (read-file-lines filename))
	 (kmer (first lines))
	 (genome (second lines)))
    (find-pattern-occurences genome kmer)))
