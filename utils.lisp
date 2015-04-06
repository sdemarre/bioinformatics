(in-package :bioinformatics)

(defun maybe-trim-eol (line)
  "depending on the platform and the git config, we might or might not have ^M at end of line, so check and remove if needed."
  (if (and (> (length line) 0) (char= #\Return (elt line (1- (length line)))))
      (subseq line 0 (1- (length line)))
      line))

(defun read-file-lines (filename)
  "list of lines from file"
  (iter (for line in-file filename using #'read-line)
	(collect (maybe-trim-eol line))))

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

(defparameter *codon-amino-data*
  '(("Isoleucine" "I" ("ATT" "ATC" "ATA"))
    ("Leucine" "L" ("CTT" "CTC" "CTA" "CTG" "TTA" "TTG"))
    ("Valine" "V" ("GTT" "GTC" "GTA" "GTG"))
    ("Phenylalanine" "F" ("TTT" "TTC"))
    ("Methionine" "M" ("ATG"))
    ("Cysteine" "C" ("TGT" "TGC"))
    ("Alanine" "A" ("GCT" "GCC" "GCA" "GCG"))
    ("Glycine" "G" ("GGT" "GGC" "GGA" "GGG"))
    ("Proline" "P"  ("CCT" "CCC" "CCA" "CCG"))
    ("Threonine" "T" ("ACT" "ACC" "ACA" "ACG"))
    ("Serine" "S" ("TCT" "TCC" "TCA" "TCG" "AGT" "AGC"))
    ("Tyrosine" "Y" ("TAT" "TAC"))
    ("Tryptophan" "W" ("TGG"))
    ("Glutamine" "Q" ("CAA" "CAG"))
    ("Asparagine" "N" ("AAT" "AAC"))
    ("Histidine" "H" ("CAT" "CAC"))
    ("Glutamic acid" "E" ("GAA" "GAG"))
    ("Aspartic acid" "D" ("GAT" "GAC"))
    ("Lysine" "K" ("AAA" "AAG"))
    ("Arginine" "R" ("CGT" "CGC" "CGA" "CGG" "AGA" "AGG"))
    ("Stop codon" "Stop" ("TAA" "TAG" "TGA"))))


(defparameter *dna-codon-amino-hash*
  (let ((result (make-hash-table :test #'equal)))
    (iter (for (name code codons) in *codon-amino-data*)
	  (iter (for codon in codons)
		(setf (gethash codon result) (list name code))))
    result))

(defparameter *protein-to-codons-hash*
  (let ((result (make-hash-table :test #'equal)))
    (iter (for (name code codons) in *codon-amino-data*)
	  (when (not (string= code "Stop"))
	    (setf (gethash (elt code 0) result) codons)))
    result))
(defun rna-to-dna (rna)
  (cl-ppcre:regex-replace-all "U" rna "T"))
(defun dna-to-rna (dna)
  (cl-ppcre:regex-replace-all "T" dna "U"))

(defparameter *rna-codon-amino-hash* (let ((result (make-hash-table :test #'equal)))
				       (iter (for (key value) in-hashtable *dna-codon-amino-hash*)
					     (setf (gethash (dna-to-rna key) result) value))
				       result))

(defun dna-codon-to-amino-acid (codon)
  (second (gethash codon *dna-codon-amino-hash*)))

(defun is-dna-stop-codon-p (codon)
  (string= "Stop" (dna-codon-to-amino-acid codon)))

(defun rna-codon-to-amino-acid (codon)
  (second (gethash codon *rna-codon-amino-hash*)))

(defun do-on-codons (input fun &optional (start 0))
  (iter (for position from start by 3)
	(while (< (+ position 2) (length input)))
	(collect (funcall fun (subseq-of-length input position 3)))))
(defun dna-to-protein (dna)
  (do-on-codons dna #'dna-codon-to-amino-acid))
(defun rna-to-protein (rna)
  (do-on-codons rna #'rna-codon-to-amino-acid))


(defun combine-strings (list-of-strings)
  (apply #'concatenate 'string list-of-strings))
(defun fasta-name-line-p (line)
  (and (> (length line) 0)
       (char= #\> (elt line 0))))
(defun fasta-combine-lines (lines)
  (let ((name)
	(dna-strings (make-hash-table :test #'equal)))
    (iter (for line in lines)
	  (if (fasta-name-line-p line)
	      (setf name (subseq line 1)) ;; get rid of initial #\>
	      (push line (gethash name dna-strings))))
    (iter (for (key value) in-hashtable dna-strings)
	  (collect (list key (combine-strings (reverse value)))))))
(defun read-fasta-lines (filename)
  "returns a list of (fasta-id string), as found in filename"
  (let ((lines (read-file-lines filename)))
    (fasta-combine-lines lines)))

(defun profile-matrix (dna-strings)
  (let ((a-counts (make-array (length (first dna-strings)) :initial-element 0))
	(t-counts (make-array (length (first dna-strings)) :initial-element 0))
	(c-counts (make-array (length (first dna-strings)) :initial-element 0))
	(g-counts (make-array (length (first dna-strings)) :initial-element 0)))
    (iter (for dna in dna-strings)
	  (iter (for index index-of-string dna)
		(for letter in-string dna)
		(cond ((char= #\A letter) (incf (elt a-counts index)))
		      ((char= #\T letter) (incf (elt t-counts index)))
		      ((char= #\C letter) (incf (elt c-counts index)))
		      ((char= #\G letter) (incf (elt g-counts index))))))
    (list a-counts t-counts c-counts g-counts)))

(defun consensus (profile-matrix)
  (destructuring-bind (a-counts t-counts c-counts g-counts) profile-matrix
    (iter (for index index-of-vector a-counts)
	  (collect
	      (iter (for letter-data in (list (cons #\A (elt a-counts index))
					      (cons #\T (elt t-counts index))
					      (cons #\C (elt c-counts index))
					      (cons #\G (elt g-counts index))))
		    (finding (car letter-data) maximizing (cdr letter-data)))))))


(defun fib (n)
  (let ((fn-1 1)
	(fn-2 0))
    (iter (for i from 1 to n)
	  (psetf fn-2 fn-1
		 fn-1 (+ fn-1 fn-2)))
    fn-2))

(defmacro destructuring-bind-integers (vars expr &body body)
  `(destructuring-bind ,vars (mapcar #'parse-integer (split-sequence:split-sequence #\Space ,expr))
     ,@body))

(defun parse-list (string fun)
  (mapcar fun (split-sequence:split-sequence #\Space string)))

(defun parse-integer-list (string)
  (parse-list string #'parse-integer))
(defun parse-integer-vector (string)
  (coerce (parse-integer-list string) 'vector))
(defun print-integer-list (list &optional (stream t))
  (format stream "~{~a~^ ~}~%" list))
(defun print-integer-vector (vector &optional (stream t))
  (print-integer-list (iter (for e in-vector vector) (collect e)) stream))

(defun parse-float-list (string)
  (let ((*read-default-float-format* 'double-float))
    (parse-list string #'parse-number:parse-real-number)))
(defun print-float-list (list &optional (stream t))
  (format stream "~{~f~^ ~}~%" list))

(defun fact (n)
  (iter (for i from 1 to n)
	(multiply i)))

(defun combinations (n p)
  (/ (fact n) (fact p) (fact (- n p))))

(defun permutations (n p)
  (/ (fact n) (fact (- n p))))

(defun probability-for-events (event-type-profile event-type event-count number-experiments &optional (max-event-count event-count))
  (if (> event-count number-experiments)
      0
      (let ((event-probability (gethash event-type event-type-profile 0)))
	(iter (for current-event-count from event-count to max-event-count)
	      (summing (* (combinations number-experiments current-event-count)
			  (expt event-probability current-event-count)
			  (expt (- 1 event-probability) (- number-experiments current-event-count))))))))

(defparameter *monoisotopic-protein-mass*
  (let ((hash (make-hash-table :test #'equal)))
    (iter (for (protein mass) in '((#\A   71.03711d0)
				   (#\C   103.00919d0)
				   (#\D   115.02694d0)
				   (#\E   129.04259d0)
				   (#\F   147.06841d0)
				   (#\G   57.02146d0)
				   (#\H   137.05891d0)
				   (#\I   113.08406d0)
				   (#\K   128.09496d0)
				   (#\L   113.08406d0)
				   (#\M   131.04049d0)
				   (#\N   114.04293d0)
				   (#\P   97.05276d0)
				   (#\Q   128.05858d0)
				   (#\R   156.10111d0)
				   (#\S   87.03203d0)
				   (#\T   101.04768d0)
				   (#\V   99.06841d0)
				   (#\W   186.07931d0)
				   (#\Y   163.06333d0)))
	  (setf (gethash protein hash) mass))
    hash))

(defun amino-acid-mass (aa)
  (gethash aa *monoisotopic-protein-mass* 0))

(defun protein-mass (protein-string)
  "input: string of aa, output: double float"
  (iter (for aa in-vector protein-string)
	(summing (amino-acid-mass aa))))

(defun count-possible-rna-sources-for-protein (protein-string)
  (iter (for protein in-vector protein-string)
	  (multiplying (length (gethash protein *protein-to-codons-hash*)))))

(defun get-uniprot-fasta (uniprot-id)
  "gets the fasta file that describes this protein from www.uniprot.org"
  (multiple-value-bind (fasta-string http-status rest)
      (drakma:http-request (format nil "http://www.uniprot.org/uniprot/~a.fasta" uniprot-id))
    (declare (ignorable rest))
    (if (= 200 http-status)
	(fasta-combine-lines (split-sequence:split-sequence #\Newline fasta-string))
	(error "http get failed, ~a" http-status))))
(defun get-uniprot-fasta-cached (uniprot-id)
  "gets the fasta file that describes this protein from local cache or from www.uniprot.org"  
  (let ((prot-filename (format nil "uniprot/~a" uniprot-id)))
    (if (cl-fad:file-exists-p prot-filename)
	(read-fasta-lines prot-filename)
	(progn
	  (with-open-file (fasta-output prot-filename :direction :output :if-exists :supersede)
	    (let ((fasta-data  (get-uniprot-fasta uniprot-id)))
	      (format fasta-output "~{~A~%~}" fasta-data)
	      fasta-data))))))

(defun protein-motif-to-regexp (protein-motif)
  (cl-ppcre:regex-replace-all "}" (cl-ppcre:regex-replace-all "{" protein-motif "[^") "]"))

(defun protein-motif-positions-wrong (protein-motif protein)
  ;; this would be ok if we didn't have to find overlapping motifs...
  (let* ((motif-rx (protein-motif-to-regexp protein-motif))
	 (matches (cl-ppcre:all-matches motif-rx protein)))
    (iter (for position on matches by #'cddr)
	  (collect (1+ (car position))))))

(defun protein-n-glycosilation-positions (protein)
  (iter (for position from 0 to (- (length protein) 4))
	(symbol-macrolet ((aa0 (elt protein position))
			  (aa1 (elt protein (1+ position)))
			  (aa2 (elt protein (+ position 2)))
			  (aa3 (elt protein (+ position 3))))
	 (when (and (char= aa0 #\N)
		    (not (char= aa1 #\P))
		    (or (char= aa2 #\S) (char= aa2 #\T))
		    (not (char= aa3 #\P)))
	   (collect (1+ position))))))


(defun position-of-substring (string1 start end string2)
  "find position of (subseq string1 start end) in string2"
  (let ((l2 (length string2))
	(l1 (length string1)))
    (flet ((substring-at-position-p (idx)
	     (iter (for idx1 from start to end)
		   (for idx2 from idx)
		   (always (char= (elt string1 idx1) (elt string2 idx2))))))
      (iter (for idx2 from 0 to (- l2 l1))
	    (when (substring-at-position-p idx2)
	      (return idx2))))))

(defun make-random-dna-string (length &optional (gc 0.5))
  (let ((result (make-array length :element-type 'character)))
    (iter (for idx index-of-string result)
	  (setf (elt result idx)
		(if (> (random 1.0d0) gc)
		    (if (zerop (random 2)) #\A #\T)
		    (if (zerop (random 2)) #\C #\G))))
    result))

(defun make-random-substrings (dna-string count length-or-fun &optional print)
  (when print
    (format t "~a~%" dna-string))
  (flet ((dots (n) (iter (repeat n) (collect #\.))))
   (iter (repeat count)
	 (let* ((length (if (numberp length-or-fun) length-or-fun (funcall length-or-fun)))
		(pos  (random (- (length dna-string) length)))
		(substr (subseq-of-length dna-string pos length)))
	   (when print
	     (format t "~{~a~}~a~{~a~}~%" (dots pos) substr (dots (- (length dna-string) pos (length substr)))))
	   (collect substr)))))
