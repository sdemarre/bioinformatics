(in-package :bioinformatics)


(defun rosalind-dna ()
  (let ((dna (first (read-file-lines "rosalind_dna.txt")))
	(count-hash (make-hash-table)))
    (iter (for letter in-vector dna)
	  (incf (gethash letter count-hash 0)))
    (format nil "~{~a ~}~%"
	    (iter (for letter in-vector "ACGT")
		  (collect (gethash letter count-hash))))))

(defun rosalind-rna ()
  (let ((dna (first (read-file-lines "rosalind_rna.txt"))))
    (cl-ppcre:regex-replace-all "T" dna "U")))

(defun rosalind-revc ()
  (let ((dna (first (read-file-lines "rosalind_revc.txt"))))
    (reverse-complement dna)))


(defun ros-fib (generations pairs-per-litter fnm2 fnm1)
  (cond ((= generations 0) fnm1)
	((= generations 1) fnm2)
	(t (ros-fib (1- generations) pairs-per-litter fnm1 (+ fnm1 (* pairs-per-litter fnm2))))))
(defun rosalind-fib ()
  (let ((data (first (read-file-lines "rosalind_fib.txt"))))
    (destructuring-bind (generations pairs-per-litter) (split-sequence:split-sequence #\Space data)
      (ros-fib (parse-integer generations) (parse-integer pairs-per-litter) 1 1))))

(defun combine-strings (list-of-strings)
  (apply #'concatenate 'string list-of-strings))
(defun fasta-name-line-p (line)
  (and (> (length line) 0)
       (char= #\> (elt line 0))))
(defun rosalind-fasta-combine (lines)
  (let ((name)
	(dna-strings (make-hash-table :test #'equal)))
    (iter (for line in lines)
	  (if (fasta-name-line-p line)
	      (setf name line)	 
	      (push line (gethash name dna-strings))))
    (iter (for (key value) in-hashtable dna-strings)
	  (collect (list key (combine-strings value))))))
(defun rosalind-gc-content (dna)
  (* 1.0 (/ (iter (for letter in-vector dna)
		  (counting (or (char= letter #\C) (char= letter #\G))))
	    (length dna))))
(defun rosalind-gc ()
  (let ((lines (read-file-lines "rosalind_gc.txt")))
    (iter (for (name dna) in (rosalind-fasta-combine lines))
	  (finding (list name (* 100 (rosalind-gc-content dna))) maximizing (rosalind-gc-content dna)))))


(defun hamming-distance (dna1 dna2)
  (iter (for letter1 in-vector dna1)
	(for letter2 in-vector dna2)
	(counting (not (char= letter1 letter2)))))
(defun rosalind-hamming ()
  (let* ((lines (read-file-lines "rosalind_hamm.txt"))
	 (dna1 (first lines))
	 (dna2 (second lines)))
    (hamming-distance dna1 dna2)))

(defun rosalind-mendelian (homozygous-dominant-count heterozygous-count homozygous-recessive-count)
  (let ((k homozygous-dominant-count)
	(l heterozygous-count)
	(m homozygous-recessive-count))
    (/ (+ (* k (1- k)) (* 2 k l) (* 2 k m) (* m l)  (/ (* 3 l (1- l)) 4)) (* (+ k l m) (+ k l m -1)))))

(defparameter *rosalind-rna-codon-hash*
  (alexandria:plist-hash-table '("UUU" "F" "CUU" "L" "AUU" "I" "GUU" "V" "UUC" "F"
				 "CUC" "L" "AUC" "I" "GUC" "V" "UUA" "L" "CUA" "L"
				 "AUA" "I" "GUA" "V" "UUG" "L" "CUG" "L" "AUG" "M"
				 "GUG" "V" "UCU" "S" "CCU" "P" "ACU" "T" "GCU" "A"
				 "UCC" "S" "CCC" "P" "ACC" "T" "GCC" "A"
				 "UCA" "S" "CCA" "P" "ACA" "T" "GCA" "A"
				 "UCG" "S" "CCG" "P" "ACG" "T" "GCG" "A"
				 "UAU" "Y" "CAU" "H" "AAU" "N" "GAU" "D"
				 "UAC" "Y" "CAC" "H" "AAC" "N" "GAC" "D"
				 "UAA" "Stop" "CAA" "Q" "AAA" "K" "GAA" "E"
				 "UAG" "Stop" "CAG" "Q" "AAG" "K" "GAG" "E"
				 "UGU" "C" "CGU" "R" "AGU" "S" "GGU" "G"
				 "UGC" "C" "CGC" "R" "AGC" "S" "GGC" "G"
				 "UGA" "Stop" "CGA" "R" "AGA" "R" "GGA" "G"
				 "UGG" "W" "CGG" "R" "AGG" "R" "GGG" "G")))
