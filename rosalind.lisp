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

(defun rosalind-gc-content (dna)
  (* 1.0 (/ (iter (for letter in-vector dna)
		  (counting (or (char= letter #\C) (char= letter #\G))))
	    (length dna))))
(defun rosalind-gc ()
  (let ((lines (read-fasta-lines "rosalind_gc.txt")))
    (iter (for (name dna) in lines)
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

(defun rosalind-rna-to-protein ()
  (let ((lines (read-file-lines "rosalind_prot.txt")))
    (coerce (mapcar #'(lambda (s) (if (= 1 (length s))
				      (elt s 0)
				      ""))
		    (butlast (rna-to-protein (first lines))))
	    'string)))

(defun rosalind-find-motifs ()
  (let* ((lines (read-file-lines "rosalind_subs.txt"))
	 (haystack (first lines))
	 (needle (second lines)))
    (mapcar #'1+ (find-pattern-occurences haystack needle))))

(defun rosalind-consensus-and-profile ()
  (let ((fasta-lines (read-fasta-lines "rosalind_cons.txt")))
    (list (consesus (mapcar #'second fasta-lines))
	  (profile-matrix (mapcar #'second fasta-lines)))))
