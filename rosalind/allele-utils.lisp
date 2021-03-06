(in-package :bioinformatics)

(defun sort-allele (allele)
  "put capitals first, e.g. aA -> Aa, AA -> AA, Aa ->Aa"
  (let (a1g1 a1g2)
    (if (char> (elt allele 0) (elt allele 1))
	(setf a1g2 (elt allele 0) a1g1 (elt allele 1))
	(setf a1g1 (elt allele 0) a1g2 (elt allele 1)))
    (format nil "~a~a" a1g1 a1g2)))
(defun number-alleles (alleles)
  (/ (length alleles) 2))
(defun select-allele (alleles number)
  (subseq alleles (* 2 number) (* 2 (1+ number))))
(defun sort-alleles (alleles)
  "put capitals first for every allele, e.g. aAbB -> AaBb"
  (format nil "~{~a~}" (iter (for i from 0 to (1- (number-alleles alleles)))
			     (while (< i (length alleles)))
			     (collect (sort-allele (select-allele alleles i))))))
(defun child-genotype (parent1-genotype parent2-genotype gene-selectors)
  "gene-selectors is a list of pairs of 0's and 1's, selecting gene 0 or 1 from the allele. you need as many pairs as there are alleles in the parents.
e.g.\"AaBb\" \"aaBB\" ((0 1) (1 0)) -> AabB"
  (format nil "~{~a~}" (iter (for allele-idx from 0 to (1- (number-alleles parent1-genotype)))
			     (for (p1-selector p2-selector) in gene-selectors)
			     (collect (format nil "~a~a"
					      (elt (select-allele parent1-genotype allele-idx) p1-selector)
					      (elt (select-allele parent2-genotype allele-idx) p2-selector))))))
(defun punnett-square (parent1-genotype parent2-genotype)
  "parent alleles: string of 4 chars, e.g.\"AaBb\""
  (let ((result (make-array (list 4 4))))
    (iter (for row from 0 to 3)
	  (iter (for column from 0 to 3)
		(let ((child-genotype (child-genotype parent1-genotype parent2-genotype
						  `((,(floor row 2) ,(floor column 2))
						    (,(mod row 2) ,(mod column 2))))))
		  (setf (aref result row column) child-genotype))))
    result))

(defun possible-unique-child-genotypes (parent1-genotype parent2-genotype)
  (let ((possible-children (possible-child-genotypes parent1-genotype parent2-genotype))
	(unique-children (make-hash-table :test #'equal)))
    (iter (for child in possible-children)
	  (setf (gethash (sort-alleles child) unique-children) 1))
    (alexandria:hash-table-keys unique-children)))
(defun get-bit (index number)
  (if (logbitp index number) 1 0))
(defun possible-child-genotypes (parent1-genotype parent2-genotype)
  (let ((num-alleles (number-alleles parent1-genotype))
	child-genotypes)
    (iter (for p1-gene-selector from 0 to (1- (expt 2 num-alleles)))
	  (iter (for p2-gene-selector from 0 to (1- (expt 2 num-alleles)))
		(let ((gene-selectors (iter (for bit from 0 to (1- num-alleles))
					    (collect (list (get-bit bit p1-gene-selector)
							   (get-bit bit p2-gene-selector))))))
		  (push (child-genotype parent1-genotype parent2-genotype gene-selectors) child-genotypes))))
    child-genotypes))

(defun child-genotype-profile (parent1-genotype parent2-genotype)
  (let ((genotype-counts (make-hash-table :test #'equal))
	(total 0))
    (iter (for genotype in (possible-child-genotypes parent1-genotype parent2-genotype))
	  (incf (gethash (sort-alleles genotype) genotype-counts 0))
	  (incf total))
    (alexandria:maphash-keys #'(lambda (k) (setf (gethash k genotype-counts) (/ (gethash k genotype-counts) total))) genotype-counts)
    genotype-counts))

(defun compute-next-generation-profile (current-generation-profile mate-genotype)
  "current-generation-profile is a hash-table with probability for every genotype"
  (let ((next-generation-profile (make-hash-table :test #'equal)))
    (iter (for (genotype probability) in-hashtable current-generation-profile)
	  (iter (for (child-genotype new-prob) in-hashtable (child-genotype-profile genotype mate-genotype))
		(incf (gethash child-genotype next-generation-profile 0) (* probability new-prob))))
    next-generation-profile))
