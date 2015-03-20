(in-package :bioinformatics)

(define-rosalind-problem :dna rosalind-dna
  "counting dna nucleotides"
  (with-single-input-line (dna)
   (let ((count-hash (make-hash-table)))
     (iter (for letter in-vector dna)
	   (incf (gethash letter count-hash 0)))
     (format nil "~{~a ~}~%"
	     (iter (for letter in-vector "ACGT")
		   (collect (gethash letter count-hash)))))))

(define-rosalind-problem :rna rosalind-rna
  "transcribing dna into rna"
  (with-single-input-line (dna)
    (cl-ppcre:regex-replace-all "T" dna "U")))

(define-rosalind-problem :revc rosalind-revc
  "complementing a strand of dna"
  (with-single-input-line (dna)
    (reverse-complement dna)))

(defun ros-fib (generations pairs-per-litter &optional (fnm2 1) (fnm1 1))
  (cond ((= generations 0) fnm1)
	((= generations 1) fnm2)
	(t (ros-fib (1- generations) pairs-per-litter fnm1 (+ fnm1 (* pairs-per-litter fnm2))))))
(define-rosalind-problem :fib rosalind-fib
  "rabbits and recurrence relations"
  (with-single-input-line (problem-data)
    (destructuring-bind-integers  (generations pairs-per-litter) problem-data
      (ros-fib generations pairs-per-litter))))

(defun rosalind-gc-content (dna)
  (* 1.0 (/ (iter (for letter in-vector dna)
		  (counting (or (char= letter #\C) (char= letter #\G))))
	    (length dna))))
(define-rosalind-problem :gc rosalind-gc
  "computing gc content"
  (with-fasta-input-lines (lines)
    (iter (for (name dna) in lines)
	  (finding (list name (* 100 (rosalind-gc-content dna))) maximizing (rosalind-gc-content dna)))))

(defun hamming-distance (dna1 dna2)
  (iter (for letter1 in-vector dna1)
	(for letter2 in-vector dna2)
	(counting (not (char= letter1 letter2)))))
(define-rosalind-problem :hamm rosalind-hamming
  "counting point mutations"
  (with-input-lines (lines)
    (let* ((dna1 (first lines))
	   (dna2 (second lines)))
      (hamming-distance dna1 dna2))))

(defun rosalind-mendelian (homozygous-dominant-count heterozygous-count homozygous-recessive-count)
  (let ((k homozygous-dominant-count)
	(l heterozygous-count)
	(m homozygous-recessive-count))
    (/ (+ (* k (1- k)) (* 2 k l) (* 2 k m) (* m l)  (/ (* 3 l (1- l)) 4)) (* (+ k l m) (+ k l m -1)))))
(define-rosalind-problem :iprb mendels-first-law
  "mendel's first law"
  (with-input-lines (lines)
    (destructuring-bind-integers (homozygous-dominant-count heterozygous-count homozygous-recessive-count) (first lines)
      (float (rosalind-mendelian homozygous-dominant-count heterozygous-count homozygous-recessive-count)))))

(define-rosalind-problem :prot rosalind-rna-to-protein
  "translating rna into protein"
  (with-input-lines (lines)
    (coerce (mapcar #'(lambda (s) (if (= 1 (length s))
				      (elt s 0)
				      ""))
		    (butlast (rna-to-protein (first lines))))
	    'string)))

(define-rosalind-problem :subs rosalind-find-motifs
  "finding a motif in dna"
  (with-input-lines (lines)
    (let* ((haystack (first lines))
	   (needle (second lines)))
      (format t "~{~a~^ ~}" (mapcar #'1+ (find-pattern-occurences haystack needle))))))

(defun rosalind-consensus-and-profile (input-filename)
  (with-fasta-input-lines (fasta-lines)
    (let ((profile (profile-matrix (mapcar #'second fasta-lines))))
      (list (consensus profile) profile))))
(defun rosalind-consensus-and-profile-print (stream input-filename)
  (destructuring-bind (consensus (a-prof t-prof c-prof g-prof)) (rosalind-consensus-and-profile input-filename)
    (format stream "~a~%" (coerce consensus 'string))
    (format stream "A: ~{~A~^ ~}~%" (coerce a-prof 'list))
    (format stream "C: ~{~A~^ ~}~%" (coerce c-prof 'list))
    (format stream "G: ~{~A~^ ~}~%" (coerce g-prof 'list))
    (format stream "T: ~{~A~^ ~}~%" (coerce t-prof 'list))))
(define-rosalind-problem :cons rosalind-consensus-and-profile-to-file
  "consensus and profile"
  (with-output-to-file (output)
    (rosalind-consensus-and-profile-print output input-filename)))

(defun print-generation (generation generation-counts)
  (format t "~a: ~{~a~^ ~}~%" generation (coerce generation-counts 'list)))
(defun mortal-rabits (generations pair-lifetime &optional (pairs-per-litter 1) (initial-pairs 1) print)
  (let ((generation-counts (make-array (list pair-lifetime) :initial-element 0)))
    (setf (elt generation-counts 0) initial-pairs)
    (iter (repeat (1- generations))
	  (for generation from 0)
	  (when print
	    (print-generation generation generation-counts))
	  (let ((new-pair-count (* pairs-per-litter (iter (for gen-count-idx from 1 to (1- pair-lifetime))
							  (summing (elt generation-counts gen-count-idx))))))
	    (iter (for gen-count-idx from (1- pair-lifetime) downto 1)
		  (setf (elt generation-counts gen-count-idx) (elt generation-counts (1- gen-count-idx))))
	    (setf (elt generation-counts 0) new-pair-count)))
    (when print
      (print-generation generations generation-counts))
    (iter (for count in-vector generation-counts)
	  (sum count))))
(define-rosalind-problem :fibd rosalind-mortal-rabits
  "mortal fibonacci rabbits"
  (with-single-input-line (problem-data)
   (destructuring-bind-integers (generations pair-lifetime) problem-data
     (mortal-rabits generations pair-lifetime))))

(defun has-common-substr-at (dna-string-1 substr-start substr-length dna-string-2 position)
  (iter (for pos-1 from substr-start)
	(for pos-2 from position)
	(repeat substr-length)
	(always (char= (elt dna-string-1 pos-1) (elt dna-string-2 pos-2)))))
(defun has-common-substr (dna-string-1 substr-start substr-length dna-string-2)
  (iter (for pos-2 from 0 to (- (length dna-string-2) substr-length))
	(thereis (has-common-substr-at dna-string-1 substr-start substr-length dna-string-2 pos-2))))
(defun all-have-common-substr (dna-string-1 substr-start substr-length other-strings)
  (iter (for dna-string-2 in other-strings)
	(always (has-common-substr dna-string-1 substr-start substr-length dna-string-2))))
(define-rosalind-problem :lcsm rosalind-longest-common-substr
  "finding a shared motif"
  (with-fasta-input-lines (fasta-lines)
   (let* ((dna-strings (mapcar #'second fasta-lines))
	  (sorted-dna-strings (sort dna-strings #'< :key #'length))
	  (shortest-string (first sorted-dna-strings)))
     (iter outer (for substr-length from (length shortest-string) downto 2)
	   (iter (for position from 0 to (- (length shortest-string) substr-length))
		 (when (all-have-common-substr shortest-string position substr-length (rest sorted-dna-strings))
		   (return-from outer (subseq-of-length shortest-string position substr-length))))))))

(defun list-permutations (list)
  (if (not (null list))
    (iter (for i in list)
	  (appending (mapcar #'(lambda (r) (cons i r)) (list-permutations (remove i list)))))
    (list nil)))
(define-rosalind-problem :perm rosalind-permutations
  "enumerating gene orders"
  (with-single-input-line (l)
   (let* ((list-length (parse-integer l))
	  (table (iter (for i from 1 to list-length)
		       (collect i))))
     (with-output-to-file (output)
       (format output "~a~%" (reduce #'* (loop for i from 1 to list-length collect i)))
       (iter (for permutation in (list-permutations table))
	     (format output "~{~a~^ ~}~%" permutation))))))


(defun is-reverse-palindrome-subseq (dna start-pos length)
  (let* ((orig-substring (subseq-of-length dna start-pos length))
	 (revc (reverse-complement orig-substring)))
    (string= orig-substring revc)))
(defun find-restriction-sites (dna-string)
  (iter outer (for substring-length from 4 to 12)
	(iter (for position from 0 to (- (length dna-string) substring-length))
	      (when (is-reverse-palindrome-subseq dna-string position substring-length)
		(in outer (collect (list position substring-length)))))))
(define-rosalind-problem :revp restriction-sites
  "locating restriction sites"
  (with-single-fasta-line (dna-string)
    (with-output-to-file (output)
      (iter (for (pos length) in (sort (find-restriction-sites dna-string) #'< :key #'first))
	    (format output "~a ~a~%" (1+ pos) length)))))

(defun fill-nth-lex-kmer (kmer monomer-vector kmer-count)
  (let ((digits (length monomer-vector)))
    (iter (for digit from (1- (length kmer)) downto 0)
	  (setf (elt kmer digit) (elt monomer-vector (mod kmer-count digits)))
	  (setf kmer-count (floor kmer-count digits)))))
(defun do-for-all-kmers (monomer-vector k fun)
  "call fun for every k-mer produced with elements from monomer-vector (in lexicographical order)"
  (let ((kmer (make-array k :initial-element (elt monomer-vector 0))))
    (iter (for kmer-count from 0 to (1- (expt (length monomer-vector) k)))
	  (fill-nth-lex-kmer kmer monomer-vector kmer-count)
	  (funcall fun kmer))))
(define-rosalind-problem :lexf rosalind-enumerate-kmers-lex
  "enumerating k-mers lexicographically"
  (with-input-lines (lines)
      (let* ((monomer-vector (remove #\Space (first lines)))
	     (k (parse-integer (second lines))))
	(with-output-to-file (s)
	  (do-for-all-kmers monomer-vector k #'(lambda (v) (format s "~a~%" (coerce v 'string))))))))

(define-rosalind-problem :grph overlap-graphs
  "overlap graphs"
  (with-fasta-input-lines (fasta-data)
    (let ((entries-with-prefix (make-hash-table :test #'equal))
	  (entries-with-postfix (make-hash-table :test #'equal)))
      (iter (for (fasta-id sequence) in fasta-data)
	    (let ((prefix (subseq-of-length sequence 0 3))
		  (postfix (subseq-of-length sequence (- (length sequence) 3) 3)))
	      (push fasta-id (gethash prefix entries-with-prefix))
	      (push fasta-id (gethash postfix entries-with-postfix))))
      (with-output-to-file (output)
	(iter (for (prefix prefix-fasta-ids) in-hashtable entries-with-prefix)
	      (iter (for prefix-fasta-id in prefix-fasta-ids)
		    (iter (for postfix-fasta-id in (gethash prefix entries-with-postfix))
			  (when (not (string= prefix-fasta-id postfix-fasta-id))
			    (format output "~a ~a~%" postfix-fasta-id prefix-fasta-id)))))))))

(define-rosalind-problem :iev expected-offspring
  "calculating expected offspring"
  (with-single-input-line (problem-data)
   (let ((number-couples (integer-list problem-data))
	 (probability-for-dominant '(1 1 1 0.75 0.5 0)))
     (* 2 (iter (for k in number-couples)
		(for p in probability-for-dominant)
		(summing (* k p)))))))

(define-rosalind-problem :lia independent-alleles
  "independent-alleles"
  (with-single-input-line (problem-data)
   (destructuring-bind-integers (k n) problem-data
     (let ((profile (child-genotype-profile "AaBb" "AaBb")))
       (float (iter (for i from n to (expt 2 k))
		    (summing (probability-for-events profile "AaBb" i (expt 2 k)))))))))

(define-rosalind-problem :prtm ros-protein-mass
  "calculating protein mass"
  (with-single-input-line (problem-data)
    (protein-mass problem-data)))

(define-rosalind-problem :mrna ros-infer-rna-from-protein
  "infering mRNA from protein"  
  (with-single-input-line (protein-string)
    ;; account for 3 possible stop codons
    (mod (* 3 (count-possible-rna-sources-for-protein protein-string)) (expt 10 6))))

(define-rosalind-problem :splc ros-splice-rna
  "rna splicing"
  (with-fasta-input-lines (fasta-data)
   (let* ((dna (second (first fasta-data)))
	  (introns (iter (for (fasta-id string) in (rest fasta-data))
			 (collect string))))
     (iter (for intron in introns)
	   (setf dna (cl-ppcre:regex-replace-all intron dna "")))
     (format nil "~{~a~}" (butlast (dna-to-protein dna))))))

(defun get-uniprot-protein (protein-id)
  (second (car (get-uniprot-fasta-cached protein-id))))
(define-rosalind-problem :mprt ros-find-protein-motif
  "finding a protein motif"
  (with-input-lines (lines)
    (with-output-to-file (stream)
      (iter (for protein-id in lines)
	    (let ((protein (get-uniprot-protein protein-id)))
	      (alexandria:when-let ((motif-positions (protein-n-glycosilation-positions protein)))
		(format stream "~a~%" protein-id)
		(print-integer-list motif-positions stream)))))))

(define-rosalind-problem :pper ros-partial-permutations
  "partial permutations"
  (with-single-input-line (problem-data)
    (destructuring-bind-integers (n p) problem-data
      (mod (permutations n p) 1000000))))

(defun string-with-matching-subseq (dna-string start end dna-string-list)
  (iter (for other-dna-string in dna-string-list)
	(alexandria:when-let (pos (position-of-substring dna-string start end other-dna-string))
	  (return (cons other-dna-string pos)))))
(define-rosalind-problem :long ros-assemble-genome
  "genome assembly as shortest superstring"
  (with-fasta-input-lines (fasta-data)
    (iter (for dna-string on (mapcar #'second fasta-data))
	  (let ((prefix-start 0)
		(prefix-end (floor (length (car dna-string)) 2)))
	    (alexandria:if-let
		((match-data (string-with-matching-subseq (car dna-string) prefix-start prefix-end (cdr dna-string))))
	      )))))

(defun apply-signs (list p)
  "uses the bits of p to decide when to change the sign of elements in list"
  (iter (for i in list)
	(if (zerop (mod p 2))
	    (collect i)
	    (collect (- i)))
	(setf p (floor p 2))))
(defun list-signed-permutations (list)
  (let ((permutations (list-permutations list))
	(l (length list)))
    (iter outer (for permutation in permutations)
	  (iter (for p from 0 to (1- (expt 2 l)))
		(in outer (collect (apply-signs permutation p)))))))
(define-rosalind-problem :sign ros-gene-oriented-orderings
  "enumerating oriented gene orderings"
  (with-single-input-line (line)
    (let* ((n (parse-integer line))
	   (signed-permutations (list-signed-permutations (iter (for i from 1 to n) (collect i)))))
      (with-output-to-file (stream)
	(format stream "~a~%" (length signed-permutations))
	(iter (for perm in signed-permutations)
	      (format stream "~{~a~^ ~}~%" perm))))))


(defun has-kmer-at-pos-p (dna-string kmer start-pos)
  (let ((dna-string-length (length dna-string))
	(kmer-length (length kmer)))
    (and (<= (+ start-pos kmer-length) dna-string-length)
	 (iter (for pos from 0 below kmer-length)
	       (let ((dna-string-pos (+ start-pos pos)))
		 (always (char= (elt dna-string dna-string-pos)
				(elt kmer pos))))))))
(defun kmer-count (dna-string kmer)
  (iter (for pos from 0 to (1- (length dna-string)))
	(counting (has-kmer-at-pos-p dna-string kmer pos))))
(defun kmer-composition (source dna-string k)
  (let (result)
    (do-for-all-kmers
	source k #'(lambda (v) (push (kmer-count dna-string (coerce v 'string)) result)))
    (reverse result)))
(define-rosalind-problem :kmer ros-kmer-composition
  "k-mer composition"
  (with-single-fasta-line (dna-string)
    (let ((result (kmer-composition "ACGT" dna-string 4)))
      (with-output-to-file (stream)
	(print-integer-list result stream)))))

(defun probability-to-create-random-string-with-gc (gc dna-string)
  (iter (for letter in-vector dna-string)
	(multiply (if (or (char= #\A letter) (char= #\T letter))
		      (/ (- 1 gc) 2)
		      (/ gc 2)))))
(define-rosalind-problem :rstr ros-match-random-motif
  (with-input-lines (lines)
    (destructuring-bind (count-str gc-str) (split-sequence:split-sequence #\Space (first lines))
      (let ((*read-default-float-format* 'double-float))
	(let* ((count (parse-integer count-str))
	       (gc (* 1.0d0 (parse-number:parse-real-number gc-str)))
	       (string-prob (probability-to-create-random-string-with-gc gc (second lines))))
	  (format t "~f~%" (- 1 (expt (- 1 string-prob) count))))))))
