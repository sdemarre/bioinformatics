(in-package :bioinformatics)

(define-rosalind-problem :dna
  "counting dna nucleotides"
  (with-single-input-line (dna)
   (let ((count-hash (make-hash-table)))
     (iter (for letter in-vector dna)
	   (incf (gethash letter count-hash 0)))
     (format nil "~{~a ~}~%"
	     (iter (for letter in-vector "ACGT")
		   (collect (gethash letter count-hash)))))))

(define-rosalind-problem :rna
  "transcribing dna into rna"
  (with-single-input-line (dna)
    (cl-ppcre:regex-replace-all "T" dna "U")))

(define-rosalind-problem :revc
  "complementing a strand of dna"
  (with-single-input-line (dna)
    (reverse-complement dna)))

(defun ros-fib (generations pairs-per-litter &optional (fnm2 1) (fnm1 1))
  (cond ((= generations 0) fnm1)
	((= generations 1) fnm2)
	(t (ros-fib (1- generations) pairs-per-litter fnm1 (+ fnm1 (* pairs-per-litter fnm2))))))
(define-rosalind-problem :fib
  "rabbits and recurrence relations"
  (with-single-input-line (problem-data)
    (destructuring-bind-integers  (generations pairs-per-litter) problem-data
      (ros-fib generations pairs-per-litter))))

(defun gc-content (dna)
  (* 1.0 (/ (iter (for letter in-vector dna)
		  (counting (or (char= letter #\C) (char= letter #\G))))
	    (length dna))))
(define-rosalind-problem :gc
  "computing gc content"
  (with-fasta-input-lines (lines)
    (iter (for (name dna) in lines)
	  (finding (list name (* 100 (gc-content dna))) maximizing (gc-content dna)))))

(defun hamming-distance (dna1 dna2)
  (iter (for letter1 in-vector dna1)
	(for letter2 in-vector dna2)
	(counting (not (char= letter1 letter2)))))
(define-rosalind-problem :hamm
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
(define-rosalind-problem :iprb
  "mendel's first law"
  (with-input-lines (lines)
    (destructuring-bind-integers (homozygous-dominant-count heterozygous-count homozygous-recessive-count) (first lines)
      (float (rosalind-mendelian homozygous-dominant-count heterozygous-count homozygous-recessive-count)))))

(define-rosalind-problem :prot
  "translating rna into protein"
  (with-input-lines (lines)
    (coerce (mapcar #'(lambda (s) (if (= 1 (length s))
				      (elt s 0)
				      ""))
		    (butlast (rna-to-protein (first lines))))
	    'string)))

(define-rosalind-problem :subs
  "finding a motif in dna"
  (with-input-lines (lines)
    (let* ((haystack (first lines))
	   (needle (second lines)))
      (with-output-to-file (stream)
	  (format stream "~{~a~^ ~}" (mapcar #'1+ (find-pattern-occurences haystack needle)))))))

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
(define-rosalind-problem :cons
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
(define-rosalind-problem :fibd
  "mortal fibonacci rabbits"
  (with-single-input-line (problem-data)
   (destructuring-bind-integers (generations pair-lifetime) problem-data
     (mortal-rabits generations pair-lifetime))))

(defun has-common-substr-at (dna-string-1 substr-start substr-length dna-string-2 position)
  (and (<= (+ substr-start substr-length) (length dna-string-1))
       (<= (+ position substr-length) (length dna-string-2))
       (>= substr-start 0)
       (>= position 0)
       (iter (for pos-1 from substr-start)
	     (for pos-2 from position)
	     (repeat substr-length)
	     (always (char= (elt dna-string-1 pos-1) (elt dna-string-2 pos-2))))))
(defun has-common-substr (dna-string-1 substr-start substr-length dna-string-2)
  (iter (for pos-2 from 0 to (- (length dna-string-2) substr-length))
	(thereis (has-common-substr-at dna-string-1 substr-start substr-length dna-string-2 pos-2))))
(defun all-have-common-substr (dna-string-1 substr-start substr-length other-strings)
  (iter (for dna-string-2 in other-strings)
	(always (has-common-substr dna-string-1 substr-start substr-length dna-string-2))))
(define-rosalind-problem :lcsm
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
(define-rosalind-problem :perm
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
(define-rosalind-problem :revp
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
  (let ((kmer (coerce (make-array k :initial-element (elt monomer-vector 0)) 'string)))
    (iter (for kmer-count from 0 to (1- (expt (length monomer-vector) k)))
	  (fill-nth-lex-kmer kmer monomer-vector kmer-count)
	  (funcall fun kmer))))
(define-rosalind-problem :lexf
  "enumerating k-mers lexicographically"
  (with-input-lines (lines)
      (let* ((monomer-vector (remove #\Space (first lines)))
	     (k (parse-integer (second lines))))
	(with-output-to-file (s)
	  (do-for-all-kmers monomer-vector k #'(lambda (v) (format s "~a~%" v)))))))

(define-rosalind-problem :grph
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

(define-rosalind-problem :iev
  "calculating expected offspring"
  (with-single-input-line (problem-data)
   (let ((number-couples (parse-integer-list problem-data))
	 (probability-for-dominant '(1 1 1 0.75 0.5 0)))
     (* 2 (iter (for k in number-couples)
		(for p in probability-for-dominant)
		(summing (* k p)))))))

(define-rosalind-problem :lia
  "independent-alleles"
  (with-single-input-line (problem-data)
   (destructuring-bind-integers (k n) problem-data
     (let ((profile (child-genotype-profile "AaBb" "AaBb")))
       (float (iter (for i from n to (expt 2 k))
		    (summing (probability-for-events profile "AaBb" i (expt 2 k)))))))))

(define-rosalind-problem :prtm
  "calculating protein mass"
  (with-single-input-line (problem-data)
    (protein-mass problem-data)))

(define-rosalind-problem :mrna
  "infering mRNA from protein"  
  (with-single-input-line (protein-string)
    ;; account for 3 possible stop codons
    (mod (* 3 (count-possible-rna-sources-for-protein protein-string)) (expt 10 6))))

(define-rosalind-problem :splc
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
(define-rosalind-problem :mprt
  "finding a protein motif"
  (with-input-lines (lines)
    (with-output-to-file (stream)
      (iter (for protein-id in lines)
	    (let ((protein (get-uniprot-protein protein-id)))
	      (alexandria:when-let ((motif-positions (protein-n-glycosilation-positions protein)))
		(format stream "~a~%" protein-id)
		(print-integer-list motif-positions stream)))))))

(define-rosalind-problem :pper
  "partial permutations"
  (with-single-input-line (problem-data)
    (destructuring-bind-integers (n p) problem-data
      (mod (permutations n p) 1000000))))

(defun string-with-matching-subseq (dna-string start end dna-string-list)
  (iter (for other-dna-string in dna-string-list)
	(alexandria:when-let (pos (position-of-substring dna-string start end other-dna-string))
	  (return (cons other-dna-string pos)))))


(defun log-prob-of-string-given-gc (dna-string gc-content)
  (iter (for letter in-vector dna-string)
	(summing (log (if (or (char= letter #\A) (char= letter #\T))
			  (/ (- 1 gc-content) 2)
			  (/ gc-content 2))
		      10))))
(define-rosalind-problem :prob
  "introduction to random strings"
  (with-input-lines (lines)
    (let* ((dna-string (first lines))
	   (gc-content-values (parse-float-list (second lines))))
      (let ((log-probs (iter (for gc-content in gc-content-values)
			     (collect (log-prob-of-string-given-gc dna-string gc-content)))))
       (with-output-to-file (stream)
	 (print-float-list log-probs stream))))))

(define-rosalind-problem :long-wrong
  "genome assembly as shortest superstring"
  (with-fasta-dna-lines (dna-strings)
    (let ((assembly-graph (make-instance 'graph :type :directed))
	  (nodes (make-array (length dna-strings))))
      (flet ((add-prefix-relation (s1 s2 idx1 idx2 substr-len substr-pos)
	       (flet ((get-node (idx s)
			(alexandria:if-let (node (elt nodes idx))
			  node
			  (let ((node (add-node assembly-graph idx)))
			    (setf (elt nodes idx) node)
			    (set-property node :dna s)
			    node))))
		 (let* ((from-node (get-node idx1 s1))
			(to-node (get-node idx2 s2))
			(edge (add-edge assembly-graph from-node to-node)))
		   (set-property edge :len substr-len)
		   (set-property edge :pos substr-pos)))))
	(iter (for dna-string-1 in dna-strings)
	      (for from-idx from 1)
	      (let* ((dna-string-1-half (floor (length dna-string-1) 2)))
		(iter (for dna-string-2 in dna-strings)
		      (for to-idx from 1)
		      (unless (eq dna-string-1 dna-string-2)
			(iter (for substr-len from (length dna-string-1) downto dna-string-1-half)
			      (for dna-string-2-pos from 0)
			      (when (has-common-substr-at dna-string-1 0 substr-len dna-string-2 dna-string-2-pos)
				(add-prefix-relation dna-string-1 dna-string-2 from-idx to-idx substr-len dna-string-2-pos)
				(return))))))))
      assembly-graph)))

(defun compute-reassembly-origin (best-prefix-matches count-strings)
  (get-a-set-element
   (rset-difference (make-full-integer-set count-strings)
		    (make-set-from-list (iter (for (k v) in-hashtable best-prefix-matches)
					      (collect (second v)))))))
(defun reassemble-string (dna-strings best-prefix-matches)
  (let* ((dna-strings (coerce dna-strings 'vector))
	 (current-string-index (compute-reassembly-origin best-prefix-matches (length dna-strings))))
    (iter (while (gethash current-string-index best-prefix-matches))
	  (let ((current-string (elt dna-strings (1- current-string-index)))
		(new-string-index (second (gethash current-string-index best-prefix-matches))))
	    (if (gethash new-string-index best-prefix-matches)
		(collect (subseq-of-length current-string 0 (- (length current-string) (third (gethash current-string-index best-prefix-matches)))))
		(collect (elt dna-strings (1- new-string-index))))
	    (setf current-string-index new-string-index)))))
(defun matching-prefix-postfix-p (string1 prefix-len string2)
  (has-common-substr-at string1 0 prefix-len string2 (- (length string2) prefix-len)))
(defun compute-best-prefix-matches (dna-strings)
  (let ((best-prefix-match (make-hash-table)))
    (iter (for dna-string-1 in dna-strings)
	  (for idx1 from 1)
	  (iter (for dna-string-2 in dna-strings)
		(for idx2 from 1)
		(unless (eq dna-string-1 dna-string-2)
		  (iter (for substr-len from (length dna-string-1) downto (floor (length dna-string-1) 2))
			(when (matching-prefix-postfix-p dna-string-1 substr-len dna-string-2)
			  (setf (gethash idx2 best-prefix-match) (list idx2 idx1 substr-len dna-string-2 dna-string-1))
			  (return))))))
    best-prefix-match))
(define-rosalind-problem :long
  "genome assembly as shortest superstring"
  (with-fasta-dna-lines (dna-strings)
    (let ((best-prefixes (compute-best-prefix-matches dna-strings)))
      (with-output-to-file (s)
	(format s "~{~a~}" (reassemble-string dna-strings best-prefixes))))))

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
(define-rosalind-problem :sign
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
(defun kmer-composition (monomers dna-string k)
  (let (result)
    (do-for-all-kmers
	monomers k #'(lambda (kmer) (push (kmer-count dna-string kmer) result)))
    (reverse result)))
(define-rosalind-problem :kmer
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
(define-rosalind-problem :rstr
  (with-input-lines (lines)
    (destructuring-bind (count-str gc-str) (split-sequence:split-sequence #\Space (first lines))
      (let ((*read-default-float-format* 'double-float))
	(let* ((count (parse-integer count-str))
	       (gc (* 1.0d0 (parse-number:parse-real-number gc-str)))
	       (string-prob (probability-to-create-random-string-with-gc gc (second lines))))
	  (with-output-to-file (stream)
	    (format stream "~f~%" (- 1 (expt (- 1 string-prob) count)))))))))

(defun enumerate-words (letters max-word-length fun &optional current-word (current-word-length 0))
  (unless (zerop current-word-length)
    (funcall fun current-word))
  (unless (= current-word-length max-word-length)
    (iter (for letter in letters)
	  (enumerate-words letters max-word-length fun (cons letter current-word) (1+ current-word-length)))))
(define-rosalind-problem :lexv
  "ordering strings ov varying length lexicographically"
  (with-input-lines (lines)
    (let ((letters (split-sequence:split-sequence #\Space (first lines)))
	  (max-word-length (parse-integer (second lines))))
      (with-output-to-file (s)
	(enumerate-words letters max-word-length
			 #'(lambda (w) (format s "~{~a~}~%" (reverse w))))))))


(defun transition-transversion-ratio (dna-string1 dna-string2)
  (let ((transitions 0)
	(transversions 0))
    (iter (for a1 in-vector dna-string1)
	  (for a2 in-vector dna-string2)
	  (unless (char= a1 a2)
	    (cond ((or (char= a1 #\A) (char= a1 #\G))
		   (if (or (char= a2 #\T) (char= a2 #\C))
		       (incf transversions)
		       (incf transitions)))
		  (t (if (or (char= a2 #\T) (char= a2 #\C))
			 (incf transitions)
			 (incf transversions))))))
    (/ transitions transversions)))
(define-rosalind-problem :tran
  "transitions and transversions"
  (with-fasta-input-lines (fasta-data)
    (let ((s1 (second (first fasta-data)))
	  (s2 (second (second fasta-data))))
      (with-output-to-file (s)
	(format s "~f~%" (transition-transversion-ratio s1 s2))))))

(defun distance-matrix (dna-strings)
  (iter (for s1 in dna-strings)
	(collect (iter (for s2 in dna-strings)
		       (collect (/ (hamming-distance s1 s2) (length s1)))))))
(define-rosalind-problem :pdst
  "creating a distance matrix"
  (with-fasta-dna-lines (dna-strings)
    (let ((matrix (distance-matrix dna-strings)))
      (with-output-to-file (s)
	(format s "~{~{~f~^ ~}~%~}" matrix)))))

(define-rosalind-problem :sset
  "counting subsets"
  (with-single-input-line (elements)
    (with-output-to-file (s)
      (format s "~a~%" (mod (expt 2 (parse-integer elements)) 1000000)))))

(define-rosalind-problem :seto
  "introduction to set operations"
  (with-input-lines (lines)
    (let* ((set-size (parse-integer (first lines)))
	   (s1 (parse-rosalind-set set-size (second lines)))
	   (s2 (parse-rosalind-set set-size (third lines))))
      (with-output-to-file (stream)
	(between-forms (terpri stream)	  
	  (print-set (rset-union s1 s2) stream)
	  (print-set (rset-intersection s1 s2) stream)
	  (print-set (rset-difference s1 s2) stream)
	  (print-set (rset-difference s2 s1) stream)
	  (print-set (rset-complement s1) stream)
	  (print-set (rset-complement s2) stream))))))

(defun spliced-motif-positions (dna-string motif)
  (let ((target-position 0))
    (iter (for letter in-vector motif)
	  (collect (setf target-position (position letter dna-string :start target-position)))
	  (incf target-position))))
(defun verify-spliced-positions (spliced-positions dna-string motif)
  (assert (= (length spliced-positions) (length motif)))
  (assert (equalp (remove-duplicates spliced-positions) spliced-positions))
  (iter (for pos in spliced-positions)
	(for motif-pos from 0)
	(assert (char= (elt dna-string pos) (elt motif motif-pos)))))
(define-rosalind-problem :sseq
  "finding a spliced motif"
  (with-fasta-dna-lines (dna-strings)
    (let ((dna-string (first dna-strings))
	  (motif (second dna-strings)))
     (let ((spliced-positions (spliced-motif-positions dna-string motif)))
       (with-output-to-file (s)
	 (print-integer-list (mapcar #'1+ spliced-positions) s))))))


