(in-package :bioinformatics)

(define-rosalind-problem :dna "rosalind_dna.txt" rosalind-dna
  "counting dna nucleotides"
  (let ((dna (first (read-file-lines input-filename)))
	(count-hash (make-hash-table)))
    (iter (for letter in-vector dna)
	  (incf (gethash letter count-hash 0)))
    (format nil "~{~a ~}~%"
	    (iter (for letter in-vector "ACGT")
		  (collect (gethash letter count-hash))))))

(define-rosalind-problem :rna "rosalind_rna.txt" rosalind-rna
  "transcribing dna into rna"
  (let ((dna (first (read-file-lines input-filename))))
    (cl-ppcre:regex-replace-all "T" dna "U")))

(define-rosalind-problem :revc "rosalind_revc.txt" rosalind-revc
  "complementing a strand of dna"
  (let ((dna (first (read-file-lines input-filename))))
    (reverse-complement dna)))

(defun ros-fib (generations pairs-per-litter &optional (fnm2 1) (fnm1 1))
  (cond ((= generations 0) fnm1)
	((= generations 1) fnm2)
	(t (ros-fib (1- generations) pairs-per-litter fnm1 (+ fnm1 (* pairs-per-litter fnm2))))))
(define-rosalind-problem :fib "rosalind_fib.txt" rosalind-fib
  "rabbits and recurrence relations"
  (destructuring-bind-integers  (generations pairs-per-litter) (first (read-file-lines input-filename))
    (ros-fib generations pairs-per-litter)))

(defun rosalind-gc-content (dna)
  (* 1.0 (/ (iter (for letter in-vector dna)
		  (counting (or (char= letter #\C) (char= letter #\G))))
	    (length dna))))
(define-rosalind-problem :gc "rosalind_gc.txt" rosalind-gc
  "computing gc content"
  (let ((lines (read-fasta-lines input-filename)))
    (iter (for (name dna) in lines)
	  (finding (list name (* 100 (rosalind-gc-content dna))) maximizing (rosalind-gc-content dna)))))

(defun hamming-distance (dna1 dna2)
  (iter (for letter1 in-vector dna1)
	(for letter2 in-vector dna2)
	(counting (not (char= letter1 letter2)))))
(define-rosalind-problem :hamm "rosalind_hamm.txt" rosalind-hamming
  "counting point mutations"
  (let* ((lines (read-file-lines input-filename))
	 (dna1 (first lines))
	 (dna2 (second lines)))
    (hamming-distance dna1 dna2)))

(defun rosalind-mendelian (homozygous-dominant-count heterozygous-count homozygous-recessive-count)
  (let ((k homozygous-dominant-count)
	(l heterozygous-count)
	(m homozygous-recessive-count))
    (/ (+ (* k (1- k)) (* 2 k l) (* 2 k m) (* m l)  (/ (* 3 l (1- l)) 4)) (* (+ k l m) (+ k l m -1)))))
(define-rosalind-problem :iprb "rosalind_iprb.txt" mendels-first-law
  "mendel's first law"
  (let ((lines (read-file-lines input-filename)))
    (destructuring-bind-integers (homozygous-dominant-count heterozygous-count homozygous-recessive-count) (first lines)
      (float (rosalind-mendelian homozygous-dominant-count heterozygous-count homozygous-recessive-count)))))

(define-rosalind-problem :prot "rosalind_prot.txt" rosalind-rna-to-protein
  "translating rna into protein"
  (let ((lines (read-file-lines input-filename)))
    (coerce (mapcar #'(lambda (s) (if (= 1 (length s))
				      (elt s 0)
				      ""))
		    (butlast (rna-to-protein (first lines))))
	    'string)))

(define-rosalind-problem :subs "rosalind_subs.txt" rosalind-find-motifs
  "finding a motif in dna"
  (let* ((lines (read-file-lines input-filename))
	 (haystack (first lines))
	 (needle (second lines)))
    (format t "~{~a~^ ~}" (mapcar #'1+ (find-pattern-occurences haystack needle)))))

(defun rosalind-consensus-and-profile (input-filename)
  (let* ((fasta-lines (read-fasta-lines input-filename))
	 (profile (profile-matrix (mapcar #'second fasta-lines))))    
    (list (consensus profile) profile)))
(defun rosalind-consensus-and-profile-print (stream input-filename)
  (destructuring-bind (consensus (a-prof t-prof c-prof g-prof)) (rosalind-consensus-and-profile input-filename)
    (format stream "~a~%" (coerce consensus 'string))
    (format stream "A: ~{~A~^ ~}~%" (coerce a-prof 'list))
    (format stream "C: ~{~A~^ ~}~%" (coerce c-prof 'list))
    (format stream "G: ~{~A~^ ~}~%" (coerce g-prof 'list))
    (format stream "T: ~{~A~^ ~}~%" (coerce t-prof 'list))))
(define-rosalind-problem :cons "rosalind_cons.txt" rosalind-consensus-and-profile-to-file
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
(define-rosalind-problem :fibd "rosalind_fibd.txt" rosalind-mortal-rabits
  "mortal fibonacci rabbits"
  (destructuring-bind-integers (generations pair-lifetime)
      (first (read-file-lines input-filename))
    (mortal-rabits generations pair-lifetime)))

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
(define-rosalind-problem :lcsm "rosalind_lcsm.txt" rosalind-longest-common-substr
  "finding a shared motif"
    (let* ((dna-strings (mapcar #'second (read-fasta-lines input-filename)))
	 (sorted-dna-strings (sort dna-strings #'< :key #'length))
	 (shortest-string (first sorted-dna-strings)))
    (iter outer (for substr-length from (length shortest-string) downto 2)
	  (iter (for position from 0 to (- (length shortest-string) substr-length))
		(when (all-have-common-substr shortest-string position substr-length (rest sorted-dna-strings))
		  (return-from outer (subseq-of-length shortest-string position substr-length)))))))

(defun list-permutations (list)
  (if (not (null list))
    (iter (for i in list)
	  (appending (mapcar #'(lambda (r) (cons i r)) (list-permutations (remove i list)))))
    (list nil)))
(define-rosalind-problem :perm "rosalind_perm.txt" rosalind-permutations
  "enumerating gene orders"
  (let* ((list-length (parse-integer (first (read-file-lines input-filename))))
	 (table (iter (for i from 1 to list-length)
		      (collect i))))
    (with-output-to-file (output)
     (format output "~a~%" (reduce #'* (loop for i from 1 to list-length collect i)))
     (iter (for permutation in (list-permutations table))
	   (format output "~{~a~^ ~}~%" permutation)))))


(defun is-reverse-palindrome-subseq (dna start-pos length)
  (let* ((orig-substring (subseq-of-length dna start-pos length))
	 (revc (reverse-complement orig-substring)))
    (string= orig-substring revc)))
(defun find-restriction-sites (dna-string)
  (iter outer (for substring-length from 4 to 12)
	(iter (for position from 0 to (- (length dna-string) substring-length))
	      (when (is-reverse-palindrome-subseq dna-string position substring-length)
		(in outer (collect (list position substring-length)))))))
(define-rosalind-problem :revp "rosalind_revp.txt" restriction-sites
  "locating restriction sites"
  (let ((dna-string (second (first (read-fasta-lines input-filename)))))
    (with-output-to-file (output)
      (iter (for (pos length) in (sort (find-restriction-sites dna-string) #'< :key #'first))
	    (format output "~a ~a~%" (1+ pos) length)))))

(defun fill-nth-lex-kmer (kmer monomer-vector kmer-count)
  (let ((digits (length monomer-vector)))
    (iter (for digit from (1- (length kmer)) downto 0)
	  (setf (elt kmer digit) (elt monomer-vector (mod kmer-count digits)))
	  (setf kmer-count (floor kmer-count digits)))))
(defun do-for-all-kmers (monomer-vector k fun)
  (let ((kmer (make-array k :initial-element (elt monomer-vector 0))))
    (iter (for kmer-count from 0 to (1- (expt (length monomer-vector) k)))
	  (fill-nth-lex-kmer kmer monomer-vector kmer-count)
	  (funcall fun kmer))))
(define-rosalind-problem :lexf "rosalind_lexf.txt" rosalind-enumerate-kmers-lex
  "enumerating k-mers lexicographically"
  (let* ((lines (read-file-lines input-filename))
	 (monomer-vector (remove #\Space (first lines)))
	 (k (parse-integer (second lines))))
    (with-output-to-file (s)
      (do-for-all-kmers monomer-vector k #'(lambda (v) (format s "~a~%" (coerce v 'string)))))))

(define-rosalind-problem :grph "rosalind_grph.txt" overlap-graphs
  "overlap graphs"
  (let ((fasta-data (read-fasta-lines input-filename))
	(entries-with-prefix (make-hash-table :test #'equal))
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
			 (format output "~a ~a~%" postfix-fasta-id prefix-fasta-id))))))))

(define-rosalind-problem :iev "rosalind_iev.txt" expected-offspring
  "calculating expected offspring"
  (let ((number-couples (integer-list (first (read-file-lines input-filename))))
	(probability-for-dominant '(1 1 1 0.75 0.5 0)))
    (* 2 (iter (for k in number-couples)
	       (for p in probability-for-dominant)
	       (summing (* k p))))))
