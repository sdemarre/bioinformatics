(in-package :bioinformatics)

(defun count-insertion-sort-swaps (data)
  (let ((swaps-count 0))
    (iter (for i from 1 to (1- (length data)))
	  (iter (for k from i downto 1)
		(symbol-macrolet ((data-k (elt data k))
				  (data-k-1 (elt data (1- k))))
		  (while (< data-k data-k-1))
		  (psetf data-k data-k-1
			 data-k-1 data-k)
		  (incf swaps-count))))
    swaps-count))
(define-rosalind-problem :ins insertion-sort-swaps
  "insertion sort"
  (let* ((data (second (read-file-lines input-filename)))
	 (unsorted-data (parse-integer-list data)))
    (count-insertion-sort-swaps (coerce unsorted-data 'vector))))

(defun binary-search (sequence compare)
  (let* ((first 0)
	 (last (length sequence))
	 (count (- last first)))
    (iter (while (> count 0))
	  (let* ((count2 (floor count 2))
		 (mid (+ first count2)))
	    (if (funcall compare (elt sequence mid))
		(progn
		  (setf first (incf mid))
		  (decf count (1+ count2)))
		(setf count count2))))
    first))
(defun position-binary-search (element sorted-array compare)
  (let ((pos (binary-search sorted-array #'(lambda (e) (funcall compare e element)))))
    (when (and (< pos (length sorted-array))
	       (= (elt sorted-array pos) element))
      pos)))
(define-rosalind-problem :bins rosalind-binary-search
  "binary search"
  (let* ((lines (read-file-lines input-filename))
	 (sorted-array (coerce (parse-integer-list (third lines)) 'vector))
	 (elements-to-find (parse-integer-list (fourth lines))))
    (with-output-to-file (stream)
      (print-integer-list
	      (iter (for element in elements-to-find)
		    (collect (let ((pos (position-binary-search element sorted-array #'<)))
		      	       (if pos (1+ pos) -1))))
	      stream))))

(defun make-graph-from-file (filename &optional (type :undirected))
  (let* ((lines (read-file-lines filename))
	 (graph (make-instance 'graph :type type)))
    (destructuring-bind-integers (count-nodes count-edges) (first lines)
      (declare (ignorable count-edges))
      (let ((nodes (make-array (1+ count-nodes))))
	(iter (for node-value from 1 to count-nodes)
	      (setf (elt nodes node-value) (add-node graph node-value)))
	(iter (for line in (rest lines))
	      (destructuring-bind-integers (source-value target-value) line
		(add-edge graph (elt nodes source-value) (elt nodes target-value))))))
    graph))
(define-rosalind-problem :deg degree-array
  "degree array"
  (let ((graph (make-graph-from-file input-filename)))
    (let ((degree-data
	   (iter (for node node-of-graph graph)
		 (collect (cons (value node) (degree node))))))
      (with-output-to-file (stream)
	(print-integer-list
	 (mapcar #'cdr (sort degree-data #'< :key #'car)) stream)))))

(define-rosalind-problem :ddeg double-degree-array
  "double-degree array"
  (let* ((graph (make-graph-from-file input-filename)))
    (let ((double-degree-data
	   (iter (for node node-of-graph graph)
		 (collect (cons (value node)
				(iter (for neighbour neighbour-of-node node)
				      (summing (degree neighbour))))))))
      (with-output-to-file (stream)
	(print-integer-list
	 (mapcar #'cdr (sort double-degree-data #'< :key #'car)) stream)))))

(defun shortest-distances (graph)
  (let ((nodes-distances (make-hash-table))
	(reachable-nodes (list (find-node graph 1))))
    (iter (while (not (null reachable-nodes)))
	  (for distance from 0)
	  (let (new-nodes)
	    (iter (for node in reachable-nodes)
		  (setf (gethash node nodes-distances) distance))
	    (iter (for node in reachable-nodes)
		  (iter (for neighbour neighbour-of-node node)
			(unless (gethash neighbour nodes-distances)
			  (pushnew neighbour new-nodes))))	    
	    (setf reachable-nodes new-nodes)))
    (let ((value-to-node (value-to-node-hash graph)))
      (iter (for value from 1 to (count-nodes graph))
	    (let ((node (gethash value value-to-node)))
	     (collect (alexandria:if-let ((distance (gethash node nodes-distances)))
			distance
			-1)))))))
(define-rosalind-problem :bfs breadth-first-search
  "breadth first search"
  (let* ((graph (make-graph-from-file input-filename :directed)))
    (with-output-to-file (stream)
      (print-integer-list (shortest-distances graph) stream))))


(defun majority-element (int-list)
  (let* ((counts (make-hash-table))
	 (list-length 0)
	 (int-with-max-count
	  (iter (for int in int-list)
		(incf list-length)
		(finding int maximizing (incf (gethash int counts 0))))))
    (if (> (gethash int-with-max-count counts) (/ list-length 2))
	int-with-max-count
	-1)))
(define-rosalind-problem :maj ros-majority-element
  "majority element"
  (iter (for line in (rest (read-file-lines input-filename)))
	(collect (majority-element (parse-integer-list line)))))

(define-rosalind-problem :mer ros-merge-sorted
  "merge sorted arrays"
  (let* ((lines (read-file-lines input-filename))
	 (list-1 (parse-integer-list (second lines)))
	 (list-2 (parse-integer-list (fourth lines))))
    (with-output-to-file (stream)
      (format stream "~{~a~^ ~}~%"
	      (iter (until (and (null list-1) (null list-2)))
		    (cond ((null list-1) (collect (pop list-2)))
			  ((null list-2) (collect (pop list-1)))
			  ((< (car list-1) (car list-2)) (collect (pop list-1)))
			  (t (collect (pop list-2)))))))))

(define-rosalind-problem :2sum ros-2sum
  "2sum"
  (let* ((lines (read-file-lines input-filename)))
    (with-output-to-file (stream)
     (iter (for line in (rest lines))
	   (for line-num from 1)
	   (let ((labeled-list (iter (for i from 1)
				     (for e in (parse-integer-list line))
				     (collect (cons e i)))))
	     (setf labeled-list (sort  labeled-list #'< :key #'(lambda (e) (abs (car e)))))
	     (let (found)
	       (iter (for e on labeled-list)
		     (while (not found))
		     (when (and (rest e)
				(= (- (car (first e))) (car (second e))))
		       (format stream "~a ~a~%" (cdr (first e)) (cdr (second e)))
		       (setf found t)))
	       (unless found
		 (format stream "-1~%"))))))))


(defun merge-sort (data)
  ;; cheating...
  (sort data #'<))
(define-rosalind-problem :ms ros-merge-sort
  "merge sort"
  (with-input-lines (lines)
   (with-output-to-file (stream)     
     (format stream "~{~a~^ ~}~%" (merge-sort (parse-integer-list (second lines)))))))


(defun transcode-open-frames (dna-string)
  (let ((start-positions (cl-ppcre:all-matches "ATG" dna-string))
	(dna-string-length (length dna-string))
	all-results)
    (iter (for start-pos on start-positions by #'cddr)
	  (let (stop
		result)
	    (iter (for codon-pos from (car start-pos) by 3)
		  (until (or  stop (> codon-pos (- dna-string-length 3))))
		  (let ((codon (subseq-of-length dna-string codon-pos 3)))
		    (if (is-dna-stop-codon-p codon)
			(progn
			  (setf stop t)
			  (unless (null result)
			    (push (coerce (nreverse result) 'string) all-results)))
			(push (elt (dna-codon-to-amino-acid codon) 0) result))))))
    all-results))
(define-rosalind-problem :orf ros-open-reading-frames
  "open reading frames"
  (with-single-fasta-line (dna-string)
    (let ((proteins (make-hash-table :test #'equal)))
      (iter (for protein in (transcode-open-frames dna-string))
	    (setf (gethash protein proteins) 1))
      (iter (for protein in (transcode-open-frames (reverse-complement dna-string)))
	    (setf (gethash protein proteins) 1))
      (with-output-to-file (stream)
	(format stream "~{~a~%~}" (alexandria:hash-table-keys proteins))))))

(defun longest-increasing-subseq-new (seq)
  (let* ((seq-length (length seq))
	 (path (make-array seq-length))
	 (max-elem (make-array (1+ seq-length)))
	 (current-best-end 0))
    (iter (for i index-of-vector seq)
	  )))
(defun longest-increasing-subseq (x)
  "as found on wikipedia"
  (let* ((n (length x))
	 (p (make-array n))
	 (m (make-array (1+ n)))
	 (l 0)
	 (new-l 0))
    (iter (for i index-of-vector x)
	  (let ((lo 1)
		(hi l)
		(mid 0))
	    (iter (while (<= lo hi))
		  (setf mid (ceiling (+ hi lo) 2))
		  (if (< (aref x (aref m mid)) (aref x i))
		      (setf lo (1+ mid))
		      (setf hi (1- mid))))
	    (setf new-l lo)
	    (setf (aref p i) (aref m (1- new-l)))
	    (setf (aref m new-l) i)
	    (setf l (max l new-l))))
    (reverse
     (let ((k (aref m l)))
       (iter (for i from (1- L) downto 0)
	     (collect (aref x k))
	     (setf k (aref p k)))))))
(define-rosalind-problem :lgis ros-longest-increasing-subseq
  "longest increasing subsequence"
  (with-input-lines (lines)
    (let* ((permutation (coerce (parse-integer-list (second lines)) 'vector)))
      (with-output-to-file (stream)
	(format stream "~{~a~^ ~}~%" (longest-increasing-subseq permutation))
	(format stream "~{~a~^ ~}~%" (reverse (longest-increasing-subseq (reverse permutation))))))))

(defun reachable-nodes-from (node)
  (let ((reachable-nodes (make-instance 'rosalind-set))
	new-nodes-reachable)
    (add-element reachable-nodes node)
    (iter (setf new-nodes-reachable nil)
	  (do-for-set-elements reachable-nodes
	    #'(lambda (node)
		(iter (for neighbour neighbour-of-node node)
		      (unless (has-element-p reachable-nodes neighbour)
			(setf new-nodes-reachable t)
			(add-element reachable-nodes neighbour)))))
	  (until (not new-nodes-reachable)))
    reachable-nodes))
(defun count-connected-components (graph)
  (let ((unprocessed-nodes (make-instance 'rosalind-set))
	(connected-components 0))
    (iter (for node node-of-graph graph)
	  (add-element unprocessed-nodes node))
    (iter (while (not (empty-set-p unprocessed-nodes)))
	  (let ((node (get-a-set-element unprocessed-nodes)))
	    (do-for-set-elements (reachable-nodes-from node)
	      #'(lambda (reached-node)
		  (remove-element unprocessed-nodes reached-node)))
	    (incf connected-components)))
    connected-components))
(define-rosalind-problem :cc ros-connected-components
  "connected components"
  (let ((graph (make-graph-from-file input-filename)))
    (with-output-to-file (stream)
      (format stream "~a~%" (count-connected-components graph)))))

(define-rosalind-problem :pmch ros-perfect-matchings
  "Perfect matchings and RNA secondary structures"
  (with-single-fasta-line (rna)
    (* (fact (count #\A rna)) (fact (count #\G rna)))))

(defun make-heap (data pred)
  (macrolet ((element (idx) `(elt data ,idx)))
    (labels ((parent-idx (idx) (floor (1- idx) 2))
	     (has-heap-property-p (idx) (funcall pred (element (parent-idx idx)) (element idx)))
	     (bubble-single-element (idx)
	       (psetf (element idx) (element (parent-idx idx))
		      (element (parent-idx idx)) (element idx)))
	     (bubble-up (idx)
	       (iter (while (not (zerop idx)))
		     (unless (has-heap-property-p idx)
		       (bubble-single-element idx))
		     (setf idx (parent-idx idx)))))
      (iter (for idx from 1 to (1- (length data)))
	    (bubble-up idx)))))
(define-rosalind-problem :hea ros-build-heap
  "building a heap"
  (with-input-lines (lines)
    (let ((data (parse-integer-vector (second lines))))
      (make-heap data #'>=)
      (with-output-to-file (s)
	(print-integer-vector data s)))))
