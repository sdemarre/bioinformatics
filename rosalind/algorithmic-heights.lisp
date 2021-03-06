(in-package :bioinformatics)

(define-rosalind-problem :fibo
    "fibonacci numbers"
  (with-single-input-line (line)
    (write-single-output-line (fib (parse-integer line)))))
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
(define-rosalind-problem :ins
  "insertion sort"
  (with-input-lines (lines)
    (let* ((data (second lines))
	   (unsorted-data (parse-integer-list data)))
      (write-single-output-line
	(count-insertion-sort-swaps (coerce unsorted-data 'vector))))))

(defun lower-bound (sequence begin end order element &key (key #'identity))
  (let* ((first begin)
	 (last end)
	 (count (- last first)))
    (iter (while (> count 0))
	  (let* ((count2 (floor count 2))
		 (mid (+ first count2)))
	    (if (funcall order (funcall key (elt sequence mid)) element)
		(progn
		  (setf first (incf mid))
		  (decf count (1+ count2)))
		(setf count count2))))
    first))
(defun position-binary-search (element sorted-array order &key (key #'identity) (begin 0) (end (length sorted-array)))
  "returns the position of element in sorted-array, which is sorted according to order, or nil if not found.
test checks if 2 elements are identical, after applying <key>.
searches through range [begin, end)"
  (let ((pos (lower-bound sorted-array  begin end order element :key key)))
    (when (and (< pos (length sorted-array))
	       (and (not (funcall order (funcall key (elt sorted-array pos)) element))
		    (not (funcall order element (funcall key (elt sorted-array pos))))))
      pos)))
(define-rosalind-problem :bins
    "binary search"
  (with-input-lines (lines)
   (let* ((sorted-array (coerce (parse-integer-list (third lines)) 'vector))
	  (elements-to-find (parse-integer-list (fourth lines))))
     (with-output-to-file (stream)
       (print-integer-list
	(iter (for element in elements-to-find)
	      (collect (let ((pos (position-binary-search element sorted-array #'<)))
			 (if pos (1+ pos) -1))))
	stream)))))

(defun make-graph-from-lines (lines &optional (type :undirected))
  (let ((graph (make-instance 'graph :type type)))
    (destructuring-bind-integers (count-nodes &optional count-edges) (first lines)
      (declare (ignorable count-edges))
      (let ((nodes (make-array (1+ count-nodes))))
	(iter (for node-value from 1 to count-nodes)
	      (setf (elt nodes node-value) (add-node graph node-value)))
	(iter (for line in (rest lines))
	      (destructuring-bind-integers (source-value target-value) line
		(add-edge graph (elt nodes source-value) (elt nodes target-value))))))
    graph))
(defun make-graph-from-file (filename &optional (type :undirected))
  (let* ((lines (read-file-lines filename)))
    (make-graph-from-lines lines type)))
(define-rosalind-problem :deg
  "degree array"
  (let ((graph (make-graph-from-file *input-filename*)))
    (let ((degree-data
	   (iter (for node node-of-graph graph)
		 (collect (cons (value node) (degree node))))))
      (with-output-to-file (stream)
	(print-integer-list
	 (mapcar #'cdr (sort degree-data #'< :key #'car)) stream)))))

(define-rosalind-problem :ddeg
  "double-degree array"
  (let* ((graph (make-graph-from-file *input-filename*)))
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
(define-rosalind-problem :bfs
  "breadth first search"
  (let* ((graph (make-graph-from-file *input-filename* :directed)))
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
(define-rosalind-problem :maj
  "majority element"
  (with-input-lines (lines)
    (with-output-to-file (s)
      (format s "~{~a~^ ~}~%"
	      (iter (for line in (rest lines))
		    (collect (majority-element (parse-integer-list line))))))))

(define-rosalind-problem :mer
  "merge sorted arrays"
  (with-input-lines (lines)
    (let* ((list-1 (parse-integer-list (second lines)))
	   (list-2 (parse-integer-list (fourth lines))))
      (with-output-to-file (stream)
	(format stream "~{~a~^ ~}~%"
		(iter (until (and (null list-1) (null list-2)))
		      (cond ((null list-1) (collect (pop list-2)))
			    ((null list-2) (collect (pop list-1)))
			    ((< (car list-1) (car list-2)) (collect (pop list-1)))
			    (t (collect (pop list-2))))))))))

(define-rosalind-problem :2sum
    "2sum"
  (with-input-lines (lines)
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
(define-rosalind-problem :ms
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
(define-rosalind-problem :orf
  "open reading frames"
  (with-single-fasta-line (dna-string)
    (let ((proteins (make-hash-table :test #'equal)))
      (iter (for protein in (transcode-open-frames dna-string))
	    (setf (gethash protein proteins) 1))
      (iter (for protein in (transcode-open-frames (reverse-complement dna-string)))
	    (setf (gethash protein proteins) 1))
      (with-output-to-file (stream)
	(format stream "~{~a~%~}" (alexandria:hash-table-keys proteins))))))

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
(define-rosalind-problem :lgis
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
(define-rosalind-problem :cc
  "connected components"
  (let ((graph (make-graph-from-file *input-filename*)))
    (write-single-output-line (count-connected-components graph))))

(define-rosalind-problem :pmch
  "Perfect matchings and RNA secondary structures"
  (with-single-fasta-line (rna)
    (write-single-output-line
      (* (fact (count #\A rna)) (fact (count #\G rna))))))

(defun check-heap-property (vec pred)
  (macrolet ((element (idx) `(elt vec ,idx))
	     (parent-element (idx) `(element (parent-idx ,idx))))
    (labels ((parent-idx (idx) (floor (1- idx) 2))
	     (has-heap-property-p (idx) (funcall pred (parent-element idx) (element idx))))
      (iter (for idx index-of-vector vec from 1)
	    (always (has-heap-property-p idx))))))
(defun make-heap (data pred)
  (macrolet ((element (idx) `(elt data ,idx))
	     (parent-element (idx) `(element (parent-idx ,idx))))
    (labels ((parent-idx (idx) (floor (1- idx) 2))
	     (has-heap-property-p (idx) (funcall pred (parent-element idx) (element idx)))
	     (bubble-single-element (idx)
	       (psetf (element idx) (parent-element idx)
		      (parent-element idx) (element idx)))
	     (bubble-up (idx)
	       (iter (while (not (zerop idx)))
		     (unless (has-heap-property-p idx)
		       (bubble-single-element idx))
		     (setf idx (parent-idx idx)))))
      (iter (for idx index-of-vector data from 1)
	    (bubble-up idx)))))
(define-rosalind-problem :hea
  "building a heap"
  (with-input-lines (lines)
    (let ((data (parse-integer-vector (second lines))))
      (make-heap data #'>=)
      (with-output-to-file (s)
	(print-integer-vector data s)))))

(defun partition (data &key (test #'<=) (key #'identity) (lo 0) (hi (1- (length data))))
  ;; from wikipedia quicksort http://en.wikipedia.org/wiki/Quicksort#Algorithm
  "partitions the data, and returns a list of (pivot-index number-of-required-swaps).
checks if the data is not already sorted."
  (macrolet ((element (idx) `(elt data ,idx)))
    (flet ((sorted-p (idx1 idx2)
	     (funcall test (funcall key (element idx1)) (funcall key (element idx2)))))
      (when (iter (for index from lo to (1- hi))
		  (always (sorted-p index (1+ index))))
	(list lo 0))
      (let* ((pivot-index 0)
	     (pivot-value (element pivot-index))
	     (store-index lo)
	     (swap-count 0))
	(flet ((swap-elements (idx1 idx2) (progn (psetf (element idx1) (element idx2)
							(element idx2) (element idx1))
						 (incf swap-count))))
	  (swap-elements pivot-index hi)
	  (iter (for idx from lo to (1- hi))
		(when (funcall test (element idx) pivot-value)
		  (swap-elements idx store-index)
		  (incf store-index)))
	  (swap-elements store-index hi)
	  (iter (for idx index-of-vector data to (1- store-index))
		(assert (sorted-p idx store-index)))
	  (iter (for idx index-of-vector data from (1+ store-index))
		(assert (not (sorted-p idx store-index))))
	  (list store-index swap-count))))))
(define-rosalind-problem :par
  "2-way partition"
  (with-input-lines (lines)
    (let ((data (parse-integer-vector (second lines))))
      (partition data)
      (with-output-to-file (s)
	(print-integer-vector data s)))))

(define-rosalind-problem :tree
  "completing a tree"
  (let ((graph (make-graph-from-file *input-filename*)))
    (write-single-output-line (1- (count-connected-components graph)))))

(define-rosalind-problem :bip
  "testing bipartiteness"
  (with-input-groups (line-groups)
    (with-output-to-file (s)
      (let ((bipartite-info
	     (iter (for graph-lines in (rest line-groups))
		   (collect (if (is-bipartite-p (make-graph-from-lines graph-lines)) 1 -1)))))
       (print-integer-list bipartite-info s)))))

(defun find-3sum-positions (line)
  (format t "line~%")
  (let* ((elements (iter (for index from 1)
			 (for value in (parse-integer-list line))
			 (collect (cons value index))))
	 (sorted-elements (coerce (sort elements #'<= :key #'car) 'vector))
	 (se-length (length sorted-elements))
	 (max-value (car (elt sorted-elements (1- se-length))))
	 result)
    (iter (for i from 0 to (- se-length 3))
	  (iter (while (and (< i (- se-length 3))
			    (= (car (elt sorted-elements i)) (car (elt sorted-elements (1+ i))))))
		(incf i))
	  (for i-value = (car (elt sorted-elements i)))
	  (until result)
	  (let ((j-start (lower-bound sorted-elements (1+ i) se-length #'< (- 0 i-value max-value) :key #'car)))
	    (iter (for j from j-start to (- se-length 2))
		  (for j-value = (car (elt sorted-elements j)))
		  (let* ((value (- 0 i-value j-value))
			 (position (position-binary-search value sorted-elements #'< :key #'car :begin (1+ j))))
		    (when position
		      (return (setf result (list (cdr (elt sorted-elements i))
						 (cdr (elt sorted-elements j))
						 (cdr (elt sorted-elements position))))))))))
    (if result
	(sort result #'<)
	(list -1))))
(defun find-3sum-positions-fast (line)
  ;; as found at http://www.ti.inf.ethz.ch/ew/courses/CG09/materials/v12.pdf
  (let* ((elements (iter (for index from 1)
			 (for value in (parse-integer-list line))
			 (collect (cons value index))))
	 (sorted-elements (coerce (sort elements #'<= :key #'car) 'vector))
	 (se-length (length sorted-elements))
	 result)
    (flet ((value (i) (car (elt sorted-elements i)))
	   (index (i) (cdr (elt sorted-elements i))))
      (iter (for i from 0 to (- se-length 3))
	    (for i-value = (value i))
	    (until result)
	    (let ((j (1+ i))
		  (k (1- se-length)))
	      (iter (while (> k j))
		    (let ((sum (+ i-value (value k) (value j))))
		      (cond ((zerop sum) (return (setf result (list (index i) (index j) (index k)))))
			    ((> sum 0) (decf k))
			    (t (incf j))))))))
    (if result
	(sort result #'<)
	(list -1))))
(define-rosalind-problem :3sum
  "3sum"
  (with-input-lines (lines)
    (with-output-to-file (s)
      (iter (for line in (rest lines))
	    (print-integer-list (find-3sum-positions-fast line) s)))))

(define-rosalind-problem :dag
  "testing acyclicity"
  (with-input-groups (input-data)
    (with-output-to-file (s)
      (print-integer-list
       (iter (for graph-data in (rest input-data))
	     (if (is-directed-acyclic-graph-p (make-graph-from-lines graph-data :directed))
		 (collect 1)
		 (collect -1)))
       s))))

(define-rosalind-problem :ts
  "topological sort"
  (with-input-lines (lines)
    (with-output-to-file (s)
      (print-integer-list
       (mapcar #'value (topological-sort (make-graph-from-lines lines :directed)))
       s))))

(define-rosalind-problem :hdag
  "hamiltonian path in DAG"
  (with-input-groups (input-data)
    (with-output-to-file (s)
      (iter (for graph-data in (rest input-data))
	    (for graph = (make-graph-from-lines graph-data :directed))
	    (for topo-sort = (topological-sort graph))
	    (for count-components = (count-connected-components graph))
	    (if (= count-components (length (nodes graph)))
		(print-integer-list (mapcar #'value topo-sort) s)
		(format s "-1~%"))))))

(define-rosalind-problem :inod
  "counting phylogenetic ancestors"
  (with-single-input-line (line)
    (write-single-output-line (- (parse-integer line) 2))))

(defun partition-count-swaps (v &key (test #'<) (key #'identity) (from 0) (to (1- (length v))))
  (let ((swap-count 0))
    (flet ((swap-elements (idx1 idx2)
	     (psetf (elt v idx1) (elt v idx2)
		    (elt v idx2) (elt v idx1))
	     (incf swap-count)))
      (let* ((pivot-index to)
	     (pivot (funcall key (elt v pivot-index)))
	     (store-index from))
	(iter (for index from from to (1- to))
	      (when (funcall test (funcall key (elt v index)) pivot)
		(swap-elements index store-index)
		(incf store-index)))
	(swap-elements store-index to)
	(list swap-count store-index)))))
(defun qsort-count-swaps (v &key (test #'<) (key #'identity) (from 0) (to (1- (length v))))
  (flet ((sorted-p (idx1 idx2) (funcall test (funcall key (elt v idx1)) (funcall key (elt v idx2)))))
    (let ((range-length (- to from -1)))
      (cond ((= range-length 1) 0)
	    ((= range-length 2) (if (sorted-p from to) 0 1))
	    (t (destructuring-bind (count pivot-position)
		   (partition-count-swaps v :test test :key key :from from :to to)
		 (+ count
		    (qsort-count-swaps v :test test :key key :from from :to pivot-position)
		    (qsort-count-swaps v :test test :key key :from (1+ pivot-position) :to to))))))))
(define-rosalind-problem :inv
    "counting inversions"
  (with-input-lines (lines)
    (write-single-output-line
      (let ((values (parse-integer-vector (second lines))))
	(iter (for v-i in-vector values with-index i from 0 below (1- (length values)))
	      (iter (for v-j in-vector values with-index j from (1+ i))
		    (in outer (count (> v-i v-j)))))))))
