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
    (with-output-to-file (output)
      (print-integer-list
	      (iter (for element in elements-to-find)
		    (collect (let ((pos (position-binary-search element sorted-array #'<)))
		      	       (if pos (1+ pos) -1))))))))

(defun make-adjacency-list (graph-lines &optional (graph-type :undirected))
  "result is a hastable with nodes as keys and list of reachable nodes as value"
  (let ((adjacency-list (make-hash-table)))
    (iter (for line in (rest graph-lines))
	  (destructuring-bind-integers (head-vertex tail-vertex) line
	    (push tail-vertex (gethash head-vertex adjacency-list))
	    (unless (eq graph-type :directed)
	      (push head-vertex (gethash tail-vertex adjacency-list)))))
    adjacency-list))
(defun make-degree-array (adjacency-list number-vertices)
  (let ((degree-array (make-array number-vertices)))
    (iter (for vertex from 1 to number-vertices)
	  (for vertex-index from 0)
	  (setf (elt degree-array vertex-index) (length (gethash vertex adjacency-list))))
    degree-array))
(define-rosalind-problem :deg degree-array
  "degree array"
  (let* ((lines (read-file-lines input-filename))
	 (adjacency-list (make-adjacency-list lines)))
    (destructuring-bind-integers (number-vertices number-edges) (first lines)
      (declare (ignorable number-edges))
      (print-integer-list
	      (iter (for degree in-vector (make-degree-array adjacency-list number-vertices))
		    (collect degree))))))

(define-rosalind-problem :ddeg double-degree-array
  "double-degree array"
  (let* ((lines (read-file-lines input-filename))
	 (adjacency-list (make-adjacency-list lines)))
    (destructuring-bind-integers (number-vertices number-edges) (first lines)
      (declare (ignorable number-edges))
      (let ((degree-array (make-degree-array adjacency-list number-vertices)))
	(print-integer-list
	 (iter (for vertex from 1 to number-vertices)
	       (collect (iter (for neighbour in (gethash vertex adjacency-list))
			      (summing (elt degree-array (1- neighbour)))))))))))

(defun shortest-distances (number-vertices adjacency-list)
  (let ((reached-vertices (make-hash-table))
	(reachable-vertices (list 1)))
    (iter (while (not (null reachable-vertices)))
	  (for distance from 0)
	  (let (new-vertices)
	    (iter (for vertex in reachable-vertices)
		  (setf (gethash vertex reached-vertices) distance))
	    (iter (for vertex in reachable-vertices)
		  (iter (for neighbour in (gethash vertex adjacency-list))
			(unless (gethash neighbour reached-vertices)
			  (pushnew neighbour new-vertices))))	    
	    (setf reachable-vertices new-vertices)))
   (iter (for vertex from 1 to number-vertices)
	 (collect (alexandria:if-let ((distance (gethash vertex reached-vertices)))
		    distance
		    -1)))))
(define-rosalind-problem :bfs breadth-first-search
  "breadth first search"
  (let* ((lines (read-file-lines input-filename))
	 (adjacency-list (make-adjacency-list lines :directed)))
    (destructuring-bind-integers (number-vertices number-edges) (first lines)
      (declare (ignore number-edges))
      (print-integer-list (shortest-distances number-vertices adjacency-list)))))

(defun make-random-graph (number-vertices number-edges &optional (stream t))
  (format stream "~a ~a~%" number-vertices number-edges)
  (let (edge-list)
    (iter (until (= number-edges (length edge-list)))
	  (let ((v1 (1+ (random number-vertices)))
		(v2 (1+ (random number-vertices))))
	    (iter (while (= v1 v2))
		  (setf v2 (1+ (random number-vertices))))
	    (pushnew (cons  v1 v2) edge-list
		    :test #'(lambda (e1 e2) (and (= (car e1) (car e2)) (= (cdr e1) (cdr e2)))))))
    (iter (for edge in edge-list)
	  (format stream "~a ~a~%" (car edge) (cdr edge)))))
(defun make-bfs-graph (number-vertices number-edges)
  (with-open-file (stream (make-input-filename :bfs) :direction :output :if-exists :supersede)
    (make-random-graph number-vertices number-edges stream)))

(defun make-bfs-dot ()
  (with-open-file (stream "rosalind_bfs.dot" :direction :output :if-exists :supersede)
    (make-dot-output (make-adjacency-list (rosalind-lines :bfs) :directed) stream)))
(defun make-dot-output (adjacency-list &optional (stream stream))
  (format stream "digraph g {~%")
  (iter (for (vertex neighbours) in-hashtable adjacency-list)
	(iter (for neighbour in neighbours)
	      (format stream "  n~a -> n~a;~%" vertex neighbour)))
  (format stream "}~%"))

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


(defun merge-sort (data &optional (start 0) (end (1- (length data))))
  )
(define-rosalind-problem :ms ros-merge-sort
  "merge sort"
  (let ((data (coerce (parse-integer-list (first (read-file-lines input-filename))) 'vector)))
    (with-output-to-file (stream)
      (merge-sort data)
      (format stream "~{~a~^ ~}~%" (coerce data 'list)))))


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

(define-rosalind-problem :lgis ros-longest-increasing-subseq
  "longest increasing subsequence"
  (with-input-lines (lines)
    (let ((permutation (coerce (parse-integer-list (second lines)) 'vector)))
      (iter (for start-pos from 0 to (1- (length permutation)))))))

(defun reachable-nodes-from (node adjacency-list)
  (let ((reachable-nodes (make-instance 'rosalind-set))
	new-nodes-reachable)
    (add-element reachable-nodes node)
    (iter (setf new-nodes-reachable nil)
	  (do-for-set-elements reachable-nodes
	    #'(lambda (node)
		(iter (for neighbour in (gethash node adjacency-list))
		      (unless (has-element-p reachable-nodes neighbour)
			(setf new-nodes-reachable t)
			(add-element reachable-nodes neighbour)))))
	  (until (not new-nodes-reachable)))
    reachable-nodes))
(defun count-connected-components (number-nodes adjacency-list)
  (let ((unprocessed-nodes (make-instance 'rosalind-set))
	(connected-components 0))
    (iter (for node from 1 to number-nodes)
	  (add-element unprocessed-nodes node))
    (iter (while (not (empty-set-p unprocessed-nodes)))
	  (let ((node (get-a-set-element unprocessed-nodes)))
	    (do-for-set-elements (reachable-nodes-from node adjacency-list)
	      #'(lambda (reached-node) (remove-element unprocessed-nodes reached-node)))
	    (incf connected-components)))
    connected-components))
(define-rosalind-problem :cc ros-connected-components
  "connected components"
  (with-input-lines (lines)
    (let ((adjacency-list (make-adjacency-list lines)))
      (destructuring-bind-integers (nodes edges) (first lines)
	(declare (ignorable edges))
	(count-connected-components nodes adjacency-list)))))

(define-rosalind-problem :pmch ros-perfect-matchings
  "Perfect matchings and RNA secondary structures"
  (with-single-fasta-line (rna)
    (* (fact (count #\A rna)) (fact (count #\G rna)))))
