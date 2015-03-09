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
(define-rosalind-problem :ins "rosalind_ins.txt" insertion-sort-swaps
  "insertion sort"
  (let* ((data (second (read-file-lines input-filename)))
	 (unsorted-data (integer-list data)))
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
(define-rosalind-problem :bins "rosalind_bins.txt" rosalind-binary-search
  "binary search"
  (let* ((lines (read-file-lines input-filename))
	 (sorted-array (coerce (integer-list (third lines)) 'vector))
	 (elements-to-find (integer-list (fourth lines))))
    (with-output-to-file (output)
      (print-integer-list
	      (iter (for element in elements-to-find)
		    (collect (let ((pos (position-binary-search element sorted-array #'<)))
		      	       (if pos (1+ pos) -1))))))))

(defun make-adjacency-list (lines &optional (graph-type :undirected))
  (let ((adjacency-list (make-hash-table)))
    (iter (for line in (rest lines))
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
(define-rosalind-problem :deg "rosalind_deg.txt" degree-array
  "degree array"
  (let* ((lines (read-file-lines input-filename))
	 (adjacency-list (make-adjacency-list lines)))
    (destructuring-bind-integers (number-vertices number-edges) (first lines)
      (declare (ignorable number-edges))
      (print-integer-list
	      (iter (for degree in-vector (make-degree-array adjacency-list number-vertices))
		    (collect degree))))))

(define-rosalind-problem :ddeg "rosalind_ddeg.txt" double-degree-array
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

(define-rosalind-problem :bfs "rosalind_bfs.txt" breadth-first-search
  "breadth first search"
  (let* ((lines (read-file-lines input-filename))
	 (adjacency-list (make-adjacency-list lines :directed))
	 (reachable-vertices (list 1))
	 (reached-vertices (make-hash-table)))
    (iter (while (not (null reachable-vertices)))
	  (for distance from 0)
	  (let (new-vertices)
	    (iter (for vertex in reachable-vertices)
		  (setf (gethash vertex reached-vertices) distance)
		  (iter (for neighbour in (gethash vertex adjacency-list))
			(unless (gethash neighbour reached-vertices)
			  (pushnew neighbour new-vertices))))	    
	    (setf reachable-vertices new-vertices)))
    (destructuring-bind-integers (number-vertices number-edges) (first lines)
      (declare (ignorable number-edges))
      (print-integer-list
       (iter (for vertex from 1 to number-vertices)
	     (collect (alexandria:if-let ((distance (gethash vertex reached-vertices)))
			distance
			-1)))))))

