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
      (format output "~{~a~^ ~}~%"
	      (iter (for element in elements-to-find)
		    (collect (let ((pos (position-binary-search element sorted-array #'<)))
		      	       (if pos (1+ pos) -1))))))))

(define-rosalind-problem :deg "rosalind_deg.txt" degree-array
  "degree array"
  (let ((vertex-data (make-hash-table))
	(lines (read-file-lines input-filename)))
    (destructuring-bind-integers (vertex-count edge-count) (first lines)
      (declare (ignorable edge-count))
      (iter (for line in (rest lines))
	    (destructuring-bind-integers (source-vertex target-vertex) line
	      (incf (gethash target-vertex vertex-data 0))
	      (incf (gethash source-vertex vertex-data 0))))
      (format t "~{~a~^ ~}~%"
	      (iter (for vertex from 1 to vertex-count)
		    (collect (gethash vertex vertex-data)))))))
