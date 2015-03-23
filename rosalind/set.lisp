(in-package :bioinformatics)

(defclass rosalind-set ()
  ((max-element :accessor max-element :initarg :max-element)
   (elements :reader elements :initarg :elements :initform (make-hash-table))))

(defmethod add-element ((this rosalind-set) el)
  (setf (gethash el (elements this)) 1))

(defmethod empty-set-p ((this rosalind-set))
  (zerop (hash-table-count (elements this))))

(defmethod remove-element ((this rosalind-set) el)
  (remhash el (elements this)))

(defun do-for-set-elements (set fun)
  (iter (for (key value) in-hashtable (elements set))
	(funcall fun key)))

(defun has-element-p (set el)
  (gethash el (elements set)))

(defun set-size (set)
  (hash-table-count (elements set)))

(defmacro with-new-set (set-name &body body)
  `(let ((,set-name (make-instance 'rosalind-set)))
     ,@body
     ,set-name))

(defmethod rset-union ((this rosalind-set) other-set)
  "rosalind-set union"
  (with-new-set union
    (do-for-set-elements this #'(lambda (el) (add-element union el)))
    (do-for-set-elements other-set #'(lambda (el) (add-element union el)))))

(defmethod rset-intersection ((this rosalind-set) other-set)
  "rosalind-set intersection"
  (with-new-set intersection
    (do-for-set-elements this #'(lambda (el) (when (has-element-p other-set el)
					       (add-element intersection el))))))

(defmethod rset-difference ((this rosalind-set) other-set)
  "rosalind-set difference"
  (with-new-set difference
      (do-for-set-elements this #'(lambda (el) (unless (has-element-p other-set el)
						 (add-element difference el))))))

(defmethod rset-equal ((this rosalind-set) (other rosalind-set))
  (and (empty-set-p (rset-difference this other))
       (empty-set-p (rset-difference other this))))

(defmethod rset-complement ((this rosalind-set))
  "rosalind-set complement"
  (with-new-set complement
    (iter (for el from 1 to (max-element this))
	  (unless (has-element-p this el)
	    (add-element complement el)))))

(defun parse-rosalind-set (max-element line)
  "{1, 2, 3} -> rosalind-set with elements 1 2 and 3"
  (with-new-set result
    (cl-ppcre:do-matches (start end "[0-9]+" line)
      (add-element result (parse-integer (subseq line start end))))
    (setf (max-element result) max-element)))

(defun print-set (rosalind-set &optional (stream t))
  (let (elements)
    (do-for-set-elements rosalind-set #'(lambda (el) (push el elements)))
    (format stream "{~{~a~^, ~}}" (nreverse elements))))

(defmethod print-object ((this rosalind-set) stream)
  (print-unreadable-object (this stream :type 'rosalind-set)
    (print-set this stream)))

(defmethod get-a-set-element ((this rosalind-set))
  (iter (for (key value) in-hashtable (elements this))
	(declare (ignorable value))
	(return key)))
