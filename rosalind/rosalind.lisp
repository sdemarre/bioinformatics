(in-package :bioinformatics)

(defparameter *rosalind-id-info* (make-hash-table))

(defmacro define-rosalind-problem (rosalind-id function-name &body body)
  `(let ((input-filename ,(make-input-filename rosalind-id))
	 (output-filename ,(make-output-filename rosalind-id)))
     (setf (gethash ,rosalind-id *rosalind-id-info*) ',function-name)
     (defun ,function-name ()
       ,@body)))

(defun rosalind-run (rosalind-id)
  (let ((problem-info (gethash rosalind-id *rosalind-id-info*)))
    (if problem-info 
	(funcall problem-info)
	(error "You didn't solve this problem yet (~a)" rosalind-id))))

(defun rosalind-find-function (rosalind-id)
  (let ((problem-info (gethash rosalind-id *rosalind-id-info*)))
    (if problem-info
	problem-info
	(error "You didn't solve this problem yet (~a)" rosalind-id))))

(defun list-solved-rosalind-problems ()
  (iter (for (id fun) in-hashtable *rosalind-id-info*)
	(collect (cons id (documentation (rosalind-find-function id) 'function)))))

(defun make-input-filename (problem-id)
  (format nil "rosalind/data/rosalind_~a.txt" (string-downcase (symbol-name problem-id))))
(defun make-output-filename (problem-id)
  (format nil "rosalind/output/rosalind_~a_output.txt" (string-downcase (symbol-name problem-id))))

(defmacro with-output-to-file ((stream-name) &body body)
  `(with-open-file (,stream-name output-filename :direction :output :if-exists :supersede)
     ,@body
     output-filename))

(defmacro with-input-lines ((lines-var-name) &body body)
  `(let ((,lines-var-name (read-file-lines input-filename)))
     ,@body))
(defmacro with-single-input-line ((line-var-name) &body body)
  (let ((gslines (gensym)))
   `(with-input-lines (,gslines)
      (let ((,line-var-name (first ,gslines)))
	,@body))))

(defmacro with-fasta-input-lines ((fasta-lines-var-name) &body body)
  `(let ((,fasta-lines-var-name (read-fasta-lines input-filename)))
     ,@body))
(defmacro with-single-fasta-line ((fasta-line-var-name) &body body)
  (let ((flines (gensym)))
    `(with-fasta-input-lines (,flines)
       (let ((,fasta-line-var-name (second (first ,flines))))
	 ,@body))))
