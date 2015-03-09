(in-package :bioinformatics)

(defparameter *rosalind-id-info* (make-hash-table))

(defmacro define-rosalind-problem (rosalind-id input-filename function-name &body body)
  `(let ((input-filename ,input-filename))
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

(defun make-output-filename (input-filename)
  (format nil "~a_output.txt" (subseq input-filename 0 (- (length input-filename) 4))))

(defmacro with-output-to-file ((stream-name) &body body)
  `(with-open-file (,stream-name (make-output-filename input-filename) :direction :output :if-exists :supersede)
     ,@body
     (make-output-filename input-filename)))


