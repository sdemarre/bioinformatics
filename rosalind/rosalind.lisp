(in-package :bioinformatics)

(defparameter *rosalind-id-info* (make-hash-table))
(defparameter *rosalind-id* nil)
(defparameter *input-filename* nil)
(defparameter *output-filename* nil)
(defparameter *rosalind-show-output* nil)
(defmacro define-rosalind-problem (rosalind-id docstring &body body)
  (let ((function-name (read-from-string (make-rosalind-function-name rosalind-id))))
    `(progn
       (setf (gethash ,rosalind-id *rosalind-id-info*) ',function-name)
       (defun ,function-name ()
	 ,docstring
	 (let ((*input-filename* ,(make-input-filename rosalind-id))
	       (*output-filename* ,(make-output-filename rosalind-id))
	       (*rosalind-id* ,rosalind-id))
	   (rosalind-maybe-get-dataset-from-website)
	   (prog1
	       ,@body
	     (rosalind-maybe-submit-solution-to-website)
	     (maybe-show-output-file)))))))
(defun maybe-show-output-file ()
  (when *rosalind-show-output*
    (iter (for line in-file *output-filename* using #'read-line)
	  (format t "~a~%" line))))

(defun rosalind-run (rosalind-id)
  (let ((problem-info (gethash rosalind-id *rosalind-id-info*)))
    (if problem-info 
	(funcall problem-info)
	(error "You didn't solve this problem yet (~a)" rosalind-id))))

(defun rosalind-run-from-web (rosalind-id)
  (let ((*rosalind-use-website* t))
    (rosalind-run rosalind-id)))

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
(defun make-rosalind-function-name (problem-id)
  (format nil "problem-~a-fun" (string-downcase (symbol-name problem-id))))

(defmacro with-output-to-file ((stream-name) &body body)
  `(with-open-file (,stream-name *output-filename* :direction :output :if-exists :supersede)
     ,@body
     *output-filename*))

(defmacro write-single-output-line (&body body)
  (let ((result (gensym))
	(s (gensym)))
    `(let ((,result (progn ,@body)))
       (with-output-to-file (,s)
	 (format ,s "~a~%" ,result)))))

(defmacro with-input-lines ((lines-var-name) &body body)
  `(let ((,lines-var-name (read-file-lines *input-filename*)))
     ,@body))
(defmacro with-single-input-line ((line-var-name) &body body)
  (let ((gslines (gensym)))
   `(with-input-lines (,gslines)
      (let ((,line-var-name (first ,gslines)))
	,@body))))
(defmacro with-input-groups ((groups-var-name) &body body)
  (let ((lines (gensym)))
    `(with-input-lines (,lines)
       (let ((,groups-var-name (split-lines-into-groups ,lines)))
	 ,@body))))

(defmacro with-fasta-input-lines ((fasta-lines-var-name) &body body)
  `(let ((,fasta-lines-var-name (read-fasta-lines *input-filename*)))
     ,@body))
(defmacro with-single-fasta-line ((fasta-line-var-name) &body body)
  (let ((flines (gensym)))
    `(with-fasta-input-lines (,flines)
       (let ((,fasta-line-var-name (second (first ,flines))))
	 ,@body))))
(defmacro with-fasta-dna-lines ((dna-strings-var-name) &body body)
  (let ((fasta-data (gensym)))
    `(with-fasta-input-lines (,fasta-data)
       (let ((,dna-strings-var-name (mapcar #'second ,fasta-data)))
	 ,@body))))

(defun problem-lines (problem-id)
  (read-file-lines (make-input-filename problem-id)))
(defun problem-fasta-lines (problem-id)
  (read-fasta-lines (make-input-filename problem-id)))
(defun problem-fasta-dna (problem-id)
  (mapcar #'second (problem-fasta-lines problem-id)))

(defun rosalind-run-all ()
  (let ((problem-ids  (mapcar #'car (list-solved-rosalind-problems))))
    (iter (for id in problem-ids)
	  (format t "running ~a~%" id)
	  (rosalind-run id))))

(defmacro between-forms (between-form &body body)
  `(progn
     ,@(iter (for form on body)
	 (collect (car form))
	 (unless (null (cdr form))
	   (collect between-form)))))


(defun split-lines-into-groups (lines)
  "splits lines into groups of lines, separated by empty lines"
  (let (temp-group
	result)
   (iter (for line in lines)
	 (if (string= "" line)
	     (progn
	       (push (nreverse temp-group) result)
	       (setf temp-group nil))
	     (push line temp-group)))
   (unless (null temp-group)
     (push (nreverse temp-group) result))
   (nreverse result)))
