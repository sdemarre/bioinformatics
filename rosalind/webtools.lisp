(in-package :bioinformatics)

(defun rosalind-url (path)
  (format nil "http://rosalind.info/~a" (string-downcase path)))

(defun rosalind-session-expired-p ()
  (or (< (drakma:cookie-expires (first (drakma:cookie-jar-cookies *rosalind-cookies*))) (get-universal-time))
      (< (drakma:cookie-expires (second (drakma:cookie-jar-cookies *rosalind-cookies*))) (get-universal-time))))

(defun maybe-init-rosalind-session ()
  (unless (and *rosalind-cookies* (not (rosalind-session-expired-p)))
    (unless (and *rosalind-username* *rosalind-password*)
      (error "You should set *rosalind-username* and *rosalind-password*"))
    (login *rosalind-username* *rosalind-password*)))

(defun logout ()
  (if *rosalind-cookies*
      (drakma:http-request (rosalind-url "accounts/logout") :cookie-jar *rosalind-cookies*)
      (error "you are not logged in, use (login)"))
  (setf *rosalind-cookies* nil))

(defun rosalind-http-request (path)
  (drakma:http-request (rosalind-url path) :cookie-jar *rosalind-cookies*))

(defun get-cookie-value (cookie-jar cookie-name)
  (let ((cookie (find cookie-name (drakma:cookie-jar-cookies cookie-jar) :test #'string= :key #'drakma:cookie-name)))
    (when cookie
      (drakma:cookie-value cookie))))

(defun login (username password)
  (format t "attempting to login to rosalind...")
  (let* ((cookies (make-instance 'drakma:cookie-jar)))
    (rosalind-http-request "accounts/login/")
    (let ((login-result
	   (multiple-value-list
	    (drakma:http-request (rosalind-url "accounts/login/")
				 :method :post
				 :cookie-jar cookies
				 :parameters `(("csrfmiddlewaretoken" . ,(get-cookie-value cookies "csrftoken"))
					       ("username" . ,username)
					       ("password" . ,password))))))
      (if (= 200 (second login-result))
	  (setf *rosalind-cookies* cookies)
	  (error "failed to log in. did you set *rosalind-username* and *rosalind-password* ?")))))

(defun rosalind-get-problem-page (problem-id)
  (with-rosalind-session
    (rosalind-http-request (format nil "problems/~a" problem-id))))

(defun rosalind-get-problem-dataset (problem-id)
  (with-rosalind-session
    (flexi-streams:octets-to-string
     (rosalind-http-request (format nil "problems/~a/dataset" problem-id)))))

(defun rosalind-maybe-get-dataset-from-website ()
  (when *rosalind-use-website*
    (with-rosalind-session
      (with-open-file (problem-input-file *input-filename* :direction :output :if-exists :supersede)
	(format problem-input-file "~a" (rosalind-get-problem-dataset *rosalind-id*))))))

(defun rosalind-resubmit (problem-id)
  (let* ((*rosalind-id* problem-id)
	 (*output-filename* (make-output-filename problem-id))
	 (*rosalind-use-website* t))
    (rosalind-maybe-submit-solution-to-website)))

(defparameter *rosalind-submit-result* nil)
(defun rosalind-maybe-submit-solution-to-website ()
  (when *rosalind-use-website*
    (with-rosalind-session
      (let ((token (get-cookie-value *rosalind-cookies* "csrftoken"))
	    (output-pathname (read-from-string (format nil "#P\"~a\"" *output-filename*))))
	(handler-case
	    (drakma:http-request (rosalind-url (format nil "problems/~a" *rosalind-id*))
				 :method :post
				 :cookie-jar *rosalind-cookies*
				 :parameters `(("csrfmiddlewaretoken" . ,token)
					       ("output_file" . ,output-pathname)))
	  (drakma:parameter-error ()
	    ;; this is a bit strange, maybe a workaround for a bug in drakma ?
	    ;; the last redirect is a :get, and drakma doesn't seem to handle a redirect which becomes a :get after initially being a :post,
	    ;;  it still tries to send the parameters, which it shouldn't (according to the http spec)
	    (setf *rosalind-submit-result*
		  (drakma:http-request (rosalind-url (format nil "problems/~a" *rosalind-id*))
				       :method :get
				       :cookie-jar *rosalind-cookies*))
	    (cond ((cl-ppcre:scan "Sorry, wrong answer" *rosalind-submit-result*)
		   (format t "Sorry, wrong answer."))
		  ((cl-ppcre:scan "Congratulations" *rosalind-submit-result*)
		   (format t "Congratulations.")))))))))

(defun phtml-node-is-text-p (node)
  (stringp node))
(defun phtml-node-children (node)
  (if (consp node)
      (rest node)))
(defun phtml-node-tag (node)
  (when (consp node)
    (if (consp (car node))
	(caar node)
	(car node))))
(defun phtml-node-attributes (node)
  (when (consp node)
    (if (consp (car node))
	(rest (car node))
	(unless (stringp (second node))
	  (rest node)))))
(defun phtml-node-get-attribute (node attribute)
  (getf (phtml-node-attributes node) attribute))

(defun find-nodes-with-tag (node tag)
  (unless (phtml-node-is-text-p node)
    (if (eq (phtml-node-tag node) tag)
	(list node)
	(iter (for child in (phtml-node-children node))
	      (appending (find-nodes-with-tag child tag))))))
(defun find-nodes (node tag attribute value)
  (let ((nodes (find-nodes-with-tag node tag)))
    (remove-if-not #'(lambda (n) (string= (phtml-node-get-attribute n attribute) value)) nodes)))

(defun list-rosalind-problems ()
  (flet ((get-problem-name (aref-node)
	   (third (split-sequence:split-sequence #\/ (phtml-node-get-attribute aref-node :href)))))
    (with-rosalind-session
      (let ((problem-html (cl-html-parse:parse-html (rosalind-http-request "problems/list-view"))))
	(flet ((problems-of-type (type)
		 (mapcar #'get-problem-name (find-nodes (second problem-html) :a :class type))))
	  (list (list :solved (problems-of-type "solved"))
		(list :accessible (problems-of-type "accessible"))
		(list :not-accessible (problems-of-type "not-accessible"))))))))

