(in-package :bioinformatics)

(defun rosalind-url (path)
  (format nil "http://rosalind.info/~a" (string-downcase path)))

(defun get-rosalind-html-lines (path)
  (mapcar #'maybe-trim-eol (split-sequence:split-sequence #\Newline (drakma:http-request (rosalind-url path)))))

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

(defun get-cookie-value (cookie-jar cookie-name)
  (let ((cookie (find cookie-name (drakma:cookie-jar-cookies cookie-jar) :test #'string= :key #'drakma:cookie-name)))
    (when cookie
      (drakma:cookie-value cookie))))

(defun login (username password)
  (format t "attempting to login to rosalind...")
  (let* ((cookies (make-instance 'drakma:cookie-jar)))
    (drakma:http-request (rosalind-url "accounts/login/") :cookie-jar cookies)
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
    (drakma:http-request (rosalind-url (format nil "problems/~a" problem-id))
			 :cookie-jar *rosalind-cookies*)))

(defun rosalind-get-problem-dataset (problem-id)
  (with-rosalind-session
    (flexi-streams:octets-to-string
     (drakma:http-request (rosalind-url (format nil "problems/~a/dataset" problem-id))
			  :cookie-jar *rosalind-cookies*))))

(defun rosalind-maybe-get-dataset-from-website ()
  (when *rosalind-use-website*    
    (with-open-file (problem-input-file *input-filename* :direction :output :if-exists :supersede)
      (format problem-input-file "~a" (rosalind-get-problem-dataset *rosalind-id*)))))

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
	  (drakma:parameter-error () "done!"))))))
