(in-package :bioinformatics)

(defmacro with-temp-nodes-binary-property ((graph property &optional initial-value get-prop set-prop clear-prop) &body body)
  "generate a new property symbol, and set it to the initial graph on all nodes of the graph.
after running body, get rid of that property on the nodes."
  `(let ((,property (gensym)))
     (initialize-property ,graph ,property ,initial-value)
     (flet ((,(if get-prop get-prop
		  (intern (format nil "HAS-~a" (symbol-name property)))) (node)
	      (get-property node ,property))
	    (,(if set-prop set-prop
		  (intern (format nil "SET-~a" (symbol-name property)))) (node)
	      (set-property node ,property t))
	    (,(if clear-prop clear-prop
		  (intern (format nil "CLEAR-~a" (symbol-name property)))) (node)
	      (set-property node ,property nil)))
       	(declare (ignorable (function ,(intern (format nil "CLEAR-~a" (symbol-name property))))))
      (unwind-protect
	   (progn ,@body)
	(clear-property ,graph ,property)))))

(defmacro with-temp-node-binary-properties ((bindings) &body body)
  (if (null bindings)
      body
      `(with-temp-nodes-binary-property ,(first bindings)
	 (with-temp-node-binary-properties ,(rest bindings)
					    ,body))))

(defmacro with-temp-nodes-property ((graph property initial-value &optional get-prop set-prop) &body body)
  "generate a new property symbol, and set it to the initial graph on all nodes of the graph.
after running body, get rid of that property on the nodes."
  `(let ((,property (gensym)))
     (initialize-property ,graph ,property ,initial-value)
     (flet ((,(if get-prop get-prop
		  (intern (format nil "GET-~a" (symbol-name property)))) (node)
	      (get-property node ,property))
	    (,(if set-prop set-prop
		  (intern (format nil "SET-~a" (symbol-name property)))) (node value)
	      (set-property node ,property value)))
      (unwind-protect
	   (progn ,@body)
	(clear-property ,graph ,property)))))


