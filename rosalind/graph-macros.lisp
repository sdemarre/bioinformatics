(in-package :bioinformatics)

(defmacro with-temp-nodes-property ((graph property initial-value) &body body)
  "generate a new property symbol, and set it to the initial graph on all nodes of the graph.
after running body, get rid of that property on the nodes."
  `(let ((,property (gensym)))
     (initialize-property ,graph ,property ,initial-value)
     (unwind-protect
	  ,@body
       (clear-property ,graph ,property))))
