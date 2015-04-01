(in-package :bioinformatics)

(defclass graph ()
  ((nodes :reader nodes :initform (make-array 10 :adjustable t :fill-pointer 0))
   (edges :reader edges :initform (make-array 10 :adjustable t :fill-pointer 0))
   (type :reader graph-type :initarg :type :initform :undirected)
   (adjacency-list :initform (make-hash-table))))

(defclass property-mixin ()
  ((properties :initform (make-hash-table))))
(defmethod set-property ((this property-mixin) symbol value)
  (setf (gethash symbol (slot-value this 'properties)) value))
(defmethod get-property ((this property-mixin) symbol)
  (gethash symbol (slot-value this 'properties)))
(defmethod has-property ((this property-mixin) symbol)
  (second (multiple-value-list (gethash (slot-value this 'properties) symbol))))
(defmethod list-properties ((this property-mixin))
  (alexandria:hash-table-keys (slot-value this 'properties)))
(defmethod remove-property ((this property-mixin) symbol)
  (remhash symbol (slot-value this 'properties)))
(defclass node (property-mixin)
  ((value :accessor value :initarg :value)
   (graph :reader graph :initarg :graph)))
(defclass edge (property-mixin)
  ((head :accessor head :initarg :source)
   (target :accessor target :initarg :target)
   (graph :reader graph :initarg :graph)))

(defmethod add-node ((this graph) node-value)
  (let ((node (make-instance 'node :graph this :value node-value)))
    (vector-push-extend node (slot-value this 'nodes))
    node))

(defmethod add-edge ((this graph) (source node) (target node))
  (vector-push-extend (make-instance 'edge :source source :target target :graph this) (slot-value this 'edges))
  (macrolet ((get-adjacency-list (s)
	       `(gethash ,s (slot-value this 'adjacency-list))))
    (symbol-macrolet ((source-adjacency-list (get-adjacency-list source))
		      (target-adjacency-list (get-adjacency-list target)))
      (unless source-adjacency-list
	(setf source-adjacency-list (make-array 10 :adjustable t :fill-pointer 0)))
      (vector-push-extend target source-adjacency-list)
      (when (eq :undirected (graph-type this))
	(unless target-adjacency-list
	  (setf target-adjacency-list (make-array 10 :adjustable t :fill-pointer 0)))
	(vector-push-extend source target-adjacency-list)))))

(defmethod count-nodes ((this graph))
  (length (slot-value this 'nodes)))

(defmethod count-edges ((this graph))
  (hash-table-count (slot-value this 'edges)))

(defmethod adjacent-nodes ((node node))
  (iter (for node in-vector (gethash node (slot-value (graph node) 'adjacency-list)))
	(collect node)))

(defmethod degree ((this node))
  (length (gethash this (slot-value (graph this) 'adjacency-list))))

(iter:defmacro-clause (for var node-of-graph g)
  "All nodes of graph"
  (let ((graph (gensym)))
    `(progn
       (with ,graph = ,g)
       (for ,var in-vector (slot-value ,graph 'nodes)))))
(iter:defmacro-clause (for var edge-of-graph g)
  "all edges of graph"
  (let ((graph (gensym)))
    `(progn
       (with ,graph = ,g)
       (for ,var in-vector (slot-value ,graph 'edges)))))
(iter:defmacro-clause (for var neighbour-of-node n)
  "neighbours of node"
  (let ((node (gensym)))
    `(progn
       (with ,node = ,n)
       (for ,var in-vector (gethash ,node (slot-value (graph ,node) 'adjacency-list))))))

(defmethod find-node ((this graph) node-value &optional (test #'eql))
  (iter (for node node-of-graph this)
	(when (funcall test node-value (value node))
	  (return node))))

(defmethod value-to-node-hash ((this graph) &optional (test 'eql))
  (let ((result (make-hash-table :test test)))
    (iter (for node node-of-graph this)
	  (setf (gethash (value node) result) node))
    result))


(defun initialize-property (graph property value)
  (iter (for node node-of-graph graph)
	(set-property node property value)))
(defun clear-property (graph property)
  (iter (for node node-of-graph graph)
	(remove-property node property)))
(defun traverse-nodes-depth-first (initial-node fun &optional loop-fun)
  "traverses nodes, depth first, with loop checking. 
calls (fun node parent-node) for every node visited.
optionally calls (loop-fun node parent-node) when a loop was detected."
  (with-temp-nodes-property ((graph initial-node) visited nil)
    (labels ((node-visited-p (node) (get-property node visited))
	     (visit-node (node) (set-property node visited t))
	     (traverse (node &optional parent)
	       (if (node-visited-p node)
		   (when loop-fun
		     (funcall loop-fun node parent))
		   (progn
		     (funcall fun node parent)
		     (visit-node node)
		     (iter (for neighbour neighbour-of-node node)
			   (unless (eq neighbour parent)
			    (traverse neighbour node)))))))
      (traverse initial-node))))

(defun is-bipartite-p (graph)
  (with-temp-nodes-property (graph color nil)
    (let (conflicting-node-found)
      (flet ((set-color (node color-name) (set-property node color color-name))
	     (get-color (node) (get-property node color))
	     (other-color (color) (case color
				    (:red :blue)
				    (:blue :red))))
	(iter (for node node-of-graph graph)
	      (if conflicting-node-found (return))
	      (unless (get-color node)
		(set-color node :red)
		(traverse-nodes-depth-first node
					    #'(lambda (node parent) (when parent
								      (set-color node (other-color (get-color parent)))))
					    #'(lambda (node parent) (when (eq (get-color node) (get-color parent))
								      (setf conflicting-node-found t)))))))
      (not conflicting-node-found))))
