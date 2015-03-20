(in-package :bioinformatics)

(defclass graph ()
  ((nodes :reader node :initform (make-array 10 :adjustable t :fill-pointer 0))
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
(defclass node (property-mixin)
  ((value :accessor value :initarg :value)
   (graph :reader graph :initarg :graph)))
(defclass edge (property-mixin)
  ((head :accessor head :initarg :source)
   (tail :accessor tail :initarg :tail)
   (graph :reader graph :initarg :graph)))

(defmethod add-node ((this graph) node-value)
  (let ((node (make-instance 'node :graph this :value node-value)))
    (vector-push-extend node (slot-value this 'nodes))
    node))

(defmethod add-edge ((this graph) (source node) (tail node))
  (vector-push-extend (make-instance 'edge :source source :tail tail :graph this) (slot-value this 'edges))
  (vector-push-extend tail (gethash source (slot-value this 'adjacency-list) (make-array 10 :adjustable t :fill-pointer 0)))
  (unless (eq :directed (graph-type this))
    (vector-push-extend source (gethash tail (slot-value this 'adjacency-list) (make-array 10 :adjustable t :fill-pointer 0)))))

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
