(in-package :bioinformatics)

(defparameter *rosalind-username* nil)
(defparameter *rosalind-password* nil)
(defparameter *rosalind-cookies* nil)
(defparameter *rosalind-use-website* nil)

(defmacro with-rosalind-session (&body body)
  `(progn
     (maybe-init-rosalind-session)
     ,@body))
