;;;; bioinformatics.asd

(asdf:defsystem #:bioinformatics
  :description "stuff developed while learning bioinformatics"
  :author "Serge De Marre <sdemarre@gmail.com>"
  :license "Specify license here"
  :serial t
  :depends-on (:iterate :alexandria :cl-ppcre :split-sequence :drakma :cl-fad :parse-number)
  :components ((:file "package")
	       (:file "utils")
               (:file "bioinformatics")
	       (:file "rosalind/rosalind")
	       (:file "rosalind/set")
	       (:file "rosalind/allele-utils")
	       (:file "rosalind/python-village")
	       (:file "rosalind/bioinformatics-stronghold")
	       (:file "rosalind/algorithmic-heights")
	       (:file "week-1")))


