;;;; bioinformatics.asd

(asdf:defsystem #:bioinformatics
  :description "Describe bioinformatics here"
  :author "Your Name <your.name@example.com>"
  :license "Specify license here"
  :serial t
  :depends-on (:iterate :alexandria :cl-ppcre :split-sequence)
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


