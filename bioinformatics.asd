;;;; bioinformatics.asd

(asdf:defsystem #:bioinformatics
  :description "Describe bioinformatics here"
  :author "Your Name <your.name@example.com>"
  :license "Specify license here"
  :serial t
  :depends-on (:iterate :alexandria)
  :components ((:file "package")
	       (:file "utils")
               (:file "bioinformatics")
	       (:file "week-1")))


