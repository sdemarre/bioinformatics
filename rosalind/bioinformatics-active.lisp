(in-package :bioinformatics)

(define-rosalind-problem :1a
  "frequent words problem"
  (with-input-lines (lines)
    (let ((dna-string (first lines))
	  (k (parse-integer (second lines))))
      (with-output-to-file (stream)
	(format stream "~{~a~^ ~}~%" (most-frequent-kmers dna-string k))))))

(define-rosalind-problem :1b
  "reverse complement problem"
  (with-single-input-line (dna-string)
    (with-output-to-file (stream)
      (format stream "~a" (reverse-complement dna-string)))))

(define-rosalind-problem :2a
  "protein translation problem"
  (with-input-lines (lines)
    (let ((proteins (do-on-codons (first lines) #'rna-to-protein)))
      (with-output-to-file (stream)
	(format stream "~{~a~}~%" (mapcar #'car (butlast proteins)))))))

