(in-package :bioinformatics)

(define-rosalind-problem :ini3
  "strings and lists"
  (with-input-lines (lines)
    (destructuring-bind-integers (start1 end1 start2 end2) (second lines)
      (with-output-to-file (s)
	  (format s "~a ~a" (subseq (first lines) start1 (1+ end1)) (subseq (first lines) start2 (1+ end2)))))))

(define-rosalind-problem :ini4
  "conditions and loops"
  (with-input-lines (lines)
    (destructuring-bind-integers (a b) (first lines)
      (write-single-output-line
	(iter (for i from a to b)
	      (when (oddp i)
		(summing i)))))))

(define-rosalind-problem :ini5
  "working with files"
  (with-input-lines (lines)
    (with-output-to-file (stream)
      (iter (for line in lines)
	    (for line-number from 1)
	    (when (evenp line-number)
	      (format stream "~a~%" line))))))

(define-rosalind-problem :ini6
  "dictionaries"
  (with-input-lines (lines)
    (let ((word-counts (make-hash-table :test #'equal)))
      (iter (for word in (split-sequence:split-sequence #\Space (first lines)))
	    (incf (gethash word word-counts 0)))
      (with-output-to-file (stream)
       (iter (for (word count) in-hashtable word-counts)
	     (format stream "~a ~a~%" word count))))))




