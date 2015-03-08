(in-package :rosalind)

(define-rosalind-problem :ini3 "rosalind_ini3.txt" multiple-substrings
  "strings and lists"
  (let ((lines (read-file-lines input-filename)))
    (destructuring-bind-integers (start1 end1 start2 end2) (second lines)
      (format nil "~a ~a" (subseq (first lines) start1 (1+ end1)) (subseq (first lines) start2 (1+ end2))))))

(define-rosalind-problem :ini4 "rosalind_ini4.txt" sum-odd-ints
  "conditions and loops"
  (let ((lines (read-file-lines input-filename)))
    (destructuring-bind-integers (a b) (first lines)
      (iter (for i from a to b)
	    (when (oddp i)
	      (summing i))))))

(define-rosalind-problem :ini5 "rosalind_ini5.txt" rosalind-even-lines-from-file
  "working with files"
  (let ((lines (read-file-lines input-filename)))
    (iter (for line in lines)
	  (for line-number from 1)
	  (when (evenp line-number)
	    (format t "~a~%" line)))))

(define-rosalind-problem :ini6  "rosalind_ini6.txt" count-words
  "dictionaries"
  (let ((lines (read-file-lines input-filename))
	(word-counts (make-hash-table :test #'equal)))
    (iter (for word in (split-sequence:split-sequence #\Space (first lines)))
	  (incf (gethash word word-counts 0)))
    (iter (for (word count) in-hashtable word-counts)
	  (format t "~a ~a~%" word count))))




