(in-package :bioinformatics)
;; from clocc-cllib/simple

(defmacro with-collect ((&rest collectors) &body forms)
  "Evaluate forms, collecting objects into lists.
Within the FORMS, you can use local macros listed among collectors,
they are returned as multiple values.
E.g., (with-collect (c1 c2) (dotimes (i 10) (if (oddp i) (c1 i) (c2 i))))
 ==> (1 3 5 7 9); (0 2 4 6 8) [2 values]
In CLISP, push/nreverse is about 1.25 times as fast as pushing into the
tail, so this macro uses push/nreverse on CLISP and push into the tail
on other lisps (which is 1.5-2 times as fast as push/nreverse there)."
  #+clisp
  (let ((ret (mapcar (lambda (cc) (gensym (format nil "~s-RET-" cc)))
                     collectors)))
    `(let (,@ret)
      (declare (list ,@ret))
      (macrolet ,(mapcar (lambda (co re) `(,co (form) `(push ,form ,',re)))
                         collectors ret)
        ,@forms
        (values ,@(mapcar (lambda (re) `(sys::list-nreverse ,re)) ret)))))
  #-clisp
  (let ((ret (mapcar (lambda (cc) (gensym (format nil "~s-RET-" cc)))
                     collectors))
        (tail (mapcar (lambda (cc) (gensym (format nil "~s-TAIL-" cc)))
                      collectors))
        (tmp (mapcar (lambda (cc) (gensym (format nil "~s-TMP-" cc)))
                     collectors)))
    `(let (,@ret ,@tail)
      (declare (list ,@ret ,@tail))
      (macrolet ,(mapcar (lambda (co re ta tm)
                           `(,co (form)
                             `(let ((,',tm (list ,form)))
                               (if ,',re (setf (cdr ,',ta) (setf ,',ta ,',tm))
                                   (setf ,',re (setf ,',ta ,',tm))))))
                         collectors ret tail tmp)
        ,@forms
        (values ,@ret)))))
