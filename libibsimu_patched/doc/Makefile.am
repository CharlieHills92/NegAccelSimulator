.PHONY: doc webdoc

doc:
	rm -Rf html
	rm -Rf latex
	doxygen

webdoc:
	rm -Rf webdoc
	doxygen Doxyfile_web

CLEANFILES = *~
