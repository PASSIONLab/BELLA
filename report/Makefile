#
#	Makefile to compile the paper
#	"make all" or "make"	create the itec99.ps file.
#	"make clean" removes all the auxiliary files.

NAME = report
ACRONYM = TCAD
LAST=1000

report:
	@rm -f *.aux *.bbl *.blg 
	@pdflatex	$(NAME)
	@if test -f $(NAME).aux &&  test `grep citation $(NAME).aux | wc -l` -ge 1; then bibtex $(NAME); fi
	@pdflatex	$(NAME)
	@pdflatex	$(NAME)
	# @dvips -P cmz -t letter -o $(NAME).ps $(NAME).dvi -l$(LAST)
	# @ps2pdf $(NAME).ps

clean:
	rm -f *.fdb_latexmk *.fls*.out *.aux *.bbl *.blg *.dvi *.log *~ *.backup $(NAME).ps $(NAME).pdf $(NAME).synctex.gz *.aux

backup:
	tar cvf $(DATA)-$(ACRONYM).tar *
	gzip $(DATA)-$(ACRONYM).tar
