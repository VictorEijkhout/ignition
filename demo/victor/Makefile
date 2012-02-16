info :
	@echo "no user makes here"
	@echo "call from generated directory with:"
	@echo "  make alg N=nnn"
alg :
	( cat ../textop.in ; cat alg${N}.in ; cat ../texbot.in ) \
	    > alg${N}.tex
	pdflatex alg${N}.tex ; mv alg${N}.pdf alg${N}.pdf
