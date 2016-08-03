#!/bin/sh
cd tex

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null
/bin/rm -rf *.idx >  /dev/null
/bin/rm -rf *.ilg >  /dev/null
/bin/rm -rf *.ind >  /dev/null

	pdflatex guide_sofi3D
	bibtex guide_sofi3D
	pdflatex guide_sofi3D
	pdflatex guide_sofi3D
	pdflatex guide_sofi3D
	pdflatex guide_sofi3D

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null
/bin/rm -rf *.idx >  /dev/null
/bin/rm -rf *.ilg >  /dev/null
/bin/rm -rf *.ind >  /dev/null
mv guide_sofi3D.pdf ../

cd ..
