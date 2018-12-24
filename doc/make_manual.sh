#!/bin/sh
clean() {
    # Clean tex directory by removing intermediate LaTeX artifacts.
    rm -rf ./*.dvi > /dev/null
    rm -rf ./*.log > /dev/null
    rm -rf ./*.out > /dev/null
    rm -rf ./*.aux > /dev/null
    rm -rf ./*.toc > /dev/null
    rm -rf ./*.blg > /dev/null
    rm -rf ./*.bbl > /dev/null
    rm -rf ./*.lof > /dev/null
    rm -rf ./*.lot > /dev/null
    rm -rf ./*.plt > /dev/null
    rm -rf ./*.fff > /dev/null
    rm -rf ./*.ttt > /dev/null
    rm -rf ./*.tit > /dev/null
    rm -rf ./*.spl > /dev/null
    rm -rf ./*.idx > /dev/null
    rm -rf ./*.ilg > /dev/null
    rm -rf ./*.ind > /dev/null
}


name=asofi3d_user_guide
fullname=${name}.pdf
cd tex

./parseREADME2tex.sh

clean

pdflatex $name
bibtex   $name
pdflatex $name
pdflatex $name
pdflatex $name
pdflatex $name

clean

cd ..

mv tex/$fullname .
