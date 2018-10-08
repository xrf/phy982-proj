.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o .pdf .png .svg

CFLAGS=-Wall -fPIC -g -O2 -std=c99
LIBS=-lm -llapack -lslatec

all: index.html paper.pdf

run: proj.py libproj.so
	./proj.py

clean:
	rm -fr dist index.html paper.pdf proj.py.cache *.o

# probably only works on GNU sed
update-makefile:
	$(CC) >Makefile.insert.tmp -MM $(objs:.o=.c)
	cp Makefile Makefile.bak
	a='^# AUTO-GENERATED BEGIN$$' && \
	z='^# AUTO-GENERATED END$$' && \
	p="/$$a/,/$$z/{/$$a/{p;r Makefile.insert.tmp\n};/$$z/p;d}" && \
	p=`printf "$$p"` && \
	sed -i "$$p" Makefile
	rm Makefile.insert.tmp

.PHONY: all clean run update-makefile

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

.svg.png:
	inkscape -D -d 300 -e $@ $<

.svg.pdf:
	inkscape -D -z -A $@ $<

paper.pdf: \
    paper.tex \
    paper.bib \
    proj-atan.pdf \
    proj-contour.pdf \
    proj-completeness.pdf \
    proj-convergence-count.pdf \
    proj-convergence-height.pdf \
    proj-l0-bound.pdf \
    proj-l0.pdf \
    proj-l2-resonance.pdf \
    proj-l2-resonance-asymptotic.pdf \
    proj-l2.pdf \
    proj-poles.pdf \
    texc
	./texc --silent $(@:.pdf=) proj-*.pdf

mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js
index.html: \
    index.md \
    template.html \
    proj-atan.png \
    proj-bessel.png \
    proj-contour.png \
    proj-completeness.png \
    proj-convergence-count.png \
    proj-convergence-height.png \
    proj-l0-bound.png \
    proj-l0.png \
    proj-l2-resonance.png \
    proj-l2.png \
    proj-poles.png \
    reveal.js/js/reveal.js
	pandoc -s -t revealjs --section-divs --no-highlight \
	--mathjax=$(mathjax)?config=TeX-AMS-MML_HTMLorMML \
	--template template.html -V theme:white -V transition:fade \
	-V revealjs-url:reveal.js/ -o $@ index.md

reveal.js/js/reveal.js: reveal.js.tar.gz
	tar xmzf reveal.js.tar.gz

objs=proj.o bessel.o gauss_kronrod.o xmalloc.o
libproj.so: $(objs)
	$(CC) $(LDFLAGS) -shared -o $@ $(objs) $(LIBS)

# AUTO-GENERATED BEGIN
proj.o: proj.c bessel.h math_defs.h gauss_kronrod.h
bessel.o: bessel.c bessel.h math_defs.h
gauss_kronrod.o: gauss_kronrod.c xmalloc.h gauss_kronrod.h math_defs.h
xmalloc.o: xmalloc.c
