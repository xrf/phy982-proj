.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o .png .svg

CFLAGS=-Wall -fPIC -g -O2 -std=c99
LIBS=-lm -lslatec

all: index.html

run: proj.py libproj.so
	./proj.py

clean:
	rm -fr dist *.o

# this sed script probably only works on GNU sed
update-makefile:
	a='^# AUTO-GENERATED BEGIN$$' && \
	z='^# AUTO-GENERATED END$$' && \
	p="/$$a/,/$$z/{/$$a/{p;r Makefile.insert.tmp\n};/$$z/p;d}" && \
	p=`printf "$$p"` && \
	gcc >Makefile.insert.tmp -MM $(objs:.o=.c) && \
	cp Makefile Makefile.bak && \
	sed -i "$$p" Makefile && \
	rm Makefile.insert.tmp

.PHONY: all clean run update-makefile

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

.svg.png:
	inkscape -D -d 300 -e $@ $<

mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js
index.html: \
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
    reveal.js/css/reveal.css
	pandoc -s -t revealjs --section-divs --no-highlight \
	--mathjax=$(mathjax)?config=TeX-AMS-MML_HTMLorMML \
	--template template.html -V theme:white -V transition:fade \
	-V revealjs-url:reveal.js/ -o $@ index.md

reveal.js/css/reveal.css: reveal.js.tar.gz
	tar xmzf reveal.js.tar.gz

objs=proj.o bessel.o gauss_kronrod.o xmalloc.o
libproj.so: $(objs)
	$(CC) $(LDFLAGS) -shared -o $@ $(objs) $(LIBS)

# AUTO-GENERATED BEGIN
proj.o: proj.c bessel.h math_defs.h gauss_kronrod.h
bessel.o: bessel.c bessel.h math_defs.h
gauss_kronrod.o: gauss_kronrod.c xmalloc.h gauss_kronrod.h math_defs.h
xmalloc.o: xmalloc.c
