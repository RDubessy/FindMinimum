all:
	cd src && make all
install:
	cp bin/findMinimum /usr/local/bin/
uninstall:
	rm /usr/local/bin/findMinimum
clean:
	cd src && make clean
	cd bin && rm -rf *
	cd doc && rm -rf html
snapshot:
	tar cvjf ~/public_html/src/findMinimum.tar.bz2 ../findMinimum --exclude-vcs --exclude=html --exclude=*.o --exclude=bin/findMinimum --exclude=*.bz2 --exclude=test
	cp -r doc/html/* ~/public_html/doc/findMinimum/
