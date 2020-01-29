# Makefile for cteq code
# K. Park 8/9/2016
#
obj= cteq-database.o cteqpdf.o cteqpdf-f77.o main.o
getinfo: $(obj) cteqpdf.h
	gcc -o cteq cteq-database.o cteqpdf.o cteqpdf-f77.o main.o header.h
main.o:main.c
	gcc -c main.c
cteq-database.o:cteq-database.c
	gcc -c cteq-database.c

cteqpdf.o:cteqpdf.c
	gcc -c	cteqpdf.c	

cteqpdf-f77.o:cteqpdf-f77.c
	gcc -c cteqpdf-f77.c	

clean:
	rm cteq $(obj)
	@echo "cteq, $(obj) files are removed."	
