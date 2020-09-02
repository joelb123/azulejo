# Makefile for linux

VERSION=11.0.667
ARCH=i86osx
BITS=32

BINARY_NAME=usearch$(VERSION)_$(ARCH)$(BITS)

default: usearch

usearch: 
	gunzip $(BINARY_NAME).gz
	mv $(BINARY_NAME) usearch
	chmod 755 usearch
	strip usearch
