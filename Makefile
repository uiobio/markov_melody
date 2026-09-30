# Makefile to build class '_template_' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = markov_melody

# input source file (class name == source file basename)
class.sources = markov.c

# all extra files to be included in binary distribution of the library

datafiles = markov-help.pd markov-meta.pd

# include Makefile.pdlibbuilder

  include ./Makefile.pdlibbuilder

# simplistic tests whether all expected files have been produced/installed
buildcheck: all
	test -e markov.$(extension)
installcheck: install
	test -e $(installpath)/markov.$(extension)
