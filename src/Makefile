#---------------------------------------------------------------------
# Define the C compiler and C flags
#---------------------------------------------------------------------

        CC = cc
        CFLAGS = -c -D_XOPEN_SOURCE=600 -std=c99

#---------------------------------------------------------------------
# Define rm & mv  so as not to return errors
#---------------------------------------------------------------------

        RM =  rm -f
        MV =  mv -f
        MKDIR = mkdir -p

#---------------------------------------------------
# Program and library names
#---------------------------------------------------

        PROGRAM = cMagTest
        LIBNAME = libcMag.a

#---------------------------------------------------------------------
# The source files that make up the application.
# The library doesn't need main
#---------------------------------------------------------------------

        SRCS = \
             maggrid.c\
             magfield.c\
             magfieldutil.c\
             mapcolor.c \
             magfielddraw.c \
             magfieldio.c \
             svg.c \
             testdata.c \
             main.c

        LIBSRCS = \
              maggrid.c\
              magfield.c\
              magfieldutil.c\
              mapcolor.c \
              magfielddraw.c \
              magfieldio.c \
              svg.c \
              testdata.c
#---------------------------------------------------------------------
# The object files (via macro substitution)
#---------------------------------------------------------------------

      OBJS = ${SRCS:.c=.o}
      LIBOBJS = ${LIBSRCS:.c=.o}

#--------------------------------------------------------------------
# required libraries
#--------------------------------------------------------------------

       LIBS = -lm -L../lib -lcMag

#---------------------------------------------------------------------
# The includes dir
#---------------------------------------------------------------------

      INCLUDES = -I../includes

#---------------------------------------------------------------------
# The object files 
#---------------------------------------------------------------------

      OBJS = ${SRCS:.c=.o}

#--------------------------------------------------------
# Default target
#--------------------------------------------------------

     all:
	$(MKDIR) ../lib
	$(MKDIR) ../bin

	$(RM) *.o
	$(RM) ../bin/$(PROGRAM)
	$(RM) ../lib/$(LIBNAME)

	$(CC) $(CFLAGS) $(INCLUDES) $(LIBSRCS)
	ar r $(LIBNAME) $(LIBOBJS)
	$(MV) $(LIBNAME) ../lib

	$(CC) $(CFLAGS) $(INCLUDES) main.c

	$(CC) -o $(PROGRAM) $(OBJS) $(LIBS) 
	$(MV) $(PROGRAM) ../bin



