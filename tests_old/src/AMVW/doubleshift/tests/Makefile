include ../environment

TOOBJS := $(patsubst $(TODIR)/%.f90,$(TODIR)/%.o,$(wildcard $(TODIR)/*.f90))
TOSRCS := $(wildcard $(TODIR)/*.f90)

OBJS := $(patsubst $(SRCDIR)/%.f90,$(SRCDIR)/%.o,$(wildcard $(SRCDIR)/*.f90))
SRCS := $(wildcard $(SRCDIR)/*.f90)

MPSRCS := $(wildcard $(MPDIR)/*.f)
MPOBJS := $(MPDIR)/mpfun90.o $(MPDIR)/mpmod90.o $(MPDIR)/mpmodm90.o $(MPDIR)/mpmodx90.o

all: spoly 
#race raceunit racejt racesp racebs

spoly: testDAMVW.f95 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o spoly testDAMVW.f95 $(TOOBJS) $(OBJS) $(LIBS)

$(SRCS):

$(TOSRCS):

$(CSRCS):	

clean:
	-rm spoly

