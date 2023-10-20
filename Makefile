# This makefile is used to compile ASTR code.
# The compiler: intel fortran compiler
#
FCFLAGS= -O3 -fbounds-check
FC=h5pfc
# INCLUDE= -I/usr/local/include

SRCDIR = src
OBJDIR = obj
BINDIR = bin
# CTRDIR = /home/fangjian/opt/cantera-2.5.1

OPTIONS1 = -J $(OBJDIR)
#OPTIONS2 = -J $(OBJDIR) -DCOMB -I$(CTRDIR)/include/cantera
#OMP = -fopenmp

EXE=PP.exe

# LIBS = -L/usr/lib/x86_64-linux-gnu -lz -lm #-L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread
LIBS= -lz -lm #-L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs=  Singleton.F90 CommVarDefine.F90 LinearAlgegra.F90 WriteTec.F90  vtkio.F90   \
      H5ReadWrite.F90 Chem.F90 Interpolation.F90 ControlFileRead.F90 ReadWrite.F90 \
      GradSolver.F90 BasicFunction.F90 FieldView.F90 BoundaryLayer.F90    \
      Statistic.F90 DataProcess.F90 UserDefine.F90 WingDesign.F90          \
      GridGen.F90 FlowAnalyse.F90 PreProcess.F90 fdm.F90 numerics.F90    \
      PostEntrance.F90 Mainpost.F90 
 
OBJS=$(srs:.F90=.o)

%.o:%.F90
	$(FC) $(INCLUDE) $(FCFLAGS) $(OPTIONS1) $(OPTIONS2) $(OMP)  -c -o $(OBJDIR)/$@  $<
default: $(OBJS)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(OMP)
	cp -v  $(TARGET) ~/bin/

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET)

