# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=nonspu - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to nonspu - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "nonspu - Win32 Release" && "$(CFG)" != "nonspu - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "nonspu.mak" CFG="nonspu - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "nonspu - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "nonspu - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
RSC=rc.exe
F90=fl32.exe

!IF  "$(CFG)" == "nonspu - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\nonspu.exe"

CLEAN : 
	-@erase ".\Release\nonspu.exe"
	-@erase ".\Release\Nonspum.obj"
	-@erase ".\Release\Cmplna.obj"
	-@erase ".\Release\Lnanon.obj"
	-@erase ".\Release\Matrixes.obj"
	-@erase ".\Release\Matset.obj"
	-@erase ".\Release\Nonmisc.obj"
	-@erase ".\Release\Cjhdmp.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/nonspu.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/nonspu.pdb" /machine:I386 /out:"$(OUTDIR)/nonspu.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Nonspum.obj" \
	"$(INTDIR)/Cmplna.obj" \
	"$(INTDIR)/Lnanon.obj" \
	"$(INTDIR)/Matrixes.obj" \
	"$(INTDIR)/Matset.obj" \
	"$(INTDIR)/Nonmisc.obj" \
	"$(INTDIR)/Cjhdmp.obj"

"$(OUTDIR)\nonspu.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "nonspu - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\nonspu.exe"

CLEAN : 
	-@erase ".\Debug\nonspu.exe"
	-@erase ".\Debug\Nonspum.obj"
	-@erase ".\Debug\Cmplna.obj"
	-@erase ".\Debug\Lnanon.obj"
	-@erase ".\Debug\Matrixes.obj"
	-@erase ".\Debug\Matset.obj"
	-@erase ".\Debug\Nonmisc.obj"
	-@erase ".\Debug\Cjhdmp.obj"
	-@erase ".\Debug\nonspu.ilk"
	-@erase ".\Debug\nonspu.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/nonspu.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/nonspu.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/nonspu.pdb" /debug /machine:I386 /out:"$(OUTDIR)/nonspu.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Nonspum.obj" \
	"$(INTDIR)/Cmplna.obj" \
	"$(INTDIR)/Lnanon.obj" \
	"$(INTDIR)/Matrixes.obj" \
	"$(INTDIR)/Matset.obj" \
	"$(INTDIR)/Nonmisc.obj" \
	"$(INTDIR)/Cjhdmp.obj"

"$(OUTDIR)\nonspu.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "nonspu - Win32 Release"
# Name "nonspu - Win32 Debug"

!IF  "$(CFG)" == "nonspu - Win32 Release"

!ELSEIF  "$(CFG)" == "nonspu - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Nonspum.for

"$(INTDIR)\Nonspum.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Cmplna.for

"$(INTDIR)\Cmplna.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Lnanon.for

"$(INTDIR)\Lnanon.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Matrixes.for

"$(INTDIR)\Matrixes.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Matset.for

"$(INTDIR)\Matset.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Nonmisc.for

"$(INTDIR)\Nonmisc.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Cjhdmp.for

"$(INTDIR)\Cjhdmp.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
