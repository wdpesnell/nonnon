# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=polyrk - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to polyrk - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "polyrk - Win32 Release" && "$(CFG)" != "polyrk - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "polyrk.mak" CFG="polyrk - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "polyrk - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "polyrk - Win32 Debug" (based on "Win32 (x86) Console Application")
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
# PROP Target_Last_Scanned "polyrk - Win32 Debug"
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "polyrk - Win32 Release"

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

ALL : "$(OUTDIR)\polyrk.exe"

CLEAN : 
	-@erase ".\Release\polyrk.exe"
	-@erase ".\Release\banbks.obj"
	-@erase ".\Release\Polyrk.obj"
	-@erase ".\Release\bandec.obj"
	-@erase ".\Release\Cjhdmp.obj"
	-@erase ".\Release\Pekeris.obj"

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)/polyrk.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/polyrk.pdb" /machine:I386 /out:"$(OUTDIR)/polyrk.exe" 
LINK32_OBJS= \
	".\Release\banbks.obj" \
	".\Release\Polyrk.obj" \
	".\Release\bandec.obj" \
	".\Release\Cjhdmp.obj" \
	".\Release\Pekeris.obj"

"$(OUTDIR)\polyrk.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "polyrk - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "polyrk__"
# PROP BASE Intermediate_Dir "polyrk__"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "debug"
# PROP Intermediate_Dir "debug"
# PROP Target_Dir ""
OUTDIR=.\debug
INTDIR=.\debug

ALL : "$(OUTDIR)\polyrk.exe"

CLEAN : 
	-@erase ".\debug\polyrk.exe"
	-@erase ".\debug\Cjhdmp.obj"
	-@erase ".\debug\bandec.obj"
	-@erase ".\debug\banbks.obj"
	-@erase ".\debug\Polyrk.obj"
	-@erase ".\debug\Pekeris.obj"
	-@erase ".\debug\polyrk.ilk"
	-@erase ".\debug\polyrk.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "polyrk__/" /c /nologo
# ADD F90 /G5 /Zi /I "debug/" /c /nologo
F90_PROJ=/G5 /Zi /I "debug/" /c /nologo /Fo"debug/" /Fd"debug/polyrk.pdb" 
F90_OBJS=.\debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/polyrk.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# SUBTRACT LINK32 /map
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/polyrk.pdb" /debug /machine:I386 /out:"$(OUTDIR)/polyrk.exe" 
LINK32_OBJS= \
	".\debug\Cjhdmp.obj" \
	".\debug\bandec.obj" \
	".\debug\banbks.obj" \
	".\debug\Polyrk.obj" \
	".\debug\Pekeris.obj"

"$(OUTDIR)\polyrk.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

# Name "polyrk - Win32 Release"
# Name "polyrk - Win32 Debug"

!IF  "$(CFG)" == "polyrk - Win32 Release"

!ELSEIF  "$(CFG)" == "polyrk - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Polyrk.for

"$(INTDIR)\Polyrk.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Cjhdmp.for

"$(INTDIR)\Cjhdmp.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Pekeris.for

"$(INTDIR)\Pekeris.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=\Msdev\Numrec\Nrf\bandec.for

"$(INTDIR)\bandec.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\Msdev\Numrec\Nrf\banbks.for

"$(INTDIR)\banbks.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
# End Project
################################################################################
