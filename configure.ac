#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([ICHyM], [1.0.0], [fabio.disante@gmail.com])
LT_INIT
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.

AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])

AX_PROG_NF_CONFIG([NC_PREFIX=`nf-config --prefix` && NFAUTO="yes"],
                  [NC_PREFIX="${NETCDF}" && NFAUTO="no"])

if eval "test x$NFAUTO = xno"; then
  AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix` && NCAUTO="yes"],
                    [NC_PREFIX="${NETCDF}" && NCAUTO="no"])
fi

AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval" && NCAUTO="no"], [FAKE="no"]
)

PARALLEL="1"
AC_ARG_ENABLE(mpiserial,
  AC_HELP_STRING([--enable-mpiserial],
                 [Use the included MPI replacement library for single processor]),
  [DONOTCHECKMPI="DONOTCHECKMPI"; PARALLEL="0"]
)

AC_ARG_ENABLE(nc4,
  AC_HELP_STRING([--enable-nc4],
                 [Supply this option if you plan on using NetCDF4/HDF5 option.]),
  [NETCDF4_HDF5="1"],[NETCDF4_HDF5="0"]
)
if eval "test x$NETCDF4_HDF5 = x1"; then
  AM_CPPFLAGS="\$(DEFINE)NETCDF4_HDF5 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

##old # Checks for libraries.
##old # FIXME: Replace `main' with a function in `-lnetcdf':
##old AC_CHECK_LIB([netcdf], [-lnetcdf])
##old # FIXME: Replace `main' with a function in `-lnetcdff':
##old AC_CHECK_LIB([netcdff], [-lnetcdff])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

if eval "test x$NFAUTO = xyes"; then
  NC_INCLUDES="`nf-config --fflags`"
  NC_LDFLAGS=""
  NC_LIBS="`nf-config --flibs`"
else
  if eval "test x$NCAUTO = xno"; then
    NC_INCLUDES="-I$NC_PREFIX/include"
    NC_LDFLAGS="-L$NC_PREFIX/lib"
    NC_LIBS="-lnetcdf"
  else
    NCVER=`nc-config --version | cut -d " " -f 2`
    if eval "test x$NCVER = x4.0.1"; then
      NC_INCLUDES="`nc-config --fflags`"
    else
      # new netcdf library has this
      NC_INCLUDES="-I`nc-config --includedir`"
    fi
    NC_LDFLAGS=""
    NC_LIBS="`nc-config --libs`"
  fi
fi

AC_SUBST(NC_PREFIX)
AC_SUBST(NC_INCLUDES)

AC_LANG_PUSH([C])
AC_C_BIGENDIAN
AC_LANG_POP([C])
__AC_FC_NAME_MANGLING

case $ac_cv_fc_mangling in
  "lower case, underscore, no extra underscore")
    name_mangling="FORTRAN_UNDERSCORE_";;
  "lower case, no underscore")
    name_mangling="FORTRAN_SAME";;
  "upper case, no underscore")
    name_mangling="FORTRAN_CAPS_";;
  "lower case, double underscore")
    name_mangling="FORTRAN_GNUF2C";;
  "lower case, no underscore, no extra underscore")
    name_mangling="FORTRAN_SAME";;
  *)
    name_mangling="FORTRAN_MANGLING_UNKNOWN"
    AC_MSG_WARN([UNKNOWN FORTRAN NAME-MANGLING SCHEME]);;
esac

RR_PATH_NETCDF
RR_PATH_NETCDF_F90
# RR_PATH_FORTRANGIS

test "$osname"=NONE && osname=`uname -s`
OSNAME=$osname
AC_SUBST(OSNAME)

test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`

# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

if echo $osname | grep -i aix >/dev/null 2>&1; then
   SYSDEF="AIX"
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   SYSDEF="DARWIN"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   SYSDEF="UNIXSYSTEMV"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
   SYSDEF="IRIX"
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   SYSDEF="IRIX64"
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   SYSDEF="LINUX"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   SYSDEF="OSF1"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   SYSDEF="SUPERUX"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   SYSDEF="SUNOS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   SYSDEF="T3E"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   SYSDEF="UNICOS"
fi
if test -z "$SYSDEF"; then
   AC_MSG_WARN([OPERATING SYSTEM UNKNOWN])
   SYSDEF="UNKNOWNOS"
fi

AM_CONDITIONAL(OSNAME_AIX, [ test x"$osname" = xAIX ] )
AM_CONDITIONAL(OSNAME_LINUX, [ test x"$osname" = xLinux ] )

AM_CONDITIONAL(COMPILER_IBM, [ test x"$FC" = xxlf2003 -o \
                                    x"$FC" = xxlf2003_r -o \
                                    x"$FC" = xxlf2008 -o \
                                    x"$FC" = xxlf2008_r -o \
                                    x"$FC" = xbgxlf2003 -o \
                                    x"$FC" = xbgxlf2003_r -o \
                                    x"$FC" = xbgxlf2008 -o \
                                    x"$FC" = xbgxlf2008_r ] )
AM_CONDITIONAL(COMPILER_GNU, [ test x"$FC" = xgfortran ] )
AM_CONDITIONAL(COMPILER_G95, [ test x"$FC" = xg95 ] )
AM_CONDITIONAL(COMPILER_SUN, [ test x"$FC" = xsunf95 ] )
AM_CONDITIONAL(COMPILER_PGI, [ test x"$FC" = xpgf95  ] )
AM_CONDITIONAL(COMPILER_INTEL, [ test x"$FC" = xifort -o \
                                      x"$FC" = ximpif90 ])
AM_CONDITIONAL(CHYM_PARALLEL, [ test x"$PARALLEL" = "x1"])

# Force Big endian I/O for known compilers
if test -z "$COMPILER_GNU_TRUE"; then
#for debugging
  #FCFLAGS="-fconvert=big-endian -fno-range-check -finline-functions -g -O0 -fcheck=all -Wall -fbacktrace"
  FCFLAGS="-fconvert=big-endian -fno-range-check -finline-functions "
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  FCFLAGS="-convert big_endian -assume byterecl \
           -heap-arrays -finline-functions"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
fi
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="-fendian=big -fmultiple-save"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="-byteswapio"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  FCFLAGS="-WF,-DIBM"
  CPRDEF="XLF"
  AC_SUBST(FCFLAGS)
fi
if test -x "$COMPILER_SUN_TRUE"; then
  CPRDEF="WORKSHOP"
fi

AC_SUBST(CPRDEF)

MCTPATH="${PWD}/external"
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="-O2 $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_GNU_TRUE"; then
  FCFLAGS="-O2 -mtune=native $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  #FCFLAGS="-O0 -fpe0 -g -check all -traceback -warn all -ftz \
  #         -fimf-arch-consistency=true -fltconsistency -mieee-fp $FCFLAGS"
  FCFLAGS="-O2 -fp-model fast=2 -ip $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  FCFLAGS="-qsuppress=1501-245 -O2 -qstrict -qxlf2003=polymorphic \
           -qmaxmem=-1 -qzerosize $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_SUN_TRUE"; then
  FCFLAGS="-O2 -xfilebyteorder=big16:%all \
           -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="-fastsse $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi

CHYM_ROOT=${PWD}
if eval "test x"$PARALLEL" = x1"; then
  ACX_MPI
else
  FORT_SIZE="real4double8"
  MPIFC=$FC
  AC_SUBST(MPIFC)
  AC_SUBST(FORT_SIZE)
  AM_CPPFLAGS="\$(DEFINE)MPI_SERIAL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi
MPEUCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF \$(DEFINE)$name_mangling"
MPEUFCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF"
if test -n "$SRKDEF"; then
  MPEUCFLAGS="$MPEUCFLAGS \$(DEFINE)$SRKDEF"
  MPEUFCFLAGS="$MPEUFCFLAGS \$(DEFINE)$SRKDEF"
fi
AC_SUBST(MPEUCFLAGS)
AC_SUBST(MPEUFCFLAGS)

F90FLAGS="$FCFLAGS"
AC_SUBST(F90FLAGS)

AC_SUBST(CHYM_ROOT)
AC_PREFIX_DEFAULT(${PWD})

if test -z "$COMPILER_SUN_TRUE"; then
  CPPFLAGS="$AM_CPPFLAGS"
  AC_SUBST(CPPFLAGS)
else
  CPPFLAGS="$AM_CPPFLAGS"
  AC_SUBST(CPPFLAGS)
fi
SYSLIBS="$LIBS"
AC_SUBST(SYSLIBS)
##

##wait AC_SUBST(NC_INCLUDES)
##wait AC_SUBST(NC_LIBS)
##wait AC_SUBST(NC_LDFLAGS)

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([preproc/Makefile])
AC_OUTPUT
