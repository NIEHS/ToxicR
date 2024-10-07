// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
//#pragma message ("Entered stdafx.h")

#include <stdio.h>
#include <limits>
#define NEAR_ZERO 1e-8
//#pragma message ("defining _USE_MATH_DEFINES in stdafx.h")
#define _USE_MATH_DEFINES

#ifndef R_COMPILATION
#pragma message("Not R_COMPILATION in stdafx.h")

#if defined WIN32 || defined _WINDOWS
#include "targetver.h"
#include <windows.h>
#include <tchar.h>
#include <math.h>
#else // i.e., not WIN32 or _WINDOWS
#include <cfloat>
#endif // defined WIN32 || defined _WINDOWS
#endif // ifndef R_COMPILATION

// Custom DEBUG and LOG macros
// - Do not define DEBUGLOG when building R because it will break the build
#ifdef DEBUGLOG
#include <fstream>
#define DEBUG_OPEN_LOG(fname, fvar)                               \
      ofstream(fvar);                                             \
      (fvar).open(fname, fstream::app);                           \
      (fvar) << __FUNCTION__ << " at line: " << __LINE__ << endl; \
      /*(fvar) << std::scientific << setprecision(17);*/
#define DEBUG_LOG(f, s) \
      (f) << s << endl; \
      flush((f));
#define DEBUG_CLOSE_LOG(fvar) (fvar).close();
#else
#define DEBUG_OPEN_LOG(fname, fvar)
#define DEBUG_LOG(f, s)
#define DEBUG_CLOSE_LOG(fvar)
#endif // DEBUGLOG
