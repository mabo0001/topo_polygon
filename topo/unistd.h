/****************************************************************************************
 * unistd.h                                                                             *
 * 2008-03-15 Last created by cheungmine.                                               *
 *  All rights reserved by cheungmine.                                                  *
 ****************************************************************************************/
#ifndef UNISTD_H__
#define UNISTD_H__

/* Standard C header files included */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*============================================================================*/

#ifndef HAVE_INT8
	#define HAVE_INT8
	typedef	signed char int8;	/* NB: non-ANSI compilers may not grok */
	typedef	unsigned char uint8, uchar, byte, BYTE;
#endif

#ifndef HAVE_INT16
	#define HAVE_INT16
	typedef	short int16;
	typedef unsigned short uint16, word_t, ushort;	/* sizeof (uint16) must == 2 */
#endif

#ifndef HAVE_INT32
	#define HAVE_INT32
	typedef	int int32;
	typedef	unsigned int uint32, size_t, dword_t;	/* sizeof (uint32) must == 4 */
	typedef	unsigned long ulong;					
#endif

typedef	long	lresult;

typedef __int64 int64, longlong;
typedef unsigned __int64 uint64, qword_t, ulonglong;

#ifndef BOOL
	typedef int     BOOL;
	#define TRUE  1
	#define FALSE 0
#endif

#ifndef RESULT
	#define RESULT		lresult
	#define _SUCCESS	0
	#define _ERROR		-1
#endif

#ifndef IN
#define IN
#endif

#ifndef OUT
#define OUT
#endif

#ifndef INOUT
#define INOUT
#endif

#ifndef OPTIONAL
#define OPTIONAL
#endif

#define SIZE_UUID	16

#ifdef _DEBUG
	#ifdef _TRACE
		#include <stdarg.h>
		static void DBG_TRACE(const char *fmt, ...) {
			va_list ap;
			va_start(ap, fmt);
			vprintf(fmt, ap);
			va_end(ap);
		}
    #else
		static void DBG_TRACE(const char *fmt, ...){}
	#endif
#else
	static void DBG_TRACE(const char *fmt, ...){}
#endif

/*============================================================================*/
#endif	/*UNISTD_H__*/
