//*****************************************************************************
#ifndef __DEBUGTOOLS_H
#define __DEBUGTOOLS_H 1

/**
@file DebugTools.h
******************
@brief DebugTools.h is a library of debugging tools.

DebugTools.h library contains various tools gathered in a single namespace
called debugtools. These tools can help debugging programs.

@see val::debugtools
*/

/*
// Uncomment this section to enable debug mode
#define _DEBUG
#define _DEBUG_WARNINGS
#define _DEBUG_ERRORS
#define _DEBUG_INSTANCES
*/

#include <iostream>
#include <cassert>

namespace val
{
/**
@namespace val::debugtools
**************************
@brief debugtools is a namespace for debuging tools

debugtools namespace contains various tools to help the debugging process of
programs. It provide tools to trace code and classes and also tools to track
the behaviour of a program. This set of tools has been designed to be executed
only in specific debug modes that the user can enable or disable depending
his/her expectations. It means that when all debug modes are turned off, these
debugging instructions will just be ignored and not compiled so the released
program will not contain them (saves speed and space). With these tools you
can leave your debug messages in your source code without worrying how you
will remove them on release time: you'll just have to disable the macros with
compiler flags! And you won't have to put them back when you'll need them
again!
Please refer to the documentation of each debug macro to learn out how to
enable or disable them and what they should be used for.

To enable every debug instruction, you have to enable their debug flags.
There are several ways to enable these flags. The "cleanest" way to activate
these flags is using the compiler command line. To do so, refer to your
compiler documentation. You have to find wich option defines macro and then
you will have to define the following macros:
- _DEBUG
- _DEBUG_INSTANCES
- _DEBUG_WARNINGS
- _DEBUG_ERRORS

For example, with linux g++, you have to add the following options:
@code
g++ -D_DEBUG -D_DEBUG_INSTANCES -D_DEBUG_WARNINGS -D_DEBUG_ERRORS
@endcode
With Visual C++, you have to go in the menu "Project"->"Properties", then
select the appropriate configuration, then select "C/C++" -> "command line"
and add in the additionnal options field the following commands:
@code
/D _DEBUG_INSTANCES /D _DEBUG_WARNINGS /D _DEBUG_ERRORS
@endcode

If you need to use the debugging tools in one file, you should try (where you
include "DebugTools.hpp" in that file):
@code
#define _DEBUG
#define _DEBUG_INSTANCES
#define _DEBUG_WARNINGS
#define _DEBUG_ERRORS
#include "DebugTools.hpp"
#undef _DEBUG_ERRORS
#undef _DEBUG_WARNINGS
#undef _DEBUG_INSTANCES
#undef _DEBUG
@endcode
You can do the opposite if you don't want to debug instances in one file:
@code
#undef _DEBUG_INSTANCES
#include "DebugTools.hpp"
#define _DEBUG_INSTANCES
@endcode

You can also edit DebugTools.h and uncomment the following lines in the
beginning of the file but that wouldn't be really clean ;-) :
@code
// Uncomment this section to enable debug mode
#define _DEBUG
#define _DEBUG_WARNINGS
#define _DEBUG_ERRORS
#define _DEBUG_INSTANCES
@endcode

@note to use debug instructions you must enable them using the debug flags as
 explained int the documentation.

@note macros are just like global instructions and doesn't really belong to
 the namespace so you must not prefix them with any "val::debugtools::".

@see DebugTools.h
@see TRACED
@see TRACEI
@see TRACEW
@see TRACER
@see ASSERT
*/
namespace debugtools
{

	//! String prefixing any trace debugging message using TRACED.
	const char SZ_TRACED_MSG_PREFIX[]		= "";
	//! String prefixing any instance debugging message using TRACEI.
	const char SZ_INSTANCES_MSG_PREFIX[]	= "INSTANCE: ";
	//! String prefixing any warning debugging message using TRACEW.
	const char SZ_WARNINGS_MSG_PREFIX[]		= "WARNING: ";
	//! String prefixing any error debugging message using TRACER.
	const char SZ_ERRORS_MSG_PREFIX[]		= "ERROR: ";
	
}; // namespace debugtools
}; // namespace val
#endif //ifndef __DEBUGTOOLS_H

//*****************************************************************************
// MACROS DEFINITIONS (ONLY)
//*****************************************************************************
// note: anything that is not a macro should be added in the "debugtools
// namespace" section, before the "#endif //ifndef __DEBUGTOOLS_H" instruction

// undefines macro before (re)defining them.
#ifdef TRACED
	#undef TRACED
#endif
#ifdef TRACEI
	#undef TRACEI
#endif
#ifdef TRACEW
	#undef TRACEW
#endif
#ifdef TRACER
	#undef TRACER
#endif

/**
@def TRACED
***********
@brief TRACED outputs debugging messages the way std::cout does.

TRACED outputs debugging messages to stdout using the standard std::cout
function. The starting "<<" of the stream can be omitted and the macro always
appends an endline code after the stream and flushes the stream to make sure it
is displayed at the right time on the standard output (std::stdout).
To enable this function, "_DEBUG" must be defined in your program
(#define _DEBUG).

Example:
 @code
 #define _DEBUG
 #include "DebugTools.hpp"
 ...
 int i = 7;
 TRACED("This message will be displayed on standard output and the value of i is " << i << "!");
 TRACED("Quite simple and flexible!");
 @endcode
 Output:
 @code
 This message will be displayed on standard output and the value of i is 7!
 Quite simple and flexible!
 @endcode

@note when "_DEBUG" is not defined, TRACED macro is seen by the compiler as no
 instruction (not compiled) and TRACED will do nothing.
@note MS Visual C++ users: this macro is different form the one provided by
 Visual C++!
@note MS Visual C++ auto-defines "_DEBUG" only in "Debug" configuration and
 not in "Release" configuration.
@note there are other tracing instructions for specific usage.

@see TRACEI
@see TRACEW
@see TRACER

@param debug_stream:
 @c debug_stream is the same kind of stream one puts after std::cout except the
 first "<<" can be omitted and the macro always appends a std::endl after the
 stream.
*/
#ifdef _DEBUG
	#define TRACED(debug_stream) std::cout << val::debugtools::SZ_TRACED_MSG_PREFIX << debug_stream << std::endl; std::cout.flush()
#else
	#define TRACED(debug_stream)
#endif

/**
@def TRACEI
***********
@brief TRACEI outputs debugging messages the way std::cout does for instances.

The way TRACEI works is just like TRACED; the two only differences are that
TRACEI prefixes its message using SZ_INSTANCES_MSG_PREFIX string and it is
enabled and disabled by defining or not "_DEBUG_INSTANCES"
(#define _DEBUG_INSTANCES).

Example:
 @code
 #define _DEBUG
 #define _DEBUG_INSTANCES
 #include "DebugTools.hpp"
 #include "MyClass.hpp"
 CMyClass::CMyClass(void)
 {
 	TRACEI("CMyClass::CMyClass(void)");
 	TRACED("The message above lets you know when the default constructor is called.");
 }
 CMyClass::~CMyClass(void)
 {
 	TRACEI("CMyClass::~CMyClass(void)");
 }
 @endcode
 Output (when _DEBUG and _DEBUG_INSTANCES are defined):
 @code
 INSTANCE: CMyClass::CMyClass(void)
 The message above lets you know when the default constructor is called.
 INSTANCE: CMyClass::~CMyClass(void)
 @endcode

@note when "_DEBUG_INSTANCES" is not defined, TRACEI macro is seen by the
 compiler as no instruction (not compiled) and TRACEI will do nothing.

@see TRACED
@see TRACEW
@see TRACER

@param debug_stream:
 @c debug_stream is the same kind of stream one puts after std::cout except the
 first "<<" can be omitted and the macro always appends a std::endl after the
 stream.
*/
#ifdef _DEBUG_INSTANCES
	#define TRACEI(debug_stream) std::cout << val::debugtools::SZ_INSTANCES_MSG_PREFIX << debug_stream << std::endl; std::cout.flush()
#else
	#define TRACEI(debug_stream)
#endif


/**
@def TRACEW
***********
@brief TRACEW outputs debugging messages the way std::cout does for warnings.

The way TRACEW works is just like TRACED; the two only differences are that
TRACEW prefixes its message using SZ_WARNINGS_MSG_PREFIX string and it is
enabled and disabled by defining or not "_DEBUG_WARNINGS"
(#define _DEBUG_WARNINGS).

Example:
 @code
 #define _DEBUG
 #define _DEBUG_WARNINGS
 #include "DebugTools.hpp"
 void MyFunc(int &iTest)
 {
	if (0 > iTest)
	{
 		TRACEW("\"iTest\" wasn't positive (" << iTest << ") and has been set to 0!");
		iTest = 0;
	}
	else
	{
		TRACED("OK: Nothing wrong with \"iTest\"=" << iTest);
	}
	//...
 }
 @endcode
 Output:\n
 if "iTest" has a positive or null value like "2":
 @code
 OK: Nothing wrong with "iTest"=2
 @endcode
 if "iTest" has a negative value like "-3":
 @code
 WARNING: "iTest" wasn't positive (-3) and has been set to 0!
 @endcode

@note when "_DEBUG_WARNINGS" is not defined, TRACEW macro is seen by the
 compiler as no instruction (not compiled) and TRACEW will do nothing.

@see TRACED
@see TRACEI
@see TRACER

@param debug_stream:
 @c debug_stream is the same kind of stream one puts after std::cout except the
 first "<<" can be omitted and the macro always appends a std::endl after the
 stream.
*/
#ifdef _DEBUG_WARNINGS
	#define TRACEW(debug_stream) std::cout << val::debugtools::SZ_WARNINGS_MSG_PREFIX << debug_stream << std::endl; std::cout.flush()
#else
	#define TRACEW(debug_stream)
#endif

/**
@def TRACER
***********
@brief TRACER outputs debugging messages the way std::cout does for errors.

The way TRACER works is just like TRACED; the two only differences are that
TRACER prefixes its message using SZ_ERRORS_MSG_PREFIX string and it is
enabled and disabled by defining or not "_DEBUG_ERRORS" (#define _DEBUG_ERRORS).

Example:
 @code
 #define _DEBUG
 #define _DEBUG_ERRORS
 #include "DebugTools.hpp"
 void MyFunc(int &iTest)
 {
	if (0 > iTest)
	{
 		TRACER("\"iTest\" wasn't positive (" << iTest << "), function aborted!");
		return;
	}
	else
	{
		TRACED("OK: Nothing wrong with \"iTest\"=" << iTest);
	}
	//...
 }
 @endcode
 Output:\n
 if "iTest" has a positive or null value like "2":
 @code
 OK: Nothing wrong with "iTest"=2
 @endcode
 if "iTest" has a negative value like "-3":
 @code
 ERROR: "iTest" wasn't positive (-3), function aborted!
 @endcode

@note when "_DEBUG_ERRORS" is not defined, TRACER macro is seen by the
 compiler as no instruction (not compiled) and TRACER will do nothing.

@see TRACED
@see TRACEI
@see TRACEW

@param debug_stream:
 @c debug_stream is the same kind of stream one puts after std::cout except the
 first "<<" can be omitted and the macro always appends a std::endl after the
 stream.
*/
#ifdef _DEBUG_ERRORS
	#define TRACER(debug_stream) std::cout << val::debugtools::SZ_ERRORS_MSG_PREFIX << debug_stream << std::endl; std::cout.flush()
#else
	#define TRACER(debug_stream)
#endif

/**
@def ASSERT
***********
@brief ASSERT is a wrapper for the standard "std::assert" instruction.

ASSERT does the same job as std::assert. However ASSERT will only do its job by
simply calling std::assert when "_DEBUG" is defined. Otherwise it will be seen
by the compiler as no instruction (not compiled). ASSERT checks if the result
of an expression is true (not NULL). If it's false, ASSERT prints a diagnostic
message and aborts the program. To enable this function, "_DEBUG" must be
defined in your program (#define _DEBUG).

Example:
 @code
 int i = 0;
 ASSERT(0 == i); // will do nothing as the test is OK
 ASSERT(i); // will stop the program execution as i==false and display a message
 @endcode

@note when "_DEBUG" is not defined, ASSERT macro is seen by the compiler as
 no instruction (not compiled) and ASSERT will do nothing.
@note if the ASSERT macro is already defined, it won't be redefined.
@note MS Visual C++ auto-defines "_DEBUG" only in "Debug" configuration and
 not in "Release" configuration.
@note MS Visual C++ users: this macro works the same way the one provided by
 Visual C++!

@param assert_expression:
 @c assert_expression is a logic expression that should return true (a non-null
 value) otherwise ASSERT stops the program.
*/
#ifdef _DEBUG
	#ifndef ASSERT
		#define ASSERT(assert_expression) std::assert(assert_expression)
	#endif
#else
	#ifndef ASSERT
		#define ASSERT(assert_expression)
	#endif
#endif
// END MACROS DEFINITIONS (ONLY)
//*****************************************************************************
