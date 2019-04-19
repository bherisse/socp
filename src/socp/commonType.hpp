/*
 * commonType.hpp
 *
 *  Created on: June 14, 2018
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#ifndef real
	#ifdef doubleType
		typedef double real;
		#define __cminpack_double__	// floating-point type is "double" for cminpack
	#elif floatType
		typedef float real;
		#define __cminpack_float__	// floating-point type is "float" for cminpack
	#else
		typedef double real; 
		#define __cminpack_double__	// floating-point type is "double" for cminpack
	#endif
#endif