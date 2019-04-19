/*
 * mThread.cpp
 *
 *  Created on: August 31, 2017
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include <iostream>

#include "mThread.hpp"

/**
* Constructor
*/
mThread::mThread(){
	#if (__cplusplus >= 201103L && !__MINGW32__)
			//std::cout << "Using std::thread (C++11)" << std::endl;
	#else
		#ifdef WIN32 // Windows
			//std::cout << "Using thread from WinAPI (< C++11)" << std::endl;
			theThread = NULL;
		#elif defined (linux) // Linux
			//std::cout << "Using pthread (< C++11)" << std::endl;
		#else // non supported plateform
			#error not defined for this platform
		#endif
	#endif
}

/**
* Destructor
*/
mThread::~mThread(){
	#if (__cplusplus >= 201103L && !__MINGW32__)
		//theThread.join();
	#else
		#ifdef WIN32 /* Windows */
			//WaitForSingleObject(theThread, INFINITE);
			if (theThread)
				CloseHandle(theThread); // to avoid memory leaks
			//FindClose(theThread);
		#elif defined (linux) /* Linux */
			//pthread_join(theThread,NULL);
		#else /* non supported plateform */
			#error not defined for this platform
		#endif
	#endif
};

/**
* Create thread
*/
void mThread::create(void (*function)(void*), void *args){
	#if (__cplusplus >= 201103L && !__MINGW32__)
		theThread = std::thread(function,args);
	#else
		#ifdef WIN32 /* Windows */
			if (theThread)
				CloseHandle(theThread); // to avoid memory leaks

			theThread = CreateThread(
					NULL,							// default security attributes
					0,								// stack size
					(LPTHREAD_START_ROUTINE) function,				// your function
					args,						// function data
					0,								// flag
					0);								// thread id
		#elif defined (linux) /* Linux */
			pthread_create(&theThread,NULL,(void* (*)(void*)) function,args);
		#else /* non supported plateform */
			#error not defined for this platform
		#endif
	#endif
};

/**
* Join
*/
void mThread::join(){
	#if (__cplusplus >= 201103L && !__MINGW32__)
		theThread.join();
	#else
		#ifdef WIN32 /* Windows */
			WaitForSingleObject(theThread, INFINITE);
		#elif defined (linux) /* Linux */
			pthread_join(theThread,NULL);
		#else /* non supported plateform */
			#error not defined for this platform
		#endif
	#endif
};
