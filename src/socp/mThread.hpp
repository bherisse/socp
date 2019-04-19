/*
 * mThread.hpp
 *
 *  Created on: August 31, 2017
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include <iostream>

#if (__cplusplus >= 201103L && !__MINGW32__)
	#include <thread>         // std::thread
	typedef std::thread hThread;
#else
	#ifdef WIN32 /* Windows */
		#include <windows.h>
		typedef HANDLE hThread;
	#elif defined (linux) /* Linux */
		#include <pthread.h>
		typedef pthread_t hThread;
	#else /* non supported plateform */
		#error not defined for this platform
	#endif
#endif

#ifndef _MTHREAD_H_
#define _MTHREAD_H_


/**************			mThread class			******************************/
class mThread
{

public:

	/**
	* Constructor
	*/
	mThread();

	/**
	* Destructor
	*/
	~mThread();

	/**
	* Create Thread
	* @param function the function to be executed by the thread
	* @param args a pointer to user arguments
	*/
	void create(void (*function)(void*), void *args);

	/**
	* Join
	*/
	void join();


private:

	hThread theThread;		///< the Thread handle

	//#if __cplusplus >= 201103L
	//	std::thread theThread;
	//#else
	//	#ifdef WIN32 /* Windows */
	//		HANDLE  theThread;
	//	#elif defined (linux) /* Linux */
	//		pthread_t theThread;
	//	#else /* non supported plateform */
	//		#error not defined for this platform
	//	#endif
	//#endif

};

#endif //_MTHREAD_H_
