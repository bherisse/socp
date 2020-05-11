/*
 * mThread.hpp
 *
 *  Created on: August 31, 2017
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include <iostream>

#if (__cplusplus >= 201103L && !__MINGW32__)
	#include <thread>         // std::thread
	#include <future>
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
	mThread() {
		#if (__cplusplus >= 201103L && !__MINGW32__)
		#else
			#ifdef WIN32 // Windows
				theThread = NULL;
			#elif defined (linux) // Linux
			#else // non supported plateform
				#error not defined for this platform
			#endif
		#endif
	};

	/**
	* Destructor
	*/
	~mThread() {
		#if (__cplusplus >= 201103L && !__MINGW32__)
		#else
			#ifdef WIN32 /* Windows */
				if (theThread)
					CloseHandle(theThread); // to avoid memory leaks
			#elif defined (linux) /* Linux */
			#else /* non supported plateform */
				#error not defined for this platform
			#endif
		#endif
	};

	/**
	* Create Thread
	* @param function the function to be executed by the thread
	* @param args a pointer to user arguments
	*/
	void create(void (*function)(void*), void *args) {
		#if (__cplusplus >= 201103L && !__MINGW32__)
			theThread = std::thread(function, args);
		#else
			#ifdef WIN32 /* Windows */
				if (theThread)
					CloseHandle(theThread); // to avoid memory leaks

				theThread = CreateThread(
					NULL,							// default security attributes
					0,								// stack size
					(LPTHREAD_START_ROUTINE)function,				// your function
					args,							// function data
					0,								// flag
					0);								// thread id
			#elif defined (linux) /* Linux */
					pthread_create(&theThread, NULL, (void* (*)(void*)) function, args);
			#else /* non supported plateform */
				#error not defined for this platform
			#endif
		#endif
	};

	/**
	* Join
	*/
	void join() {
		#if (__cplusplus >= 201103L && !__MINGW32__)
			theThread.join();
		#else
			#ifdef WIN32 /* Windows */
				WaitForSingleObject(theThread, INFINITE);
			#elif defined (linux) /* Linux */
				pthread_join(theThread, NULL);
			#else /* non supported plateform */
				#error not defined for this platform
			#endif
		#endif
	};

	/**
	* Join with timeout
	*/
	int join_for(int timeout) {
		#if (__cplusplus >= 201103L && !__MINGW32__)
			mFuture = std::async(std::launch::async, &std::thread::join, &theThread);
			if (mFuture.wait_for(std::chrono::milliseconds(timeout)) == std::future_status::timeout) {
				std::cout << std::endl << " ! The thread has not terminated within " << timeout << " ms" << std::endl;
								
				return -1;
			}
		#else
			#ifdef WIN32 /* Windows */
				if (WaitForSingleObject(theThread, (DWORD)timeout) == WAIT_TIMEOUT) {
					std::cout << std::endl << " ! The thread has not terminated within " << timeout << " ms" << std::endl;

					return -1;
				}
			#elif defined (linux) /* Linux */
				timespec timeStruct;
				clock_gettime(CLOCK_REALTIME, &timeStruct);
				int timeout_sec = timeout/1000;
				timeStruct.tv_sec += timeout_sec;
				timeStruct.tv_nsec = (timeout-timeout_sec*1000)*1e6;
				if (pthread_timedjoin_np(theThread, NULL, &timeStruct) == ETIMEDOUT) {
					std::cout << std::endl << " ! The thread has not terminated within " << timeout << " ms" << std::endl;

					return -1;
				}
			#else /* non supported plateform */
				#error not defined for this platform
			#endif
		#endif

		return 0;
	};


private:

	hThread theThread;		///< the Thread handle
	#if (__cplusplus >= 201103L && !__MINGW32__)
		std::future<void> mFuture;
	#endif

};

#endif //_MTHREAD_H_
