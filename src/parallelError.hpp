#ifndef PARALLELERROR_HPP
#define PARALLELERROR_HPP

#include <string>
#include <stdexcept>
#include <atomic>

using namespace std;

/*
* errors inside an OpenMP parallel loop: an exception must not leave a parallel region (otherwise the program aborts),
* so each iteration catches it and stores the first error here; the remaining iterations are skipped (failed()),
* and rethrow() after the loop stops SNEEP with this error (main: message and exit code 1)
*
* usage:
*	ParallelError parallelError;
*	#pragma omp parallel for
*	for (...){
*		if (parallelError.failed()) continue;
*		try{ ... }catch (const exception& e){ parallelError.set(e.what()); }
*	}
*	parallelError.rethrow();
*/
class ParallelError{

	public:
	bool failed(){
		return failed_;
	}
	void set(const string& message){
		#pragma omp critical (parallelError)
		{
			if (!failed_){
				message_ = message;
				failed_ = true;
			}
		}
	}
	void rethrow(){
		if (failed_){
			throw runtime_error(message_);
		}
	}

	private:
	atomic<bool> failed_{false};
	string message_ = "";
};

#endif/*PARALLELERROR_HPP*/
