#pragma once


#include <thread>
#include <vector>
#include <atomic>

#include "Options.h"

class EPseudopotentialFrame;

class EPThread
{
public:
	EPThread(const Options& options, EPseudopotentialFrame* frame, unsigned int startPoint, unsigned int endPoint);
	~EPThread();


	void Start();
	void join();

	void Terminate();


	const Options& m_options;

	EPseudopotentialFrame* m_frame;

	unsigned int m_startPoint;
	unsigned int m_endPoint;


	std::vector<std::vector<double>> results;

protected:
	void Calculate();

	std::thread mThread;
	std::atomic_bool terminate;
};

