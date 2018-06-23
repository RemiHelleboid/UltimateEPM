#include "EPThread.h"


#include "EPseudopotentialFrame.h"

EPThread::EPThread(const Options& options, EPseudopotentialFrame* frame, unsigned int startPoint, unsigned int endPoint)
	: m_options(options), m_frame(frame), m_startPoint(startPoint), m_endPoint(endPoint), terminate(false)
{
}


EPThread::~EPThread()
{
	join();
}

void EPThread::Start()
{
	mThread = std::thread([this]() {
		Calculate();		
	});
}

void EPThread::join()
{
	if (mThread.joinable()) mThread.join();
}


void EPThread::Calculate()
{
	EmpiricalPseudopotential::Material &mat = m_frame->bandStructure.materials.materials[std::string(m_options.materialName.c_str())];

	results = m_frame->bandStructure.Compute(mat, m_startPoint, m_endPoint, m_options.nrLevels, terminate);

	--m_frame->runningThreads;
}

void EPThread::Terminate()
{
	terminate = true;
}