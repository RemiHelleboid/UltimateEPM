#pragma once

#define _MATH_DEFINES_DEFINED

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
	#pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers
#ifndef WX_PRECOMP
	#include "wx/wx.h"
#endif

#include <atomic>
#include <list>


//#include "Options.h"

#include "wxVTKRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"

#include "vtkGraph.h"
#include "vtkGraphItem.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"

#include "vtkTable.h"
#include "vtkArray.h"
#include "VtkFloatArray.h"
#include "vtkNew.h"
#include "vtkChart.h"
#include "vtkChartXY.h"

#include "vtkStringArray.h"
#include "vtkDoubleArray.h"
#include "vtkPlot.h"
#include "vtkAxis.h"

#include "Options.h"
#include "BandStructure.h"
#include "EPThread.h"

class EPseudopotentialFrame : public wxFrame
{
public:
	EPseudopotentialFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
	~EPseudopotentialFrame();

private:
	void OnCalculate(wxCommandEvent& event);
	void OnUpdateCalculate(wxUpdateUIEvent& event);

	void OnExit(wxCommandEvent& event);
	void OnOptions(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnTimer(wxTimerEvent& event);

protected:
	void ConstructVTK();
	void ConfigureVTK(const std::string& name, const std::vector<std::vector<double>>& results, std::vector<unsigned int>& symmetryPointsPositions, std::vector<std::string>& symmetryPointsLabels);
	void DestroyVTK();


	bool isFinished() const;
	void StopThreads(bool cancel = false);

private:
	wxVTKRenderWindowInteractor *m_pVTKWindow;

	// vtk classes
	vtkRenderer     *pRenderer;
	vtkContextView	*pContextView;

	vtkChartXY *pChart;

	wxTimer timer;

	std::list<std::unique_ptr<EPThread>> threadsList;

	Options computeOptions; // what's actually displayed
public:
	std::atomic_int runningThreads;

	EmpiricalPseudopotential::BandStructure bandStructure;

	Options currentOptions; // what's edited

private:

	wxDECLARE_EVENT_TABLE();
public:
	void Compute();
};
