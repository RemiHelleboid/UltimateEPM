#include "EPseudopotentialFrame.h"

#include "OptionsFrame.h"


#include <vtkAutoInit.h>


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
//VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);



#define MY_VTK_WINDOW 102

#define ID_CALCULATE 105

wxBEGIN_EVENT_TABLE(EPseudopotentialFrame, wxFrame)
EVT_MENU(ID_CALCULATE, EPseudopotentialFrame::OnCalculate)
EVT_UPDATE_UI(ID_CALCULATE, EPseudopotentialFrame::OnUpdateCalculate)
EVT_MENU(wxID_EXIT, EPseudopotentialFrame::OnExit)
EVT_MENU(wxID_PREFERENCES, EPseudopotentialFrame::OnOptions)
EVT_MENU(wxID_ABOUT, EPseudopotentialFrame::OnAbout)
EVT_TIMER(101, EPseudopotentialFrame::OnTimer)
wxEND_EVENT_TABLE()


EPseudopotentialFrame::EPseudopotentialFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame(NULL, wxID_ANY, title, pos, size),
	timer(this, 101), runningThreads(0)
{
	wxMenu *menuFile = new wxMenu;

	menuFile->Append(ID_CALCULATE, "C&alculate\tCtrl+a", "Starts computing");
	menuFile->Append(wxID_SEPARATOR);
	menuFile->Append(wxID_EXIT);

	wxMenu *menuView = new wxMenu;
	menuView->Append(wxID_PREFERENCES);

	wxMenu *menuHelp = new wxMenu;
	menuHelp->Append(wxID_ABOUT);

	wxMenuBar *menuBar = new wxMenuBar;
	menuBar->Append(menuFile, "&File");
	menuBar->Append(menuView, "&View");
	menuBar->Append(menuHelp, "&Help");

	SetMenuBar(menuBar);

	CreateStatusBar();
	SetStatusText("Welcome to Empirical Pseudopotential!");

	m_pVTKWindow = new wxVTKRenderWindowInteractor(this, MY_VTK_WINDOW);
	m_pVTKWindow->UseCaptureMouseOn();
	//m_pVTKWindow->DebugOn();
	m_pVTKWindow->DebugOff();
	ConstructVTK();

	std::vector<std::vector<double>> empty_results;
	std::vector<unsigned int> empty_pos;
	std::vector<std::string> empty_strings;
	ConfigureVTK("", empty_results, empty_pos, empty_strings);
	
	currentOptions.Open();
	currentOptions.Load();

	//Compute();
}


EPseudopotentialFrame::~EPseudopotentialFrame()
{
	StopThreads(true);
	DestroyVTK();
	if (m_pVTKWindow) m_pVTKWindow->Delete();
}


void EPseudopotentialFrame::OnOptions(wxCommandEvent& /*event*/)
{
	OptionsFrame *optionsFrame = new OptionsFrame("Options", this);
	optionsFrame->options = currentOptions;
	if (wxID_OK == optionsFrame->ShowModal())
	{
		currentOptions = optionsFrame->options;
		currentOptions.Save();
	}

	delete optionsFrame;
}


void EPseudopotentialFrame::OnExit(wxCommandEvent& /*event*/)
{
	currentOptions.Save();
	currentOptions.Close();
	Close(true);
}

void EPseudopotentialFrame::OnAbout(wxCommandEvent& /*event*/)
{
	wxMessageBox("Empirical Pseudopotential ver 1.0", "About EPseudopotential", wxOK | wxICON_INFORMATION);
}

void EPseudopotentialFrame::OnCalculate(wxCommandEvent& /*event*/)
{
	Compute();
}

void EPseudopotentialFrame::OnUpdateCalculate(wxUpdateUIEvent& event)
{
	event.Enable(isFinished());
}

void EPseudopotentialFrame::ConstructVTK()
{
	pRenderer = vtkRenderer::New();
	pContextView = vtkContextView::New();

	vtkRenderWindow *pRenderWindow = m_pVTKWindow->GetRenderWindow();
	pRenderWindow->AddRenderer(pRenderer);
	pContextView->SetInteractor(pRenderWindow->GetInteractor());
	//pContextView->GetInteractor()->Initialize();

	pChart = vtkChartXY::New();
	pChart->SetRenderEmpty(true);		
	pContextView->GetScene()->AddItem(pChart);
}

void EPseudopotentialFrame::ConfigureVTK(const std::string& name, const std::vector<std::vector<double>>& results, std::vector<unsigned int>& symmetryPointsPositions, std::vector<std::string>& symmetryPointsLabels)
{
	pChart->ClearPlots();

	if (name.empty()) pChart->SetTitle("Band");
	else
	{
		std::string Name = name;
		Name += " Band";
		pChart->SetTitle(Name.c_str());
	}


	pChart->SetAutoAxes(false);

	pChart->GetAxis(vtkAxis::BOTTOM)->SetTitle("k - Symmetry Points");
	pChart->GetAxis(vtkAxis::LEFT)->SetTitle("Energy (eV)");

	if (results.empty()) return;

	int numPoints = results.size();


	// set up the data table

	vtkNew<vtkTable> table;

	vtkNew<vtkFloatArray> arrX;
	arrX->SetName("X");
	table->AddColumn(arrX.GetPointer());

	for (int i = 0; i < results[0].size(); ++i)
	{
		vtkNew<vtkFloatArray> arrC;
		arrC->SetName((std::string("Y") + std::to_string(i)).c_str());
		table->AddColumn(arrC.GetPointer());
	}

	table->SetNumberOfRows(numPoints);


	// set values for X axis column
	for (int j = 0; j < numPoints; ++j)
		table->SetValue(j, 0, j);


	// now the lines in the chart, each has its own column in the table
	for (int j = 0; j < numPoints; ++j)
		for (int i = 0; i < results[j].size(); ++i)
		  table->SetValue(j, i + 1, results[j][i]);

	// add the lines to the chart
	for (int i = 0; i < results[0].size(); ++i)
	{
	  vtkPlot *line = pChart->AddPlot(vtkChart::LINE);
	  // Use columns 0 and 1 for x and y
	  line->SetInputData(table.GetPointer(), 0, i + 1);
	  
	  // Make the plot blue, with a width of 2.0 pixels
	  line->SetColor(0, 0, 255, 255);
	  line->SetWidth(2.0);	
	}

	vtkAxis *xAxis = pChart->GetAxis(vtkAxis::BOTTOM);

	vtkNew<vtkDoubleArray> posArray;
	vtkNew<vtkStringArray> labelsArray;

	

	posArray->SetName("X");
	labelsArray->SetName("X");

	posArray->SetNumberOfValues(symmetryPointsPositions.size() + 1);
	labelsArray->SetNumberOfValues(symmetryPointsPositions.size() + 1);

	unsigned int index = 0;
	for (unsigned int pos : symmetryPointsPositions)
	{
		posArray->SetValue(index, pos);
		labelsArray->SetValue(index, symmetryPointsLabels[index]);
		++index;
	}
	// add the last

	posArray->SetValue(index, numPoints - 1);
	labelsArray->SetValue(index, symmetryPointsLabels[index]);

	xAxis->SetCustomTickPositions(posArray.GetPointer(), labelsArray.GetPointer());	
}

void EPseudopotentialFrame::DestroyVTK()
{
	if (pChart) pChart->Delete();
	if (pRenderer) pRenderer->Delete();
	if (pContextView) pContextView->Delete();
}


void EPseudopotentialFrame::OnTimer(wxTimerEvent& WXUNUSED(event))
{
	if (isFinished())
	{
		timer.Stop();
		StopThreads();

		SetTitle("Finished - EPseudopotential");

		Refresh();
	}
}

void EPseudopotentialFrame::Compute()
{
	if (!isFinished()) return;

	wxBeginBusyCursor();

	computeOptions = currentOptions;

	SetTitle("Computing - EPseudopotential");

	bandStructure.symmetryPoints.path = computeOptions.paths[computeOptions.pathNo];
	bandStructure.Initialize(computeOptions.nrPoints, computeOptions.nearestNeighbors);

	unsigned int nrThreads = computeOptions.nrThreads;
	if (0 == nrThreads) computeOptions.nrThreads = nrThreads = 1;

	runningThreads = nrThreads;

	const unsigned int nrPoints = bandStructure.GetPointsNumber();
	const unsigned int interval =  nrPoints / nrThreads;

	unsigned int startPoint = 0;
	for (unsigned int i = 0; i < nrThreads; ++i)
	{
		threadsList.push_back(std::make_unique<EPThread>(computeOptions, this, startPoint, std::min(startPoint + interval, nrPoints)));
		startPoint += interval;

		threadsList.back()->Start();
	}

	timer.Start(100);
}

bool EPseudopotentialFrame::isFinished() const
{
	return 0 == runningThreads;
}


void EPseudopotentialFrame::StopThreads(bool cancel)
{
	if (cancel)
		for (auto &thrd : threadsList) thrd->Terminate();

	for (auto &thrd : threadsList)
	{
		thrd->join();

		if (!cancel) bandStructure.results.insert(bandStructure.results.end(), thrd->results.begin(), thrd->results.end());
	}
	threadsList.clear();

	SetTitle("Finished - EPseudopotential");


	if (!cancel)
	{
		bandStructure.AdjustValues();		
		ConfigureVTK(std::string(computeOptions.materialName.c_str()), bandStructure.results, bandStructure.symmetryPointsPositions, bandStructure.symmetryPoints.path);
	}

	if (wxIsBusy()) wxEndBusyCursor();
}