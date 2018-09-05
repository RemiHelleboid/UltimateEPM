#include "EPseudopotentialApp.h"
#include "EPseudopotentialFrame.h"


wxIMPLEMENT_APP(EPseudopotentialApp);


bool EPseudopotentialApp::OnInit()
{
	if (!wxApp::OnInit())
		return false;
	
	frame = new EPseudopotentialFrame("EPseudopotential", wxPoint(50, 50), wxSize(1024, 800));
	frame->Show(true);

	return true;
}
