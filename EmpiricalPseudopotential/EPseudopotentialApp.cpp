#include "EPseudopotentialApp.h"
#include "EPseudopotentialFrame.h"


wxIMPLEMENT_APP(EPseudopotentialApp);


bool EPseudopotentialApp::OnInit()
{
	frame = new EPseudopotentialFrame("EPseudopotential", wxPoint(50, 50), wxSize(1024, 800));
	frame->Show(true);

	return true;
}
