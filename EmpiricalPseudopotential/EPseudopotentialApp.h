#pragma once



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


class EPseudopotentialFrame;

class EPseudopotentialApp : public wxApp
{
public:
	EPseudopotentialFrame * frame = nullptr;

	virtual bool OnInit();
};

