#include "Options.h"

#include <wx/stdpaths.h> 


Options::Options()
	: nrThreads(4), materialName("Si"), nrPoints(600), nearestNeighbors(10), nrLevels(10), pathNo(0),
	paths{ { {"K", "W", "X", "G", "L", "W"}, {"L", "G", "X", "K", "G" }, {"W", "G", "X", "W", "L", "G"},  {"L", "G", "X", "U", "K", "G"} } },
	m_fileconfig(nullptr)
{
}


Options::~Options()
{
}

void Options::Open()
{
	if (m_fileconfig) return;

	wxString dir = wxStandardPaths::Get().GetConfigDir() + wxFileName::GetPathSeparator();

	if(!wxFileName::DirExists(dir))
		wxFileName::Mkdir(dir, 0777, wxPATH_MKDIR_FULL);

	wxString iniFilePath = dir + "EPseudopotential.ini";

	m_fileconfig = new wxFileConfig("EPseudopotential", wxEmptyString, iniFilePath);

	wxConfigBase::Set(m_fileconfig);
}

void Options::Close()
{
	delete m_fileconfig;
	m_fileconfig = NULL;
	wxConfigBase::Set(NULL);
}

void Options::Load()
{
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		nrThreads = conf->ReadLong("/nrThreads", 4);
		materialName = conf->Read("/materialName", "Si");
		nrPoints = conf->ReadLong("/nrPoints", 600);
		nearestNeighbors = conf->ReadLong("/nearestNeighbors", 10);
		nrLevels = conf->ReadLong("/nrLevels", 10);
		pathNo = conf->ReadLong("/pathNo", 0);
	}
}

void Options::Save()
{
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		conf->Write("/nrThreads", static_cast<long int>(nrThreads));
		conf->Write("/materialName", materialName);
		conf->Write("/nrPoints", static_cast<long int>(nrPoints));
		conf->Write("/nearestNeighbors", static_cast<long int>(nearestNeighbors));
		conf->Write("/nrLevels", static_cast<long int>(nrLevels));
		conf->Write("/pathNo", static_cast<long int>(pathNo));
	}

	if (m_fileconfig)
		m_fileconfig->Flush();
}
