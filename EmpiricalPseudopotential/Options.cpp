#include "Options.h"

#include <iostream>
// #include <wx/stdpaths.h>

Options::Options()
    : nrThreads(4),
      materialName("Si"),
      nrPoints(400),
      nearestNeighbors(10),
      nrLevels(10),
      pathNo(1),
      paths{{
          {"K", "W", "X", "G", "L", "W"},
          {"W", "G", "X", "W", "L", "G"},
          {"W", "L", "G", "X", "W", "K"},
          {"L", "G", "X", "K", "G"},
          {"L", "G", "X", "W", "K", "G"},
          {"L", "G", "X", "U", "K", "G"},
          {"L", "G", "X", "U", "G"},
          {"L", "K", "W", "G", "X", "W", "L", "G", "K"},
          {"G", "X", "K", "G", "L", "K", "W", "X"},
          {"G", "X", "W", "L", "G", "K", "W", "U"},
          {"G", "X", "W", "L", "G", "K"},
          {"G", "X", "W", "L", "G", "K", "X"},
          {"G", "X", "W", "L", "G", "X"},
          {"G", "X", "W", "G", "U", "X"},
          {"G", "X", "W", "K", "L", "G"},
          {"G", "X", "W", "K", "G", "L", "U", "W", "L", "K"},
          {"G", "X", "U", "K", "G", "L", "W", "X"},
      }},
      m_fileconfig(nullptr) {}

void Options::print_options() {
    std::cout << "Options:" << std::endl;
    std::cout << "nrThreads: " << nrThreads << std::endl;
    std::cout << "materialName: " << materialName << std::endl;
    std::cout << "nrPoints: " << nrPoints << std::endl;
    std::cout << "nearestNeighbors: " << nearestNeighbors << std::endl;
    std::cout << "nrLevels: " << nrLevels << std::endl;
    std::cout << "pathNo: " << pathNo << std::endl;
    std::cout << "paths: " << std::endl;
    for (auto& path : paths) {
        for (auto& point : path) {
            std::cout << point << " ";
        }
        std::cout << std::endl;
    }
}

// void Options::Open()
// {
// 	if (m_fileconfig) return;

// 	wxString dir = wxStandardPaths::Get().GetConfigDir() + wxFileName::GetPathSeparator();

// 	if(!wxFileName::DirExists(dir))
// 		wxFileName::Mkdir(dir, 0777, wxPATH_MKDIR_FULL);

// 	wxString iniFilePath = dir + "EPseudopotential.ini";

// 	m_fileconfig = new wxFileConfig("EPseudopotential", wxEmptyString, iniFilePath);

// 	wxConfigBase::Set(m_fileconfig);
// }

// void Options::Close()
// {
// 	delete m_fileconfig;
// 	m_fileconfig = NULL;
// 	wxConfigBase::Set(NULL);
// }

// void Options::Load()
// {
// 	Open();
// 	wxConfigBase *conf=wxConfigBase::Get(false);
// 	if (conf)
// 	{
// 		nrThreads = conf->ReadLong("/nrThreads", 4);
// 		materialName = conf->Read("/materialName", "Si");
// 		nrPoints = conf->ReadLong("/nrPoints", 600);
// 		nearestNeighbors = conf->ReadLong("/nearestNeighbors", 10);
// 		nrLevels = conf->ReadLong("/nrLevels", 10);
// 		pathNo = conf->ReadLong("/pathNo", 0);
// 	}
// 	Close();
// }

// void Options::Save()
// {
// 	Open();
// 	wxConfigBase *conf=wxConfigBase::Get(false);
// 	if (conf)
// 	{
// 		conf->Write("/nrThreads", static_cast<long int>(nrThreads));
// 		conf->Write("/materialName", materialName);
// 		conf->Write("/nrPoints", static_cast<long int>(nrPoints));
// 		conf->Write("/nearestNeighbors", static_cast<long int>(nearestNeighbors));
// 		conf->Write("/nrLevels", static_cast<long int>(nrLevels));
// 		conf->Write("/pathNo", static_cast<long int>(pathNo));
// 	}

// 	if (m_fileconfig)
// 		m_fileconfig->Flush();
// 	Close();
// }
