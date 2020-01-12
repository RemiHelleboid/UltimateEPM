#pragma once

#include <string>
#include <vector>

#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/fileconf.h>

class Options
{
public:
	Options();
	~Options()
	{
		delete m_fileconfig;
	}

	// avoid double deletion of m_fileconfig at destruction if copied
	Options(const Options& other)
		: 
		nrThreads(other.nrThreads),
		materialName(other.materialName),
		nrPoints(other.nrPoints),
		nearestNeighbors(other.nearestNeighbors),
		nrLevels(other.nrLevels),
		pathNo(other.pathNo),
		paths(other.paths),
		m_fileconfig(nullptr)
	{
	}

	Options& operator=(const Options& other)
	{
		nrThreads = other.nrThreads;
		materialName = other.materialName;
		nrPoints = other.nrPoints;
		nearestNeighbors = other.nearestNeighbors;
		nrLevels = other.nrLevels;
		pathNo = other.pathNo;
		paths = other.paths;
		m_fileconfig = nullptr;

		return *this;
	}

	void Load();
	void Save();

	int nrThreads;

	wxString materialName;
	int nrPoints;
	int nearestNeighbors;
	int nrLevels;

	int pathNo;

	std::vector<std::vector<std::string>> paths;

protected:
	void Open();
	void Close();

	wxFileConfig *m_fileconfig;
};

