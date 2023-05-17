#pragma once
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cstdint>
#include <stdlib.h>
#include <stdio.h>
#include <json.hpp>
#include <chrono>
#define EXPORT extern "C" __declspec(dllexport)
using namespace std;
using namespace Eigen;
using json = nlohmann::json;
constexpr char PATH[] = "D:\\statistics\\samples\\114557_869.bin"; //anno 없는데 wgs있음
constexpr char PATH2[] = "D:\\statistics\\183512_090_abnormal.bin";


struct _ST_RCP_FFT
{
	int fftWindowIdx;
	int viewMode;
	int dbvRange;
	double freqStart;
	double freqEnd;
};

struct _ST_WINDOW
{
	int windowDataIdx;
	int statistic;
	std::vector<double> params;
	_ST_RCP_FFT fftInfos;

	std::string startType;
	double startTime;
	std::string startVS;

	std::string endType;
	double endTime;
	double endLength;
	std::string endVS;

	int period;

	double* datas;
};

struct _ST_HARDWARE
{
	bool enable;
	std::string daqModel;
	int couplingType;
	int voltageRange;
	int srcCnt;
};

struct _ST_SOFTWARE
{
	bool bRawSave;
	int rawDataStoragePeriod;
	std::string rawDataPath;

	bool bWindowSave;
	bool bWindowRawSave;
	int windowStoragePeriod;
	std::string windowDataPath;
};
enum class ENUM_WINDOWS_INDEX
{
	rectangle,
	Hann,
	Hamming,
	Blackman,
	Flat_top
};

enum class ENUM_VIEW_MODE
{
	Linear,
	DBV
};
enum class EStatistics
{
	RMS = 0, 
	Mean,
	MeanH,
	MeanG,
	StDev,
	Skew,
	Kurt,
	Mode,
	Median,
	Q1,
	Q3,
	IQR,
	Min,
	Max,
	Range,
	Percentile,
	SampleRate,
	Sum,
	Count,
	HitCountAbove,
	HitCountBelow,
	PointCount,
	TimeCountAbove,
	TimeCountBelow,
	MeanT,
	Intercept,
	Slope,
	SlopeR2,
	Area,
	MaxGap,
	MinGap,
	Trend,
	Cp,
	Cpk,
	Fft,
	CNT
};

class CMain
{
public:
	json ReadJsonFile(const char* path);
	int GetWindowsIdx(string input);
	int GetViewMode(string input);
	int GetCouplingType(std::string input);
	int GetVoltageRange(std::string input);

	int GetStatIndex(string sName);
	
	bool MappingRecipe(const char* path);
	vector<_ST_WINDOW> rcpWindows;

	bool MappingSetting(const char* path);
	_ST_HARDWARE _hardware;
	_ST_SOFTWARE _software;
};
