#include "main.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cstdint>
#include <stdlib.h>
#include <stdio.h>
#include <json.hpp>
#include <chrono>

using json = nlohmann::json;
using namespace std;
using namespace Eigen;
//constexpr char PATH[] = "D:\\statistics\\samples\\115052_624.bin"; //wgs두개임
constexpr char PATH[] = "D:\\statistics\\samples\\114557_869.bin"; //anno 없는데 wgs있음
//constexpr char PATH[] = "D:\\statistics\\samples\\031713_164.bin";



double simpleRMS(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();
	int nCnt = arr.size();
	double const rms = numext::sqrt((arr.squaredNorm()) / nCnt);
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;

	cout << "-----[simpleRMS]------" << endl;
	cout << "value : " << rms << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;

	return rms;
}

double simpleMEAN(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();
	const double mean = arr.mean();
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;

	cout << "-----[simpleMEAN]------" << endl;
	cout << "value : " << mean << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return mean;
}

double simpleMEANH(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();
	int nCnt = arr.size();
	double sum = 0;
	for (int i = 0; i < nCnt; ++i) {
		if (arr[i] == 0) {
			arr[i] = numeric_limits<double>::infinity();
		}
		sum += 1.0 / arr[i];
	}
	double meanH =  nCnt / sum;
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleMEANH]------" << endl;
	cout << "value : " << meanH << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return meanH;
}

double simpleMEANG(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();
	auto meanG = pow(arr.prod(), 1.0 / arr.size());
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleMEANG]------" << endl;
	cout << "value : " << meanG << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return meanG;
}

double simpleSTDEV(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	double stdev = numext::sqrt((arr.array() - arr.mean()).square().sum() / (arr.size() -1));

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleSTDEV]------" << endl;
	cout << "value : " << stdev << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return stdev;
}

double simpleSKEWNESS(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	auto stdev = simpleSTDEV(arr);
	double skewness = ((arr.array() - arr.mean()).pow(3)).sum() / ((arr.size() - 1) * pow(stdev,3));

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleSKEWNESS]------" << endl;
	cout << "value : " << skewness << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return skewness;
}

double simpleKURTOSIS(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	double std_dev = simpleSTDEV(arr);
	double sum = ((arr.array() - arr.mean()).pow(4)).sum();
	double kurtosis = sum / (arr.size() * pow(std_dev,4)) - 3;

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleKURTOSIS]------" << endl;
	cout << "value : " << kurtosis << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return kurtosis;
}

double* simpleMODES(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	unordered_map<double, int> freqMap;
	for (size_t i = 0; i < arr.size(); i++)
	{
		freqMap[arr(i)]++;
	}
	int maxFreq = 0;
	for (auto it = freqMap.begin(); it != freqMap.end(); ++it) {
		if (it->second > maxFreq) {
			maxFreq = it->second;
		}
	}

	vector<double> result;
	for (auto it = freqMap.begin(); it != freqMap.end(); ++it) {
		if (it->second == maxFreq) {
			result.resize((result.size()) + 1);
			result[result.size() - 1] = it->first;
		}
	}

	double* toPtr = result.data();
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleMODES]------" << endl;
	for (size_t i = 0; i < result.size(); i++)
	{
		cout << "value[" <<i << "] : " << result[i] << endl;
	}
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return toPtr;
}

double simpleMODE(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	sort(arr.begin(), arr.end());
	double mode = 0, max_freq = 0, freq = 1;
	for (size_t i = 1; i < arr.size(); i++) {
		if (arr[i] == arr[i - 1]) {
			freq++;
		}
		else {
			if (freq > max_freq) {
				max_freq = freq;
				mode = arr[i - 1];
			}
			freq = 1;
		}
	}
	if (freq > max_freq) {
		mode = arr[arr.size() - 1];
	}

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleMODE]------" << endl;
	cout << "value : " << mode << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return mode;
}

double simpleMEDIAN(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	int nCnt = arr.size();
	double medianValue = 0;
	if (nCnt % 2 == 0)
	{
		int medianIdx1 = (nCnt / 2)-1;
		nth_element(arr.data(), arr.data() + medianIdx1, arr.data() + nCnt);
		double medianValue1 = arr(medianIdx1);

		int medianIdx2 = nCnt / 2;
		nth_element(arr.data(), arr.data() + medianIdx2, arr.data() + nCnt);
		double medianValue2 = arr(medianIdx2);

		medianValue = (medianValue1 + medianValue2) / 2.0;

	}
	else
	{
		int medianIdx = nCnt / 2;

		nth_element(arr.data(), arr.data() + medianIdx, arr.data() + nCnt);
		medianValue = arr(medianIdx);
	}

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleMEDIAN]------" << endl;
	cout << "value : " << medianValue << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return medianValue;
}

double simpleQ1(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	VectorXd lower_half = arr.segment(0, arr.size() / 2);
	sort(lower_half.data(), lower_half.data() + lower_half.size());
	double q1 = 0;
	if (lower_half.size() % 2 == 0)
	{
		q1 = (lower_half(lower_half.size() / 2 - 1) + lower_half(lower_half.size() / 2)) / 2.0;
	}
	else
	{
		q1 = lower_half(lower_half.size() / 2);
	}

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleQ1]------" << endl;
	cout << "value : " << q1 << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return q1;
}

double simpleQ3(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	double q3;
	Eigen::VectorXd upper_half = arr.segment((arr.size() + 1) / 2, arr.size() / 2); 
	sort(upper_half.data(), upper_half.data() + upper_half.size()); 
	if (upper_half.size() % 2 == 0) 
	{
		q3 = (upper_half(upper_half.size() / 2 - 1) + upper_half(upper_half.size() / 2)) / 2.0;
	}
	else
	{
		q3 = upper_half(upper_half.size() / 2);
	}

	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleQ3]------" << endl;
	cout << "value : " << q3 << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return q3;
}

double simpleIQR(VectorXd arr)
{
	auto start = chrono::high_resolution_clock::now();

	double q1 = simpleQ1(arr);
	double q3 = simpleQ3(arr);
	double iqr = q3 - q1;

	auto end = chrono::high_resolution_clock::now();

	chrono::duration<double, milli> duration = end - start;
	cout << "-----[simpleIQR]------" << endl;
	cout << "value : " << iqr << endl;
	cout << "duration : " << duration.count() << "ms" << endl;
	cout << "----------------------" << endl;
	return iqr;
}

vector<vector<double>> rawDatas;
vector<vector<double>> fftDatas;
int crtPosition = 0;

bool readBinFile2()
{
	try
	{
		ifstream file(PATH, ios::binary);
		const char nullStr = '\0';
		if (!file || !file.good()) return false;
		
		//GET FILE LENGTH
		file.seekg(0, ios::end);
		const int fileSize = (int)file.tellg();

		//SET DATA BUFFER
		char* dataBuffer = new char[fileSize];
		memset(dataBuffer, 0x00, fileSize);

		//GET DATA
		file.seekg(0, ios::beg);
		file.read(dataBuffer, fileSize);

		//PARSING JSON COUNT
		int* jsonCntBin = new int;
		memcpy(jsonCntBin, dataBuffer, sizeof(int));
		crtPosition += sizeof(int);
		const int jsonCntDec = static_cast<int>(*jsonCntBin);
		//cout << "jsonCntDec : " << jsonCntDec << endl;

		//PARSING JSON DATA
		char* jsonData = new char[jsonCntDec + 1];
		memcpy(jsonData, dataBuffer + crtPosition, jsonCntDec);
		crtPosition += jsonCntDec;
		jsonData[jsonCntDec] = nullStr;
		//cout << jsonData << endl;

		//CONVERT TO JSON
		//json파일 안에서 enable 개수 => src cnt, smapleRate => 데이터 수
		json j = json::parse(jsonData);
		
		int chSize = j["ch"].size();
		//cout << "chSize : " << chSize << endl;
		
		//PUT ENABLE SOURCE INTO JSON ARRAY
		json daqArr = json::array();
		for (int i = 0; i < chSize; i++)
		{
			if (j["ch"][i]["enabled"] == 1)
			{
				json daq;
				daq["coupling"] = j["ch"][i]["coupling"];
				daq["sampleRate"] = j["ch"][i]["sampleRate"];
				daq["scale"] = j["ch"][i]["scale"];
				daqArr.push_back(daq);
			}
		}

		//PARSING ANNOTATION
		char* wgs = new char[4];
		int wgsByte = 3;
		memcpy(wgs, dataBuffer + fileSize - wgsByte, wgsByte);
		//int annoLenByte, annoCnt, annoTotalByte;
		int* annoLenByte = new int;
		int annoCnt = 0;
		wgs[3] = nullStr;
		if (strcmp(wgs, "wgs") != 0)
		{
			cout << "no annotation" << endl;
		}
		else
		{
			cout << "there are annotations" << endl;
			memcpy(annoLenByte, dataBuffer + fileSize - wgsByte - 4, 4);
			annoCnt = static_cast<int>(*annoLenByte);
		}


		const int totalDataByte = fileSize - sizeof(int) - jsonCntDec - wgsByte - sizeof(int) - (annoCnt *8);
		int oneCycle = sizeof(double);
		int oneCycleLength = 1;
		const int fftRate = 2048;
 		for (size_t i = 0; i < daqArr.size(); i++)
		{
			int rate = (daqArr[i]["sampleRate"]);
			//oneCycle += (rate * sizeof(double)) + fftRate * sizeof(double);
			oneCycleLength += rate + fftRate;
			oneCycle += sizeof(double) * (rate + fftRate);
		}
		oneCycleLength += fftRate;
		oneCycle += fftRate * sizeof(double); //공통 X Axis

		const int cycleCnt = totalDataByte / oneCycle;


		//PARSING CHART DATA

		const int nSrcCnt = daqArr.size();
		rawDatas.resize(nSrcCnt);
		fftDatas.resize(nSrcCnt);

		crtPosition += 8; //TIME => 8BYTE
		for (int i = 0; i < nSrcCnt; i++)
		{
			int dataSize = (int)daqArr[i]["sampleRate"];
			rawDatas[i].resize(dataSize);
			memcpy(&rawDatas[i][0], dataBuffer + crtPosition, dataSize * sizeof(double));
			crtPosition += dataSize * 8;
		}


		for (int i = 0; i < nSrcCnt; i++)
		{
			int fftSize = 2048;
			fftDatas[i].resize(fftSize);
			memcpy(&fftDatas[i][0], dataBuffer + crtPosition, fftSize * sizeof(double));
			crtPosition += fftSize * sizeof(double);
		}
		
		file.close();
		//DELETE MEMORY
		delete[] dataBuffer;
		delete jsonCntBin;
		delete[] jsonData;
		delete[] wgs;
		dataBuffer = NULL;
		jsonCntBin = NULL;
		jsonData = NULL;
		wgs = NULL;
		return true;
	}
	catch (const exception& e)
	{
		cout << "ERROR : " << e.what() << endl;
		return false;
	}
}


bool toCsv()
{
	try
	{
		auto start = chrono::high_resolution_clock::now();
		auto now = chrono::system_clock::now();
		time_t time = chrono::system_clock::to_time_t(now);
		stringstream ss;
		tm tm;
		gmtime_s(&tm, &time);
		ss << "test_" << put_time(&tm, "%Y%m%d_%H%M%S") << ".csv";
		string fileName = ss.str();
		ofstream dataFile(fileName, ios_base::app);

		for (size_t j = 0; j < rawDatas[0].size(); j++)
		{
			dataFile << rawDatas[0][j] << "\n";
		}

		dataFile.close();
		return true;
	}
	catch (const exception&)
	{
		return false;
	}
}

double* results;

void getResults()
{
	int cycleCnt = 10;
	int crtPosition = 0;
	int oneCycle = 0;
	results = new double[cycleCnt];
	for (int i = 0; i < cycleCnt; i++)
	{
		//getData() 해서 가져온데이터를
		crtPosition += oneCycle;
		//results[i] = //getData해서 가져온 데이터
	}
}

struct _ST_FILE_INFO
{
	int nFileSize = 0;
	int nSrcCnt = 0;
	int nAnnoCnt = 0;
	int nCycleCnt = 0;
	//int nCycleByte = 0;
	int nFftSize = 0;
	int nJsonCnt = 0;
	int nCrtPosition = 0;
	int nRawPos = 0;
	int nFftPos = 0;
	json daqArray = json::array();
	char* dataBuffer = nullptr;
};

void openFile(const char* path, _ST_FILE_INFO* fileInfo)
{
	try
	{
		//Open File
		std::ifstream file(path, std::ios::binary);
		if (!file || !file.good())  throw std::runtime_error("Failed to open file.");

		//Get File Size
		file.seekg(0, std::ios::end);
		fileInfo->nFileSize = (int)file.tellg();

		//Set Data Buffer
		fileInfo->dataBuffer = new char[fileInfo->nFileSize];
		memset(fileInfo->dataBuffer, 0x00, fileInfo->nFileSize);
		file.seekg(0, std::ios::beg);
		file.read(fileInfo->dataBuffer, fileInfo->nFileSize);
	}
	catch (const std::exception& e)
	{
		throw std::runtime_error("Failed to open file.");
	}
}
void getJsonData(_ST_FILE_INFO* fileInfo)
{
	try
	{
		int* jsonCntBin = new int;
		memcpy(jsonCntBin, fileInfo->dataBuffer, sizeof(int));
		fileInfo->nCrtPosition += sizeof(int);
		fileInfo->nJsonCnt = static_cast<int>(*jsonCntBin);

		char* jsonData = new char[fileInfo->nJsonCnt + 1];
		memcpy(jsonData, fileInfo->dataBuffer + fileInfo->nCrtPosition, fileInfo->nJsonCnt);
		fileInfo->nCrtPosition += fileInfo->nJsonCnt;
		//추가 : 여기서 rawCrtPosition = nCrtPostion
		fileInfo->nRawPos = fileInfo->nCrtPosition;

		jsonData[fileInfo->nJsonCnt] = '\0';
		cout << "JSON DATA : " << jsonData << endl;

		const json parsedJson = json::parse(jsonData);

		for (size_t i = 0; i < parsedJson["ch"].size(); i++)
		{
			if (parsedJson["ch"][i]["enabled"] == 1)
			{
				json daq;
				daq["coupling"] = parsedJson["ch"][i]["coupling"];
				daq["sampleRate"] = parsedJson["ch"][i]["sampleRate"];
				daq["scale"] = parsedJson["ch"][i]["scale"];
				fileInfo->daqArray.push_back(daq);
			}
		}

		fileInfo->nSrcCnt = fileInfo->daqArray.size();

		delete jsonCntBin;
		delete[] jsonData;
		jsonCntBin = NULL;
		jsonData = NULL;
	}
	catch (const std::exception&)
	{
		throw std::runtime_error("Failed to get Json data.");
	}
}

const int WGS_BYTE = 3;
void getAnnotation(_ST_FILE_INFO* fileInfo)
{
	try
	{
		char* wgs = new char[WGS_BYTE + 1];
		memcpy(wgs, fileInfo->dataBuffer + fileInfo->nFileSize - WGS_BYTE, WGS_BYTE);
		int* annoLenByte = new int;
		wgs[WGS_BYTE] = '\0';
		if (strcmp(wgs, "wgs") != 0)
		{
			std::cout << "no annotation" << std::endl;
		}
		else
		{
			std::cout << "there are annotations" << std::endl;
			memcpy(annoLenByte, fileInfo->dataBuffer + fileInfo->nFileSize - WGS_BYTE - sizeof(int), sizeof(int));
			fileInfo->nAnnoCnt = static_cast<int>(*annoLenByte);
		}

		delete[] wgs;
		delete annoLenByte;
		wgs = NULL;
		annoLenByte = NULL;
	}
	catch (const std::exception&)
	{
		throw std::runtime_error("Failed to get annotation.");
	}
}

void getCycleCnt(_ST_FILE_INFO* fileInfo)
{
	try
	{
		const int totalDataByte = fileInfo->nFileSize - sizeof(int) - fileInfo->nJsonCnt - WGS_BYTE - sizeof(int) - (fileInfo->nAnnoCnt * 8);
		int oneCycle = sizeof(double); //Time
		//int oneCycleLength = 1; //Time
		const int fftRate = 2048;
		for (size_t i = 0; i < fileInfo->daqArray.size(); i++)
		{
			int rate = (fileInfo->daqArray[i]["sampleRate"]);
			//oneCycleLength += rate + fftRate;
			oneCycle += sizeof(double) * (rate + fftRate);
			fileInfo->nFftSize += sizeof(double) * fftRate;
		}
		fileInfo->nFftSize += sizeof(double) * fftRate;
		//oneCycleLength += fftRate; //공통 X Axis
		oneCycle += fftRate * sizeof(double); //공통 X Axis
		//fileInfo->nCycleByte = oneCycle;
		fileInfo->nCycleCnt = totalDataByte / oneCycle; //int가 아닐 때 에러처리 해야됨
		cout << "nCycleCnt : " << fileInfo->nCycleCnt << " / " << "totalDataByte : " << totalDataByte << " / " << "oneCycle : " << oneCycle << endl;
	}
	catch (const std::exception&)
	{
		throw std::runtime_error("Failed to get cycle count.");

	}
}

std::vector<std::vector<double>> getSrcDataPerCycle(_ST_FILE_INFO* fileInfo)
{
	std::vector<std::vector<double>>  rawData(fileInfo->nSrcCnt);
	fileInfo->nRawPos += sizeof(double); //Time
	//fileInfo->nCrtPosition += sizeof(double);//Time 
	for (int i = 0; i < fileInfo->nSrcCnt; i++)
	{
		cout << "fileInfo->nRawPos : " << fileInfo->nRawPos << endl;

		int dataSize = (int)fileInfo->daqArray[i]["sampleRate"];
		rawData[i].resize(dataSize);
		//memcpy(&rawData[i][0], fileInfo->dataBuffer + fileInfo->nCrtPosition, dataSize * sizeof(double));
		memcpy(&rawData[i][0], fileInfo->dataBuffer + fileInfo->nRawPos, dataSize * sizeof(double));
		//fileInfo->nCrtPosition += dataSize * sizeof(double);
		fileInfo->nRawPos += dataSize * sizeof(double);
		//cout << "nRawPose : " << fileInfo->nRawPos << " / " << "nCrtPos : " << fileInfo->nCrtPosition << endl;
	}
	fileInfo->nFftPos = fileInfo->nRawPos;
	//fileInfo->nRawPos = fileInfo->nCrtPosition;
	return rawData;
}

std::vector<std::vector<double>> getFftDataPerCycle(_ST_FILE_INFO* fileInfo)
{
	std::vector<std::vector<double>>  fftData(fileInfo->nSrcCnt);
	for (int i = 0; i < fileInfo->nSrcCnt; i++)
	{
		int fftSize = 2048;
		fftData[i].resize(fftSize);
		//memcpy(&fftData[i][0], fileInfo->dataBuffer + fileInfo->nCrtPosition, fftSize * sizeof(double));
		memcpy(&fftData[i][0], fileInfo->dataBuffer + fileInfo->nFftPos, fftSize * sizeof(double));
		//fileInfo->nCrtPosition += fftSize * sizeof(double);
		fileInfo->nRawPos += fftSize * sizeof(double);
	}
	return fftData;
}


int main(){
	//VectorXd arr(12);
	//arr <<  65, 70.5 ,70.5, 75, 80,80, 82, 85, 90.5,90.5, 95, 100;

	//readBinFile2();
	//VectorXd v = Map<VectorXd, Unaligned>(rawDatas[0].data(), rawDatas[0].size());
	//simpleRMS(v);
	//simpleMEAN(v);
	//simpleMEANH(v);
	//simpleMEANG(v);
	//simpleMODE(v);
	//simpleMEDIAN(v);
	//simpleKURTOSIS(v);
	//simpleSTDEV(v);
	//simpleSKEWNESS(v);
	//simpleIQR(v);
	//simpleMODES(v);
	//cout << "최대 ? " << v.maxCoeff() << endl;


	_ST_FILE_INFO* fileInfo = new _ST_FILE_INFO();
	openFile(PATH, fileInfo);
	getJsonData(fileInfo);
	getAnnotation(fileInfo);
	getCycleCnt(fileInfo);
	//getSrcDataPerCycle(fileInfo)[0].size();
		 
	for (int i = 0; i < fileInfo->nCycleCnt; i++)
	{
		std::vector<std::vector<double>> data = getSrcDataPerCycle(fileInfo);
		cout << "첫번째 데이터 " << i << " : " << data[0][0] << endl;
		fileInfo->nRawPos += fileInfo->nFftSize;
		VectorXd tr = Map<VectorXd, Unaligned>(data[0].data(), data[0].size());
		simpleRMS(tr);
		//일단 여기까지한게 X rawData까지만 한거
	}

	delete fileInfo;
}

