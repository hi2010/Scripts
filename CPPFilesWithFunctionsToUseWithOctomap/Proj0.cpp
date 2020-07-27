//============================================================================
// Name        : Proj0.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Testprogram for octreelib and other things
// OS		   : Linux / Ubuntu 18.10
// Dependencies: octomap, octovis, binvox(in path loc (/usr/local/bin/, o.e.)) (, linux), C++17
//
// deps: boost, gnuplot iostram, octomap
//		+ alle dependencies der genannten packages
//============================================================================

#include "Proj0.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "octomap/octomap.h"
#include "octomap/OcTree.h"
//#include "QtCore/QVariant"
#include "octomap/ColorOcTree.h"
#include "octomap/math/Vector3.h"
#include "octomap/OcTreeNode.h"
#include "octomap/OcTreeKey.h"
#include "octomap/Pointcloud.h"
#include "octomap/OcTreeIterator.hxx"

#include "gnuplot-iostream.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "boost/multi_array.hpp"

//#include "proto/caffe2.pb.h"
//#include "plaidml/plaidml++.h"
//#include "caffe2/core/init.h"
//#include "caffe2/core/tensor.h"
//#include <torch/all.h>
//#include "plaidml/plaidml++.h"
// may probably be excluded.
#include "keras_model.h"

#include <chrono>

#include "omp.h"

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <filesystem>
namespace fs = std::filesystem;


#include <limits>


using namespace std;
using namespace octomap;

//TODO tackle size problem: size of octree, at least from octovis is wrong / too large vll faktor 1000? (mm als m)


//----------------------------------------------------------------------------------------
// definitions / constants
//----------------------------------------------------------------------------------------

// if 1: select file dialog opens
// if 0: default file is used
static const int SELECT_FILE = 1;

// if 1: binvox is used to voxelize the selected file
// if 0: the object is not going to be voxelized
static const int VOXELIZE_FILE = 1;

//
static const int CREATE_BT_FROM_BINVOX = 1;

// 1: opens bt file in octovis
// 0: does not open octovis
static const int OPEN_FILE_IN_OCTOVIS = 1;

// number of voxel per dimension for voxeliser
static const int VOXEL_RESOLUTION = 256;

static const int VOXELIZE_EXACT = 1;

static const int CONVERT_BT_TO_OT = 1;

// if 1: time for raycasting in areaToHeightmap is calculated and output
static const int BENCH_HEIGHTMAP_GEN = 0;

static const int BENCH_DEPTHMAP_GEN = 0;

// if 1: time for state generation for one part in testWorkflow is calculated and output
static const int BENCH_MACHINE_STATE_GEN = 1;

// file extension for statefiles
static const string STATEFILE_EXTENSION = ".state";
// after eacht value of the state / value separator in matrix line
static const string SEPARATOR = ",";
// at the end of each state line / matrix line
static const string LINE_SEPARATOR = "\n";
// at the end of each state reperesentation
static const string STATE_SEPARATOR = "\n\n";

static const int PRINT_FOUND_FILES = 1;

//----------------------------------------------------------------------------------------
// typedefs
//----------------------------------------------------------------------------------------

// boost multiarray 2D used for heightmap
typedef boost::multi_array<float, 2> fMultiArr2D;
// index for boost multiarray 2D
//typedef fMultiArr2D::index fMA2DIndex;

//----------------------------------------------------------------------------------------
// program start
//----------------------------------------------------------------------------------------

int main() {
//	modelToOctovis();
//	for(int i = 0; i <= 99; i++){
//		string filePath = "/home/louis/Documents/Models/data/stlfiles/testFile"+to_string(i)+".stl";
//		stlToOt(filePath);
//	}
//	string filePath = "/home/louis/Documents/Models/data/stlfiles/testFile2.stl";
//	stlToOt(filePath);
//	openOtFileWithOctomap();

//	testWorkflow3();
//	testWorkflow4();
//	testWorkflow5();
//TODO
	testWorkflow7();
//	printHello();

//	testWorkflow7();

//	testWorkflow8();

//	string result = exec("binvox2bt -h");
//	cout << result << endl;
//	int sl1 = 256;
//	int sl2 = 256;
//	float* maxVals = new float[sl1][sl2];
//	typedef	boost::multi_array<float, 2> arr2D;
//	typedef arr2D::index index;
//	arr2D maxVals(boost::extents[sl1][sl2]);
//	multiArray2DToVector(maxVals);
	return 0;
}

//----------------------------------------------------------------------------------------
// workflow methods
//----------------------------------------------------------------------------------------

/**
 * This method is used for testing
 */
void testWorkflow(void){
	// test prog is working at all
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	// test function calling
	printHello();

	// test octomap lib / header is working
	//octomap::
	point3d offset(0.0, 0.0, 0.0);
	point3d offs2(1, 1, 1);
	cout << "created point 3d using octomap" << endl;

	// test file open dialog
	// only if file select is active, else use default file
	string fileToOpen = "";
	// got to add a method for file extension extraction
	string fileExtension = ".stl";

	if (SELECT_FILE){
		fileToOpen =  openFileRawName();
	} else {
		// return of sprintf is not important, as there is no variable
		// and if the path is bad the programmer -> me probbably fucked up :D
		fileToOpen = "/home/louis/Documents/Models/CSG_1";
		cout << fileToOpen << endl;
	}

	// define stringstream, used to add file extensions to the fileToOpen raw path
	stringstream ss;

	if (VOXELIZE_FILE){
		ss << fileToOpen << fileExtension;
		voxelizeFile(ss.str() , VOXEL_RESOLUTION, 0);
		ss.str(string()); // clear ss
	}

	if (CREATE_BT_FROM_BINVOX){
		ss << fileToOpen << ".binvox";
		createBtFromBinvoxfile(ss.str());
		ss.str(string());
	}

	if (OPEN_FILE_IN_OCTOVIS){
		ss << fileToOpen << ".binvox.bt";
		openFileInOctovis(ss.str());
		ss.str(string());
	}

}

void testWorkflow2(){
		OcTree btOt("/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt");
		Pointcloud cylinderTool = createCircularPointCloud(1,btOt.getResolution());
		OcTree machinedOt(btOt.getResolution());
	//	insertTool(&machinedOt,cylinderTool,pose6d(0,0,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(1,0,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(2,0,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,0,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,1,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,2,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,3,0,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,3,1,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,3,2,0,0,0));
	//	insertTool(&machinedOt,cylinderTool,pose6d(3,3,3,0,0,0));

		// create infinite loop
		insertToolInLine(&machinedOt,cylinderTool,point3d(0,0,0),point3d(3,0,0));
		insertToolInLine(&machinedOt,cylinderTool,point3d(3,0,0),point3d(3,3,0));
		insertToolInLine(&machinedOt,cylinderTool,point3d(3,3,0),point3d(3,3,3));
		insertToolInLine(&machinedOt,cylinderTool,point3d(3,3,3),point3d(0,0,0));

		saveOcTree(&machinedOt, true);
}

void testWorkflow3(void){

	OcTree btOt("/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt");
	Pointcloud cylinderTool = createCircularPointCloud(1,btOt.getResolution());
	OcTree machinedOt(btOt.getResolution());

	double btOtMin[3];
	double btOtMax[3];
	btOt.getMetricMin(btOtMin[0],btOtMin[1],btOtMin[2]);
	btOt.getMetricMax(btOtMax[0],btOtMax[1],btOtMax[2]);

	// these two files should look the same, but dont
	// still looks way better. investigate if gcode scaling is screwed up
//	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGCodeFiles/scScaledtestFile0.ngc");
//	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGcode2/scScaledtestFile0.ngc");
//	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scledGCode001mm/scScaledtestFile0.ngc");
	// scaled gcode und somit py2 enthält fehler
//	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGCodes3_001mm/scScaledtestFile0.ngc");
	// rewo3 of py2 looks good again hopefully this is good one
//	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGCodes3_001mm/scScaledtestFile0.ngc");
	// rewo3 with 1mm
	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGCode1mmWithPy2/scScaledtestFile0.ngc");

	ofstream stateFile = ofstream("/home/louis/Documents/Models/testStatesFile0_1.state");

	point3d origin;
	point3d destination;

	if (file.peek() != EOF){
		origin = getNextPointFromScaledGCode(file);
		insertTool(&machinedOt, cylinderTool, pose6d(origin.x(), origin.y(), origin.z(), 0, 0, 0) );
	}

	float firstZ = origin.z();

	chrono::microseconds usBeforeStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	int stateNum = 0;
	while (file.peek() != EOF){
		destination = getNextPointFromScaledGCode(file);

		double min[3];
		double max[3];
		machinedOt.getMetricMin(min[0],min[1],min[2]);
		machinedOt.getMetricMax(max[0],max[1],max[2]);
		cout << "mins: " << min[0] << " : " << min[1] << " : " << min[2] << endl;
		cout << "maxs: " << max[0] << " : " << max[1] << " : " << max[2] << endl;


		vector<vector<float>> tempDm = areaToDepthmapWithBBXOptimised(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
											,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});


		write2DVectorToFile(tempDm, stateFile);

		// round up looks better
		insertToolInLineRoundUp(&machinedOt, cylinderTool, origin, destination);
//		insertToolInLine(&machinedOt, cylinderTool, origin, destination);

		origin.x() = destination.x();
		origin.y() = destination.y();
		origin.z() = destination.z();

		stateNum ++;
		if (stateNum % 100 == 0){
			cout << "state: " << stateNum << endl;
			cout << "dm Dim: " << tempDm.size() << " : " << tempDm[0].size() << endl;
		}
	}

	// write last state
	double min[3];
	double max[3];
	machinedOt.getMetricMin(min[0],min[1],min[2]);
	machinedOt.getMetricMax(max[0],max[1],max[2]);
	vector<vector<float>> tempDm = areaToDepthmapWithBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
												,0 ,(float) min[0], (float) min[1], (float) max[0], (float) max[1]);
	write2DVectorToFile(tempDm, stateFile);


	cout << "all states created. Count: " << stateNum << endl;

	if (BENCH_MACHINE_STATE_GEN){
			chrono::microseconds usAfterStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
			cout << "Machine state gen. took: " << (usAfterStateGen.count()-usBeforeStateGen.count()) << " us or: " << (usAfterStateGen.count()-usBeforeStateGen.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

	file.close();

	saveOcTree(&machinedOt, true);
}


void testWorkflow4(void){
	OcTree btOt("/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt");
	double btOtMin[3];
	double btOtMax[3];
	btOt.getMetricMin(btOtMin[0],btOtMin[1],btOtMin[2]);
	btOt.getMetricMax(btOtMax[0],btOtMax[1],btOtMax[2]);

	float xStart = (float)btOtMin[0];
	float yStart = (float)btOtMin[1];
	float resolution = btOt.getResolution();
	int xSize = ceil((btOtMax[0]-btOtMin[0])/resolution);
	int ySize = ceil((btOtMax[1]-btOtMin[1])/resolution);

//	auto myMap = areaToHeightmapWithLeafIterator(&btOt, 0, 0,
//			xSize, ySize, resolution, 0);
//	printOtSpecs(&btOt);
//	auto myMap2 = areaToHeightmapWithLeafIteratorAndBBX(&btOt,17,17,40,40,resolution,0,12,bbx2D{0,0,200,200});
	auto myMap2 = areaToDepthmapWithLeafIteratorAndBBX(&btOt,17,17,40,40,resolution,0,12,bbx2D{0,0,200,200});
	cout << "tw4: Ot num nodes: " << btOt.getNumLeafNodes() << endl;

//	plot2DVectorWithGnuPlot(myMap);
	plot2DVectorWithGnuPlot(myMap2);
	ofstream file = ofstream("/home/louis/Documents/Models/leafDmGen");
	write2DVectorToFile(myMap2, file);

}


void testWorkflow5(void){

	//TODO implement new insertToolInLine
	//TODO use files list to create trainings data
	auto filesList = retrieveFilesListFromPath("/home/louis/Documents/Models/",".state");

//	return;

	// open targetfile as octree, needed to retrieve parameters and generate states
	OcTree btOt("/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt");
	Pointcloud cylinderTool = createCircularPointCloud(1,btOt.getResolution());
	OcTree machinedOt(btOt.getResolution());

	double btOtMin[3];
	double btOtMax[3];
	btOt.getMetricMin(btOtMin[0],btOtMin[1],btOtMin[2]);
	btOt.getMetricMax(btOtMax[0],btOtMax[1],btOtMax[2]);

	//------------------------------------------------------------
	// create files
	//------------------------------------------------------------

	// rewo3 with 1mm
	ifstream file = readGCodeFile("/home/louis/Documents/Models/data/scaledGCode1mmWithPy2/scScaledtestFile0.ngc");

	// file path for states of the machined part
	string stateFilePath = "/home/louis/Documents/Models/testStatesFile0_isStates.state";
	// if file does not end with statefile extension add statefile extension
	if (!(stateFilePath.substr(stateFilePath.find_last_of("."))==STATEFILE_EXTENSION)) stateFilePath += STATEFILE_EXTENSION;
	ofstream stateFile = ofstream(stateFilePath);

	// file path for states of the target file
	string stateFilePathTarget = "/home/louis/Documents/Models/testStatesFile0_targetStates.state";
	// if file does not end with statefile extension add statefile extension
	if (!(stateFilePathTarget.substr(stateFilePathTarget.find_last_of("."))==STATEFILE_EXTENSION)) stateFilePathTarget += STATEFILE_EXTENSION;
	ofstream stateFileTarget = ofstream(stateFilePathTarget);

	//------------------------------------------------------------

	//------------------------------------------------------------
	// create training files
	//------------------------------------------------------------



	point3d origin;
	point3d destination;

	if (file.peek() != EOF){
		origin = getNextPointFromScaledGCode(file);
		insertTool(&machinedOt, cylinderTool, pose6d(origin.x(), origin.y(), origin.z(), 0, 0, 0) );
	}

	float firstZ = origin.z();

	chrono::microseconds usBeforeStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

//	vector<vector<float>> btOtdm = areaToDepthmapWithLeafIteratorAndBBX(&btOt, btOtMin[0], btOtMin[0], 256, 256, .1
//													,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});
	float btOtDmResol = .1;
	int btOtDmXLen = ceil((btOtMax[0]-btOtMin[0])/btOtDmResol);
	int btOtDmYLen = ceil((btOtMax[1]-btOtMin[1])/btOtDmResol);
	vector<vector<float>> btOtdm = areaToHeightmapWithLeafIterator(&btOt,btOtMin[0],btOtMin[1],btOtDmXLen,btOtDmYLen,btOtDmResol,0);
//	plot2DVectorWithGnuPlot(btOtdm);

	int stateNum = 0;
	while (file.peek() != EOF){
		destination = getNextPointFromScaledGCode(file);

		double min[3];
		double max[3];
		machinedOt.getMetricMin(min[0],min[1],min[2]);
		machinedOt.getMetricMax(max[0],max[1],max[2]);


		// only create new state if points are different
		if (!( (origin.x()==destination.x()) &&
				(origin.y()==destination.y()) &&
				(origin.z()==destination.z()) )){
			if (euclideanDistance(origin,destination) < machinedOt.getResolution()){
				cout << "\tdistance is smaller than resol, distance: " << euclideanDistance(origin,destination)
						<< " stateNum is: " << stateNum << endl;
			}
			// -2 cause window has size 4
			vector<vector<float>> tempDm = areaToDepthmapWithLeafIteratorAndBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
														,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});
			vector<vector<float>> targetDm = depthMapToPartlyDepthMap(btOtdm,.1,btOtMin[0],btOtMin[1],origin.x()-2,origin.y()-2,40,40,0);

			stateFile << "clPoint: " << origin.x() << "," << origin.y() << "," << origin.z() << endl;
			stateFileTarget << "clPoint: " << origin.x() << "," << origin.y() << "," << origin.z() << endl;
			write2DVectorToFile(tempDm, stateFile);
			write2DVectorToFile(targetDm, stateFileTarget);
			stateNum ++;
			if (stateNum % 100 == 0){
				cout << "state: " << stateNum << endl;
				cout << "dm Dim: " << tempDm.size() << " : " << tempDm[0].size() << endl;
			}
		}

		// round up looks better
//		insertToolInLineRoundUp(&machinedOt, cylinderTool, origin, destination);
//		insertToolInLine(&machinedOt, cylinderTool, origin, destination);
		insertToolInLineUsingOctomapComputeRayKeys(&machinedOt, cylinderTool, origin, destination);

		origin.x() = destination.x();
		origin.y() = destination.y();
		origin.z() = destination.z();

	}

	// write last state
	double min[3];
	double max[3];
	machinedOt.getMetricMin(min[0],min[1],min[2]);
	machinedOt.getMetricMax(max[0],max[1],max[2]);
	vector<vector<float>> tempDm = areaToDepthmapWithBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
												,0 ,(float) min[0], (float) min[1], (float) max[0], (float) max[1]);
	write2DVectorToFile(tempDm, stateFile);


	cout << "all states created. Count: " << stateNum << endl;

	if (BENCH_MACHINE_STATE_GEN){
			chrono::microseconds usAfterStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
			cout << "Machine state gen. took: " << (usAfterStateGen.count()-usBeforeStateGen.count()) << " us or: " << (usAfterStateGen.count()-usBeforeStateGen.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

	file.close();

	saveOcTree(&machinedOt, true);
}


/**
 * This Method was used to create trainingdata
 */
void testWorkflow6(void){
	// create fileslist
	vector<string> filesList = retrieveFilesListFromPath("/home/louis/Documents/Models/data/btFiles/",".bt");

	chrono::microseconds usBeforeAllStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

#pragma omp parallel for
	for (int i = 0; i < filesList.size(); i++){
		string file = filesList[i];
//	}
//	for (string file : filesList){
		string filename = getFileNameFromPathWithoutFileExtension(file,2);
		// inputs
		string gCodePath = "/home/louis/Documents/Models/data/gCodeFixedAndScaled/scScaled" + filename + ".ngc";
		string btPath = "/home/louis/Documents/Models/data/btFiles/" + filename + ".binvox.bt";
		// outputs
		string machinedBtPath = "/home/louis/Documents/Models/data/machinedBtFiles4/" + filename;
		string stateFilePath = "/home/louis/Documents/Models/data/isStateFiles4/" + filename + STATEFILE_EXTENSION;
		string targetStateFilePath = "/home/louis/Documents/Models/data/targetStateFiles4/" + filename + STATEFILE_EXTENSION;
		createTrainingdata(gCodePath, btPath, machinedBtPath, stateFilePath, targetStateFilePath);
	}

	if (BENCH_MACHINE_STATE_GEN){
			chrono::microseconds usAfterAllStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
			cout << "\n\nCreated all machie states." << endl;
			cout << "Gen. of all machine states took: " << (usAfterAllStateGen.count()-usBeforeAllStateGen.count()) << " us or: " << (usAfterAllStateGen.count()-usBeforeAllStateGen.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

}


// test predictor network
void testWorkflow7(void){
	string gCodeFilePath = "/home/louis/Documents/eclipse-workspace/Proj0/src/knnScripts/testGCode.gcode";
	ofstream gCodeFile = ofstream(gCodeFilePath);
	gCodeFile << "G91" << "\n" << "G01" << endl;
	string btOtPath = "/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt";
	OcTree btOt = OcTree(btOtPath);
//	string machinedBtPath = "/home/louis/Documents/Models/data/machinedBtFiles3/testFile0.bt";
//	OcTree mOt = OcTree(machinedBtPath);
	double resol = btOt.getResolution();
	OcTree mOt = OcTree(resol);

	double mOtMin[3];
	double mOtMax[3];

	mOt.getMetricMin(mOtMin[0], mOtMin[1], mOtMin[2]);
	mOt.getMetricMax(mOtMax[0], mOtMax[1], mOtMax[2]);
	printOtSpecs(&mOt);

	double btOtMin[3];
	double btOtMax[3];

	btOt.getMetricMin(btOtMin[0], btOtMin[1], btOtMin[2]);
	btOt.getMetricMax(btOtMax[0], btOtMax[1], btOtMax[2]);

	point3d origin;

	float firstZ = btOtMax[2] + 1;

	printOtSpecs(&btOt);
	cout << "firstZ: " << firstZ << endl;

	origin.x() = btOtMin[0];
	origin.y() = btOtMin[1];
	origin.z() = btOtMax[2];


	string pathToVenv = "/home/louis/plaidml-venv/bin/python3.6";
	string currentPath = GetCurrentWorkingDir();
	cout << "working in: " << currentPath << endl;
	string cmd = pathToVenv + " " + currentPath + "/src/knnScripts/genPredictions.py";

	for (int i = 0; i < 400; i++){

		vector<vector<float>> tempDm = areaToDepthmapWithLeafIteratorAndBBX(
				&mOt, origin.x() - 2, origin.y() - 2, 40, 40, .1, 0, firstZ,
				bbx2D { (float) btOtMin[0], (float) btOtMin[1],
						(float) btOtMax[0], (float) btOtMax[1] });
		vector<vector<float>> targetDm = areaToHeightmapWithLeafIteratorAndBBX(
				&btOt, origin.x() - 2, origin.y() - 2, 40, 40, .1, 0, 0, bbx2D {
						(float) btOtMin[0], (float) btOtMin[1],
						(float) btOtMax[0], (float) btOtMax[1] });

		ofstream gCodeFile = ofstream(gCodeFilePath);
		gCodeFile << "X" << origin.x() << " Y" << origin.y() << " Z" << origin.z() << endl;
		gCodeFile.flush();
		gCodeFile.close();

		ofstream isFile = ofstream(
				currentPath + "/src/knnScripts/isState.state");
		ofstream targetFile = ofstream(
				currentPath + "/src/knnScripts/targetState.state");

		isFile << "clPoint: " << origin.x() << "," << origin.y() << ","
				<< origin.z() << endl;
		targetFile << "clPoint: " << origin.x() << "," << origin.y() << ","
				<< origin.z() << endl;

		write2DVectorToFile(tempDm, isFile);
		write2DVectorToFile(targetDm, targetFile);

		isFile.flush();
		targetFile.flush();
		isFile.close();
		targetFile.close();

		plot2DVectorWithGnuPlot(tempDm);
		plot2DVectorWithGnuPlot(targetDm);

		cout << "gen pred." << endl;

		string res = exec(cmd.c_str());
		cout << "pred: " << res << endl;

		// extract point
		string filteredRes = res.substr(res.find("[") + 2,
				(res.find("]")) - (res.find("[") + 2));
		// if begins with space remove that
		cout << "filtered res: " << filteredRes << endl;
		string strPos[3];
		stringstream ss(filteredRes);
		ss >> strPos[0] >> strPos[1] >> strPos[2];

		cout << "x: " << strPos[0] << " y: " << strPos[1] << " z: " << strPos[2]
				<< endl;
		float fPos[3];
		fPos[0] = min(stof(strPos[0]) * 1.25, 1.);
		fPos[1] = min(stof(strPos[1]) * 1.25, 1.);
		fPos[2] = min(stof(strPos[2]) * 1.25, 1.);

		point3d destination(origin.x() + fPos[0], origin.y() + fPos[1],
				origin.z() + fPos[2]);

		cout << origin << ":::::" << destination << endl;
		cout << "poinNum: " << i << endl;

		Pointcloud cylinderCutter = createCircularPointCloud(1, resol);

		insertToolInLineUsingOctomapComputeRayKeys(&mOt, cylinderCutter, origin,
				destination);

		origin.x() = destination.x();
		origin.y() = destination.y();
		origin.z() = destination.z();

		saveOcTree(&mOt, true, "/home/louis/Documents/Models/data/testRuns2/run/testOt"+to_string(i)+".bt");
	}

}


/**
 * This Method was used to create trainingdata
 */
void testWorkflow8(void){
	// create fileslist
	// vector<string> filesList = retrieveFilesListFromPath("/home/louis/Documents/Models/data/btFiles/",".bt");

	chrono::microseconds usBeforeAllStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	string filename = "testFile65";
		// inputs
		string gCodePath = "/home/louis/Documents/Models/data/gCodeFixedAndScaled/scScaled" + filename + ".ngc";
		string btPath = "/home/louis/Documents/Models/data/btFiles/" + filename + ".binvox.bt";
		// outputs
		string machinedBtPath = "/home/louis/Documents/Models/data/machinedBtFiles4/" + filename;
		string stateFilePath = "/home/louis/Documents/Models/data/isStateFiles4/" + filename + "IsState" + STATEFILE_EXTENSION;
		string targetStateFilePath = "/home/louis/Documents/Models/data/targetStateFiles4/" + filename + "TargetState" + STATEFILE_EXTENSION;
		createTrainingdata(gCodePath, btPath, machinedBtPath, stateFilePath, targetStateFilePath);


	if (BENCH_MACHINE_STATE_GEN){
			chrono::microseconds usAfterAllStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
			cout << "\n\nCreated all machie states." << endl;
			cout << "Gen. of all machine states took: " << (usAfterAllStateGen.count()-usBeforeAllStateGen.count()) << " us or: " << (usAfterAllStateGen.count()-usBeforeAllStateGen.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

}


/**
 * This method is used to create a bt file from a 3D-Model with supported format.
 * The created octree is then shown with octovis
 */
void modelToOctovis(void){
	// get filePath
	string fileToOpen = "";
	if (SELECT_FILE){
		fileToOpen = openFile();
	} else {
		fileToOpen = "/home/louis/Documents/Models/CSG_1";
	}

	cout << fileToOpen << endl;

	// get raw path and fileextension
	size_t lastindex = fileToOpen.find_last_of(".");
	cout << lastindex << endl;
	string filePath = fileToOpen.substr(0, lastindex);
	string fileExtension = fileToOpen.substr(lastindex);
	cout << "pat: " << filePath << " ext: " << fileExtension;

	if (VOXELIZE_FILE){
		// create binvox file from model with res 512
		voxelizeFile(fileToOpen, VOXEL_RESOLUTION, VOXELIZE_EXACT);
	}

	if (CREATE_BT_FROM_BINVOX){
		// create bonsai tree (bt) file from binvox file
		stringstream ss;
		ss << filePath << ".binvox";
		createBtFromBinvoxfile(ss.str());
		ss.clear();
	}

	if (CONVERT_BT_TO_OT){
		// convert ot file to bt file
		stringstream ss;
		ss << filePath << ".binvox.bt";
		btToOt(ss.str());
		ss.clear();
	}


	if (OPEN_FILE_IN_OCTOVIS){
		// open file in octovis
		stringstream ss;
		ss << filePath << ".binvox.bt.ot";
		openFileInOctovis(ss.str());
		ss.clear();
	}
}


/**
 * This method is used to create training data from a gcode file, a bt file of the target body
 * and paths for the stateFile of the machined and the target file.
 * inputs:
 * 		<input files>
 * 			gCodePath
 * 		<output files>
 */
void createTrainingdata(string gCodePath, string btPath, string machinedBtPath,
		string stateFilePath, string targetStateFilePath){

	// open targetfile as octree, needed to retrieve parameters and generate states
	OcTree btOt(btPath);
	Pointcloud cylinderTool = createCircularPointCloud(1,btOt.getResolution());
	OcTree machinedOt(btOt.getResolution());

	btOt.setOccupancyThres(.000001);
	machinedOt.setOccupancyThres(.000001);

	double btOtMin[3];
	double btOtMax[3];
	btOt.getMetricMin(btOtMin[0],btOtMin[1],btOtMin[2]);
	btOt.getMetricMax(btOtMax[0],btOtMax[1],btOtMax[2]);

	//------------------------------------------------------------
	// create files
	//------------------------------------------------------------

	// rewo3 with 1mm
	ifstream file = readGCodeFile(gCodePath);

	// if file does not end with statefile extension add statefile extension
	if (!(stateFilePath.substr(stateFilePath.find_last_of("."))==STATEFILE_EXTENSION)) stateFilePath += STATEFILE_EXTENSION;
	ofstream stateFile = ofstream(stateFilePath);

	// if file does not end with statefile extension add statefile extension
	if (!(targetStateFilePath.substr(targetStateFilePath.find_last_of("."))==STATEFILE_EXTENSION)) targetStateFilePath += STATEFILE_EXTENSION;
	ofstream stateFileTarget = ofstream(targetStateFilePath);

	//------------------------------------------------------------

	//------------------------------------------------------------
	// create training files
	//------------------------------------------------------------



	point3d origin;
	point3d destination;

	if (file.peek() != EOF){
		origin = getNextPointFromScaledGCode(file);
		insertTool(&machinedOt, cylinderTool, pose6d(origin.x(), origin.y(), origin.z(), 0, 0, 0) );
	}

	// TODO fix g-code -> first z fixed
	//float firstZ = origin.z();
	float firstZ = btOtMax[2] + 1;

	chrono::microseconds usBeforeStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

//	vector<vector<float>> btOtdm = areaToDepthmapWithLeafIteratorAndBBX(&btOt, btOtMin[0], btOtMin[0], 256, 256, .1
//													,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});

	// create dm for dpethMapToPartlyDepthmap
//	float btOtDmResol = .1;
//	int btOtDmXLen = ceil((btOtMax[0]-btOtMin[0])/btOtDmResol);
//	int btOtDmYLen = ceil((btOtMax[1]-btOtMin[1])/btOtDmResol);
//	vector<vector<float>> btOtdm = areaToHeightmapWithLeafIterator(&btOt,btOtMin[0],btOtMin[1],btOtDmXLen,btOtDmYLen,btOtDmResol,0);
//	plot2DVectorWithGnuPlot(btOtdm);

	int stateNum = 0;
	while (file.peek() != EOF){
		destination = getNextPointFromScaledGCode(file);

		double min[3];
		double max[3];
		machinedOt.getMetricMin(min[0],min[1],min[2]);
		machinedOt.getMetricMax(max[0],max[1],max[2]);


		// only create new state if points are different
		if (!( (origin.x()==destination.x()) &&
				(origin.y()==destination.y()) &&
				(origin.z()==destination.z()) )){
//			if (euclideanDistance(origin,destination) < machinedOt.getResolution()){
//				cout << "\tdistance is smaller than resol, distance: " << euclideanDistance(origin,destination)
//						<< " stateNum is: " << stateNum << endl;
//			}
			// -2 cause window has size 4
			// seems to got n "error" in this method
			// "error" seems to bee taht at the upper limit the machiens ot dimensions are used
			// maybe only last state was real bs, other may be kk.
			// -> in that case the error should be fixed and not matter
			// for other states max .1 difference -> bad, but not that shabby
			// -> thought maybe the reason is, that the method with leaf iterator is "higher resolution"
			// 		-> solution would be to create the dm from which the partly dm is gereated with higher resolution
			vector<vector<float>> tempDm = areaToDepthmapWithLeafIteratorAndBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
														,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});
//			vector<vector<float>> targetDm = depthMapToPartlyDepthMap(btOtdm,.1,btOtMin[0],btOtMin[1],origin.x()-2,origin.y()-2,40,40,0);
			vector<vector<float>> targetDm = areaToHeightmapWithLeafIteratorAndBBX(&btOt, origin.x()-2, origin.y()-2, 40, 40, .1
																	,0, 0 , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});

			stateFile << "clPoint:, " << origin.x() << "," << origin.y() << "," << origin.z() << " ,State:, " << stateNum << " ,,dir:, " << to_string(destination.x()-origin.x()) << "," << to_string(destination.y()-origin.y()) << "," << to_string(destination.z()-origin.z()) << endl;
			stateFileTarget << "clPoint:, " << origin.x() << "," << origin.y() << "," << origin.z() << " ,State:, " << stateNum << " ,,dir:, " << to_string(destination.x()-origin.x()) << "," << to_string(destination.y()-origin.y()) << "," << to_string(destination.z()-origin.z()) << endl;
			write2DVectorToFile(tempDm, stateFile);
			write2DVectorToFile(targetDm, stateFileTarget);
			stateNum ++;
			if (stateNum % 100 == 0){
				cout << "state: " << stateNum << endl;
				cout << "dm Dim: " << tempDm.size() << " : " << tempDm[0].size() << endl;
			}
		}

		// round up looks better
//		insertToolInLineRoundUp(&machinedOt, cylinderTool, origin, destination);
//		insertToolInLine(&machinedOt, cylinderTool, origin, destination);
		insertToolInLineUsingOctomapComputeRayKeys(&machinedOt, cylinderTool, origin, destination);

		origin.x() = destination.x();
		origin.y() = destination.y();
		origin.z() = destination.z();

//		saveOcTree(&machinedOt, 1, machinedBtPath+"state"+to_string(stateNum));

	}

	// write last state
	double min[3];
	double max[3];
	machinedOt.getMetricMin(min[0],min[1],min[2]);
	machinedOt.getMetricMax(max[0],max[1],max[2]);

	stateFile << "clPoint: " << origin.x() << "," << origin.y() << "," << origin.z() << endl;
	stateFileTarget << "clPoint: " << origin.x() << "," << origin.y() << "," << origin.z() << endl;

//	vector<vector<float>> tempDm = areaToDepthmapWithBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
//												,0 ,(float) min[0], (float) min[1], (float) max[0], (float) max[1]);
	vector<vector<float>> tempDm = areaToDepthmapWithLeafIteratorAndBBX(&machinedOt, origin.x()-2, origin.y()-2, 40, 40, .1
															,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});
//	vector<vector<float>> targetDm = depthMapToPartlyDepthMap(btOtdm,.1,btOtMin[0],btOtMin[1],origin.x()-2,origin.y()-2,40,40,0);
	vector<vector<float>> targetDm = areaToHeightmapWithLeafIteratorAndBBX(&btOt, origin.x()-2, origin.y()-2, 40, 40, .1
																,0, firstZ , bbx2D{(float)btOtMin[0], (float)btOtMin[1], (float)btOtMax[0], (float)btOtMax[1]});
	write2DVectorToFile(tempDm, stateFile);
	write2DVectorToFile(targetDm, stateFileTarget);


	cout << "all states created. Count: " << stateNum << endl;

	if (BENCH_MACHINE_STATE_GEN){
			chrono::microseconds usAfterStateGen = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
			cout << "Machine state gen. took: " << (usAfterStateGen.count()-usBeforeStateGen.count()) << " us or: " << (usAfterStateGen.count()-usBeforeStateGen.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

	file.close();

	saveOcTree(&machinedOt, true, machinedBtPath);
}


//TODO comment
string stlToOt(void){
	// get filePath
	string fileToOpen = openFile();


	cout << fileToOpen << endl;

	// get raw path and fileextension
	size_t lastindex = fileToOpen.find_last_of(".");
	cout << lastindex << endl;
	string filePath = fileToOpen.substr(0, lastindex);
	string fileExtension = fileToOpen.substr(lastindex);
	cout << "pat: " << filePath << " ext: " << fileExtension;

	// create binvox file from model with res 512
	voxelizeFile(fileToOpen, VOXEL_RESOLUTION, VOXELIZE_EXACT);

	// create bonsai tree (bt) file from binvox file
	stringstream ss;
	ss << filePath << ".binvox";
	createBtFromBinvoxfile(ss.str());
	ss.clear();
	ss.str("");

	// convert ot file to bt file
	ss << filePath << ".binvox.bt";
	btToOt(ss.str());
	ss.clear();
	ss.str("");

	return "btToOt run";
}

//TODO comment
string stlToOt(string path){
	// get filePath
	string fileToOpen = path;


	cout << fileToOpen << endl;

	// get raw path and fileextension
	size_t lastindex = fileToOpen.find_last_of(".");
	cout << lastindex << endl;
	string filePath = fileToOpen.substr(0, lastindex);
	string fileExtension = fileToOpen.substr(lastindex);
	cout << "pat: " << filePath << " ext: " << fileExtension;

	// create binvox file from model with res 512
	voxelizeFile(fileToOpen, VOXEL_RESOLUTION, VOXELIZE_EXACT);

	// create bonsai tree (bt) file from binvox file
	stringstream ss;
	ss << filePath << ".binvox";
	createBtFromBinvoxfile(ss.str());
	ss.clear();
	ss.str("");

	// convert ot file to bt file
	ss << filePath << ".binvox.bt";
	btToOt(ss.str());
	ss.clear();
	ss.str("");

	return "btToOt run";
}

/**
 * This method creates a bt file from a stl file.
 * The filepath can be choosen with the file open dialogue.
 * input:
 * 		void
 * return:
 * 		int		: 0 if the whole code ran (maybe adding later 0 if succesfull)
 */
string stlToBt(void){
	// open file open dialogue and save filepath as string fileToOpen
	string fileToOpen = openFileRawName();
	// convert the file to a bt file using overloaded method stlToBt
	string stlToBtReturn = stlToBt(fileToOpen);

	return stlToBtReturn;
}

/**
 * This method creates a bt file from a stl file.
 * The filepath can be choosen with the file open dialogue.
 * input:
 * 		string	: path, file path in raw form without file extension
 * 					(this function may be depreacated if i add second string for file extension
 * 					-> smarter -> more formats)
 * return:
 * 		int		: 0 if the whole code ran (maybe adding later 0 if succesfull)
 */
string stlToBt(string path){
	string fileToOpen = path;
	string fileExtension = ".stl";
	// define stringstream, used to add file extensions to the fileToOpen raw path
	stringstream ss;

	// voxelize stl
	ss << "\"" << fileToOpen << "\"" << fileExtension;
	voxelizeFile(ss.str() , 0, 0);
	ss.str(string()); // clear ss

	// create bonsai tree from binvox
	ss << fileToOpen << ".binvox";
	string createBtReturn = createBtFromBinvoxfile(ss.str());
	ss.str(string());

	return createBtReturn;
}

//TODO
/**
 * path complete path with extension
 */
string btToOt(string path){
	stringstream ss;
	ss << "convert_octree " << "\"" << path << "\"";
	cout << "exec: " << ss.str() << endl;
	string convertOctreeReturn = exec(ss.str().c_str());
	cout << convertOctreeReturn << endl;
	ss.clear();
	return convertOctreeReturn;
}

/**
 * This method converts a stl file to a bt file and opens this in octovis.
 * The file is selected with the file open dialogue.
 * input:
 * 		void
 * return:
 * 		int		: 0 if all code ran
 */
int stlToOctovis(void){
	string fileToOpen = openFileRawName();
	string fileExtension = ".stl";
	// define stringstream, used to add file extensions to the fileToOpen raw path
	stringstream ss;

	// create bt from stl
	stlToBt(fileToOpen);

	// open file in binvox
	ss << fileToOpen << ".binvox.bt";
	openFileInOctovis(ss.str());
	ss.str(string());
	return 0;
}

/**
 * This method converts a stl file to a bt file and opens this in octovis.
 * The file is passed with the path vaiable, where the path is raw path,
 * without file extension.
 * input:
 * 		string	: path, file path without extension.
 * return:
 * 		int		: 0 if all code ran
 */
int stlToOctovis(string path){
	string fileToOpen = path;
	string fileExtension = ".stl";
	// define stringstream, used to add file extensions to the fileToOpen raw path
	stringstream ss;

	// create bt from stl
	stlToBt(fileToOpen);

	// open file in binvox
	ss << "\"" << fileToOpen << "\"" << ".binvox.bt";
	openFileInOctovis(ss.str());
	ss.str(string());

	return 0;
}

//OcTree openBtFileWithOctomap(string filePath){
//	return OcTree(filePath);
//}

//TODO was auch immer fehler war falsche einbindung der lib :|
//TODO add conversion from bt to ot to workflow
/**
 * This Method is used to create a OcTree from a ot file.
 * Terefore this method can be used to read a ot file and edit it.
 * 	!!! input and output not jet implemented !!!
 * input:
 * 		string	: path to the .ot file
 * return:
 * 		OcTree*	: OcTree object from .ot file
 */
int openOtFileWithOctomap(void){
	//string pathToFile = openFile();
	//double res = 0.05;  // create empty tree with resolution 0.05 (different from default 0.1 for test)
	//const string& fn = pathToFile;
//	string fn = "/home/louis/Documents/Models/f5-tiger-scale-model-1-64-by-guaro3d/Fuselaje1_1.binvox.bt.ot";
//	string fn = "/home/louis/Documents/Models/f5-tiger-scale-model-1-64-by-guaro3d/Aladerecha.binvox.bt.ot";
//	string fn = "/home/louis/Documents/Models/f5-tiger-scale-model-1-64-by-guaro3d/Fuselaje1.binvox.bt.ot";
//	string fn = "/home/louis/Documents/Models/turbocharger-with-electric-motor-by-bulgakova-tanya/IMPELLER_1.binvox.bt.ot";

	string fn = "/home/louis/Documents/Models/data/testFile82.binvox.bt.ot";
//	double res = .05;
//	OcTree tr(res);
	octomap::OcTree* tree = dynamic_cast<OcTree*>(octomap::OcTree::read(fn));
	OcTree* octrees = dynamic_cast<OcTree*>(tree);

	OcTree btOt("/home/louis/Documents/Models/data/btFiles/testFile0.binvox.bt");
//	btOt.leaf_iterator();

	// create ot with same resolution. should later be used to generate the nodified surface
	double resolutionOfMainOt = octrees->getResolution();
	OcTree otm(resolutionOfMainOt);
	// pointer conversion not sure if good idea,
	// but in this case probably not as bad as object will later on life till the end
	OcTree* ot = &otm;

	compareOtSpecs(octrees, ot);

	//vector<vector<float>> myMap = getZFromBbx(point3d(min[0],min[1],min[2]), point3d(max[0],max[1],max[2]), octrees, octrees->getResolution());

	// get metric bbx (values are in mm despite the octomap doc :D)
	double min[3];
	double max[3];
	octrees->getMetricMin(min[0], min[1], min[2]);
	octrees->getMetricMax(max[0], max[1], max[2]);
//
//	// gen. heightmap for whole body
//	vector<vector<float>> myMap = areaToHeightmap(octrees, 0, 0, 256, 256, octrees->getResolution());
//
//
////	ot->insertPointCloud(keysFromOctrees,p3d, pose6d(0,0,0,0,0,0),maxRange,0,1);
//
//	//TODO find scaling factor
//	//cout << "maxvals: " << myMap[128][110] << endl; // this is max should be approx 37.9683
//    Gnuplot gp;
//    gp << "unset key\n";
//    gp << "set pm3d\n";
//    gp << "set hidden3d\n";
//    gp << "set view map\n";
//    gp << "set xrange [ -20 : 270 ] \n";
//    gp << "set yrange [ -20 : 270 ] \n";
//    //gp << "set zrange [  31000 : 37400  ] \n";
//    //gp << "set zrange [  -10 : -9.8  ] \n";
//    gp << "splot '-'\n";
//    gp.send2d(myMap);
//    gp.flush();
//	cout << "gnuplot ended?" << endl;
//
//	addVector2dToOctree(myMap, &otm, point3d(0,0,0));

	OcTree otree = createCircularOctree(4,.1,point3d(5,5,5));
	Pointcloud poicl = createCircularPointCloud(8,.1);
	otree.insertPointCloud(poicl,point3d(1,2,2),pose6d(3,3,3,0,0,0),-1,0,0);
	// insert a second circle
	KeySet free;
	KeySet occp;
//	otree.computeUpdate(poicl,point3d(2,2,2),free,occp,-1);
//	(&occp)->
//	otree.computeUpdate(poicl,point3d(4,4,0),free,occp,-1);
	// (poicl,point3d(2,2,2),pose6d(0,0,0,0,0,0),-1,0,0)
//	otree.setNodeValue(free.begin(),1,0);

	otree.updateInnerOccupancy();
	string otmPath = getOutputFilePath();
	otree.write(otmPath);
	cout << "otm written to: " << otmPath << endl;

//	otm.updateInnerOccupancy();
//	string otmPath = getOutputFilePath();
//	otm.write(otmPath);
//	cout << "otm written to: " << otmPath << endl;


    //onode-
	cout << "metric min " << doubleArrToString(min) << " metric max " << doubleArrToString(max) << endl;
	cout << "octree loaded" << endl;
	return 0;
}

/**
 * Only works if empty nodes are non existant / unknown
 * Inserts points from ray to pointcloud
 */
void addKeysFromRay(KeyRay ray, Pointcloud& keysFromOctrees, OcTree* octree){
	auto iter = ray.begin();
	auto kk = *(iter);
	float x,y,z;
	for (int i = 0; i < ray.size(); i++){
		kk = *(iter);
		if (octree->search(kk)){
			x = octree->keyToCoord(kk.k[0]);
			y = octree->keyToCoord(kk.k[1]);
			z = octree->keyToCoord(kk.k[2]);
#pragma omp critical
			keysFromOctrees.push_back(x, y, z);
		}
		iter ++;
	}
}

/**
 * This method is used to create a octree representation of the inserted height map.
 *	input:
 *		vector<vector<float>>	: heightMap: the height or depth map, which values are inserted
 *									in to the octree.
 *		float					: resolution: the metric difference between steps / indices in the depth map
 *	return:
 *		OcTree					: OcTree object with values from heightMap inserted
 */
OcTree vector2dToOctree(vector<vector<float>> heightMap, float resolution, point3d org = point3d(0,0,0)){
	OcTree ot(resolution);
	Pointcloud pc;
	for (double i = 0; i < heightMap.size(); i++){
		for(double j = 0; j < heightMap.size(); j++){
			// without the *.999999 the resulting octree contains gaps
			pc.push_back(point3d(ot.getResolution()*i*.999999,ot.getResolution()*j*.999999,heightMap[i][j]));
		}
	}
	ot.insertPointCloud(pc, org, pose6d(0,0,0,0,0,0), -1, 0, 0);
	return ot;
}

/**
 * This method is used to insert a depth map in an existing octree.
 * As this method is somewhat also one kind of wrapper of octomaps inserPointCloud method, this may corrupt
 * existing point entries.
 *	input:
 *		vector<vector<float>>	: heightMap: the height or depth map, which values are inserted
 *									in to the octree.
 *		float					: resolution: the metric difference between steps / indices in the depth map
 *		point3d					: org: starting position of the octree structure, defualt: (0,0,0)
 */
void addVector2dToOctree(vector<vector<float>> heightMap, OcTree* ot, point3d org = point3d(0,0,0)){
	Pointcloud pc;
	for (double i = 0; i < heightMap.size(); i++){
		for(double j = 0; j < heightMap.size(); j++){
			pc.push_back(point3d(ot->getResolution()*i*.999999,ot->getResolution()*j*.999999,heightMap[i][j]));
		}
	}
	ot->insertPointCloud(pc, org, pose6d(0,0,0,0,0,0), -1, 0, 0);
}

/**
 * Testmethod for creating a circular body in point cloud form
 */
octomap::Pointcloud createCircularPointCloud(float radius, float resolution) {
// circular function: x^2+y^2 <= r
	Pointcloud pc;
	for (float x = -radius; x < radius; x += resolution) {
		for (float y = -radius; y < radius; y += resolution) {
			if ((x * x) + (y * y) <= radius) {
				pc.push_back(point3d(x, y, 0));
			}
		}
	}
	return pc;
}

/**
 * Testmethod for creating a circular body in octree form
 */
OcTree createCircularOctree(float radius, float resolution, point3d center){
	// circular function: x^2+y^2 <= r
	OcTree ot(resolution);
	Pointcloud pc;
	for (float x=-radius; x < radius; x+= resolution){
		for(float y=-radius; y < radius; y+= resolution){
			if ( (x*x) + (y*y) <= radius ){
				pc.push_back(point3d(x, y, 0));
			}
		}
	}
	ot.insertPointCloud(pc, center, pose6d(0,0,0,0,0,0), -1, 0, 0);
	return ot;
}

// TODO insert some sort of linerainterpolation between new and old point like insert line
void insertTool(OcTree* ot, Pointcloud pc, pose6d position){
	ot->insertPointCloud(pc,point3d(),position,-1,0,0);
}

// TODO linerinterpolate between two points and insert toll in steps of octree resolution
/**
 * Method to insert a tool represented as a pointcloud in an octree.
 * Number of steps needed from beginning to the end is rounded down as the variable is truncated to a int.
 */
void insertToolInLine(OcTree* ot, Pointcloud pc, point3d startPosition,
		point3d endPosition) {
	// calculate differences
	auto dx = endPosition.x() - startPosition.x();
	auto dy = endPosition.y() - startPosition.y();
	auto dz = endPosition.z() - startPosition.z();
	// scale differences to min raster diffference / ot resolution
	auto mx = max(max(abs(dx), abs(dy)), abs(dz));
	if (abs(mx) > 0) {
		// point cloud steps need to be a little bit finer than ot resolution
		// otherwise there will be gaps somehow
		auto scalingFactor = ot->getResolution()*.9999999 / mx;
		auto dxs = dx * scalingFactor;
		auto dys = dy * scalingFactor;
		auto dzs = dz * scalingFactor;
//		cout << dxs << " :: " << dys << " :: " << dzs << endl;
		auto mxs = max(max(abs(dxs), abs(dys)), abs(dzs));
		int numOfSteps = mx / mxs;
//		cout << "num of steps: " << numOfSteps << endl;
		// if no steps are needed simply add the tool
		if (numOfSteps == 0) {
			cout << "no stesps" << endl;
			ot->insertPointCloud(pc, point3d(),
					pose6d(startPosition.x(), startPosition.y(),
							startPosition.z(), 0, 0, 0), -1, 0, 0);
		} else {
			// else insert the tool "as line"
			for (int i = 0; i < numOfSteps; i++) {
				ot->insertPointCloud(pc, point3d(),
						pose6d(startPosition.x() + dxs * i,
								startPosition.y() + dys * i,
								startPosition.z() + dzs * i, 0, 0, 0), -1, 0,
						0);
			}
		}
	}
}

/**
 * Method to insert a tool represented as a pointcloud in an octree.
 * Number of steps needed from beginning to the end is rounded up.
 * This way all voxel hit at least slightly are "machined".
 * This Method got implemented as other function for speed reasons, saving
 * one mode variable check.
 */
void insertToolInLineRoundUp(OcTree* ot, Pointcloud pc, point3d startPosition,
		point3d endPosition) {
	// calculate differences
	auto dx = endPosition.x() - startPosition.x();
	auto dy = endPosition.y() - startPosition.y();
	auto dz = endPosition.z() - startPosition.z();
	// scale differences to min raster diffference / ot resolution
	auto mx = max(max(abs(dx), abs(dy)), abs(dz));
	if (abs(mx) > 0) {
		// point cloud steps need to be a little bit finer than ot resolution
		// otherwise there will be gaps somehow
//		auto scalingFactor = ot->getResolution()*.9999999 / mx;
		// smaller steps may reduce uncertainty -> does so, but increases calculation time
		auto scalingFactor = ot->getResolution()*.4 / mx;
		auto dxs = dx * scalingFactor;
		auto dys = dy * scalingFactor;
		auto dzs = dz * scalingFactor;
//		cout << dxs << " :: " << dys << " :: " << dzs << endl;
		auto mxs = max(max(abs(dxs), abs(dys)), abs(dzs));
		int numOfSteps = ceil(mx / mxs);
//		cout << "num of steps: " << numOfSteps << endl;
		// if no steps are needed simply add the tool
		if (numOfSteps == 0) {
			ot->insertPointCloud(pc, point3d(),
					pose6d(startPosition.x(), startPosition.y(),
							startPosition.z(), 0, 0, 0), -1, 0, 0);
		} else {
			// else insert the tool "as line"
			for (int i = 0; i < numOfSteps; i++) {
				ot->insertPointCloud(pc, point3d(),
						pose6d(startPosition.x() + dxs * i,
								startPosition.y() + dys * i,
								startPosition.z() + dzs * i, 0, 0, 0), -1, 0,
						0);
			}
		}
	}
}

/**
 * Method to insert a tool represented as a pointcloud in an octree.
 * Number of steps needed from beginning to the end is rounded up.
 * This way all voxel hit at least slightly are "machined".
 * This Method got implemented as other function for speed reasons, saving
 * one mode variable check.
 * This Method uses a smaler stepsize.
 * This may results in longer calculation time, but increases certainty in the created octree nodes
 * and in case of diagonal paths the chance of slightly cut nodes is increased.
 */
void insertToolInLineRoundUpMiniSteps(OcTree* ot, Pointcloud pc, point3d startPosition,
		point3d endPosition) {
	// calculate differences
	auto dx = endPosition.x() - startPosition.x();
	auto dy = endPosition.y() - startPosition.y();
	auto dz = endPosition.z() - startPosition.z();
	// scale differences to min raster diffference / ot resolution
	auto mx = max(max(abs(dx), abs(dy)), abs(dz));
	if (abs(mx) > 0) {
		// point cloud steps need to be a little bit finer than ot resolution
		// otherwise there will be gaps somehow
//		auto scalingFactor = ot->getResolution()*.9999999 / mx;
		// smaller steps may reduce uncertainty -> does so, but increases calculation time
		auto scalingFactor = ot->getResolution()*.4 / mx;
		auto dxs = dx * scalingFactor;
		auto dys = dy * scalingFactor;
		auto dzs = dz * scalingFactor;
//		cout << dxs << " :: " << dys << " :: " << dzs << endl;
		auto mxs = max(max(abs(dxs), abs(dys)), abs(dzs));
		int numOfSteps = ceil(mx / mxs);
//		cout << "num of steps: " << numOfSteps << endl;
		// if no steps are needed simply add the tool
		if (numOfSteps == 0) {
			ot->insertPointCloud(pc, point3d(),
					pose6d(startPosition.x(), startPosition.y(),
							startPosition.z(), 0, 0, 0), -1, 0, 0);
		} else {
			// else insert the tool "as line"
			for (int i = 0; i < numOfSteps; i++) {
				ot->insertPointCloud(pc, point3d(),
						pose6d(startPosition.x() + dxs * i,
								startPosition.y() + dys * i,
								startPosition.z() + dzs * i, 0, 0, 0), -1, 0,
						0);
			}
		}
	}
}


void insertToolInLineUsingOctomapComputeRayKeys(octomap::OcTree* ot, octomap::Pointcloud pc, octomap::point3d startPosition, octomap::point3d endPosition){
	KeyRay ray;
	const point3d pt1 = startPosition;
	const point3d pt2 = endPosition;
	vector<point3d> pointsList;
	ot->computeRayKeys(pt1, pt2, ray);
	for (auto key : ray){
		point3d rayPt = ot->keyToCoord(key);
		ot->insertPointCloud(pc, point3d(),
				pose6d(rayPt.x(),
						rayPt.y(),
						rayPt.z(), 0, 0, 0), -1, 0, 0);
	}
	// the last point needs to be made less accurate to match the rest of the insert,
	// otherwise start and endpoint offset calculation and linerinterpolation
	// for every point would be needed to keep the "original" positions
	auto newP = ot->keyToCoord(ot->coordToKey(pt2));
	// computeRayKeys excludes end point, therefore insert also at endPosition.
	ot->insertPointCloud(pc, point3d(),
			pose6d(newP.x(),
				   newP.y(),
				   newP.z(), 0, 0, 0), -1, 0, 0);
}


/**
 * Only works if empty nodes are non existant / unknown
 */
float getMaxZFromRay(KeyRay ray, OcTree* octree){
	auto iter = ray.begin();
	auto kk = *(iter);
	double min[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	float zMax=min[2];
//	float x,y,z;
	for (int i = 0; i < ray.size(); i++){
		kk = *(iter);
		if (octree->search(kk)!=NULL){
			point3d coord = octree->keyToCoord(kk);
			if(coord.z()>zMax){
				// test TODO
//				return zMax;
				zMax = coord.z();
			}
		}
		iter ++;
	}
	return zMax;
}

float getMaxZFromRayWithAlternateMin(KeyRay ray, OcTree* octree, double minZ){
	auto iter = ray.begin();
	auto kk = *(iter);
	float zMax=minZ;
//	float x,y,z;
	for (int i = 0; i < ray.size(); i++){
		kk = *(iter);
		if (octree->search(kk)!=NULL){
			point3d coord = octree->keyToCoord(kk);
			if(coord.z()>zMax){
				// test TODO
//				return zMax;
				zMax = coord.z();
			}
		}
		iter ++;
	}
	return zMax;
}

float getMinZFromRay(KeyRay ray, OcTree* octree){
	auto iter = ray.begin();
	auto kk = *(iter);
	double max[3];
	octree->getMetricMax(max[0],max[1],max[2]);
	float zMin=max[2];
//	float x,y,z;
	for (int i = 0; i < ray.size(); i++){
		kk = *(iter);
		if (octree->search(kk)!=NULL){
			point3d coord = octree->keyToCoord(kk);
			if(coord.z()<zMin){
				// test TODO
//				return zMax;
				zMin= coord.z();
			}
		}
		iter ++;
	}
	return zMin;
}

/**
 * Returns on first node found -> may increase speed.
 */
float getMinZFromRayOptimised(KeyRay ray, OcTree* octree, float minZ, float inBBXZ, bbx2D bbx){
	double max[3];
	octree->getMetricMax(max[0], max[1], max[2]);
	float zMin = max[2];
	bool foundNode = false;

	if (octree->keyToCoord(*(ray.end())).z()
			< octree->keyToCoord(*(ray.begin())).z()) {

		// if ray end is below ray start invert ray
		auto iter = ray.end();
		auto kk = *(iter);
		//	float x,y,z;
		auto rb = ray.begin();
			while (iter != rb) {
				kk = *(iter);
				if (octree->search(kk) != NULL) {
					point3d coord = octree->keyToCoord(kk);
					if (coord.z() < zMin) {
						// test TODO
		//				return zMax;
						// return
						zMin = coord.z();
						foundNode = true;
					}
				}
				iter--;
			}
	} else {
		auto iter = ray.begin();
		auto kk = *(iter);
		//	float x,y,z;
		auto re = ray.end();
			while (iter != re) {
				kk = *(iter);
				if (octree->search(kk) != NULL) {
					point3d coord = octree->keyToCoord(kk);
					if (coord.z() < zMin) {
						// test TODO
		//				return zMax;
						// return
						zMin = coord.z();
						foundNode = true;
					}
				}
				iter++;
			}
	}
	if (foundNode) {
		return zMin;
	}

	// if no node was found check if point is in bbx or not and return coresponding value
	// otk is the starting point of the ray
	auto otk = octree->keyToCoord((*ray.begin()));
	if ((otk.x() >= bbx.xMin) && (otk.x() <= bbx.xMax) && (otk.y() >= bbx.yMin)
			&& (otk.y() <= bbx.yMax)) {
		return inBBXZ;
	} else {
		// if point is outside of bbx set to min z
		return minZ;
	}

//	return zMin;
}

/**
 * This method is used to create a 2D-heightmap from an octree for a defined area
 * input:
 * 		octree	:	octree from which to create the heightmap
 * 		float	:	xStartpos: x start position of area to create heightmap
 * 						this position value is relative to min of octree / object -> 0 = begin of obj
 * 		flaot	:	yStartpos: y start position equal to xStartpos but for y
 * 		int		:	xSize: size of heightmap in x direction. Equals to number of steps taken in x direction from
 * 						xStartpos
 * 		int		:	ySize: same as xSize, but for y dimension.
 * 		float	:	stepSize: size per step in every direction (x and y equally)
 * 	return:
 * 		vector<vector<float>>	:	heightmap in 2d-vector format (height is in absolute coordinate)
 */
vector<vector<float>> areaToHeightmap(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
//	float xoffs = min[0] + xStartpos;
//	float yoffs = min[1] + yStartpos;
	// calc max vas for num of rays
//	float* maxVals = new float[sl1][sl2];
	typedef	boost::multi_array<float, 2> arr2D;
	typedef arr2D::index index;
	arr2D maxVals(boost::extents[xSize][ySize]);

	chrono::microseconds usBeforeRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	// calc zvals for rays in parallel (using maxVals array for parallel op)
#pragma omp parallel for
	for (index i = 0; i < xSize; i++){
#pragma omp parallel for
		for (index j = 0; j < ySize; j++){
			KeyRay ray;
			octree->computeRayKeys(point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,min[2]-1),point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,max[2]+1), ray);
			maxVals[i][j] = getMaxZFromRay(ray, octree);
			ray.reset();
		}
	}

	if (BENCH_HEIGHTMAP_GEN){
		chrono::microseconds usAfterRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: " << (usAfterRayCast.count()-usBeforeRayCast.count()) << " us or: " << (usAfterRayCast.count()-usBeforeRayCast.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}
	return multiArray2DToVector(maxVals);
}

/**
 * This method is used to create a 2D-heightmap from an octree for a defined area
 * input:
 * 		octree	:	octree from which to create the heightmap
 * 		float	:	xStartpos: x start position of area to create heightmap
 * 						this position value is relative to min of octree / object -> 0 = begin of obj
 * 		flaot	:	yStartpos: y start position equal to xStartpos but for y
 * 		int		:	xSize: size of heightmap in x direction. Equals to number of steps taken in x direction from
 * 						xStartpos
 * 		int		:	ySize: same as xSize, but for y dimension.
 * 		float	:	stepSize: size per step in every direction (x and y equally)
 * 	return:
 * 		vector<vector<float>>	:	heightmap in 2d-vector format (height is in absolute coordinate)
 */
vector<vector<float>> areaToHeightmapAlternateMinZ(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize, double minZ){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
//	float xoffs = min[0] + xStartpos;
//	float yoffs = min[1] + yStartpos;
	// calc max vas for num of rays
//	float* maxVals = new float[sl1][sl2];
	typedef	boost::multi_array<float, 2> arr2D;
	typedef arr2D::index index;
	arr2D maxVals(boost::extents[xSize][ySize]);

	chrono::microseconds usBeforeRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	// calc zvals for rays in parallel (using maxVals array for parallel op)
#pragma omp parallel for
	for (index i = 0; i < xSize; i++){
#pragma omp parallel for
		for (index j = 0; j < ySize; j++){
			KeyRay ray;
			octree->computeRayKeys(point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,min[2]-1),point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,max[2]+1), ray);
			maxVals[i][j] = getMaxZFromRayWithAlternateMin(ray, octree, minZ);
			ray.reset();
		}
	}

	if (BENCH_HEIGHTMAP_GEN){
		chrono::microseconds usAfterRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: " << (usAfterRayCast.count()-usBeforeRayCast.count()) << " us or: " << (usAfterRayCast.count()-usBeforeRayCast.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}
	return multiArray2DToVector(maxVals);
}

vector<vector<float>> areaToHeightmapRelativeToBody(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
	float xoffs = min[0] + xStartpos;
	float yoffs = min[1] + yStartpos;
	return areaToHeightmap(octree, xoffs, yoffs, xSize, ySize, stepSize);
}

/**
 * This method is used to create a 2D-heightmap from an octree for a defined area
 * input:
 * 		octree	:	octree from which to create the heightmap
 * 		float	:	xStartpos: x start position of area to create heightmap
 * 						this position value is relative to min of octree / object -> 0 = begin of obj
 * 		flaot	:	yStartpos: y start position equal to xStartpos but for y
 * 		int		:	xSize: size of heightmap in x direction. Equals to number of steps taken in x direction from
 * 						xStartpos
 * 		int		:	ySize: same as xSize, but for y dimension.
 * 		float	:	stepSize: size per step in every direction (x and y equally)
 * 	return:
 * 		vector<vector<float>>	:	heightmap in 2d-vector format (height is in absolute coordinate)
 */
vector<vector<float>> areaToDepthmap(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
//	float xoffs = min[0] + xStartpos;
//	float yoffs = min[1] + yStartpos;
	// calc max vas for num of rays
//	float* maxVals = new float[sl1][sl2];
	typedef	boost::multi_array<float, 2> arr2D;
	typedef arr2D::index index;
	arr2D maxVals(boost::extents[xSize][ySize]);

	chrono::microseconds usBeforeRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	// calc zvals for rays in parallel (using maxVals array for parallel op)
#pragma omp parallel for
	for (index i = 0; i < xSize; i++){
#pragma omp parallel for
		for (index j = 0; j < ySize; j++){
			KeyRay ray;
			octree->computeRayKeys(point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,min[2]-1),point3d(stepSize*i+xStartpos,stepSize*j+yStartpos,max[2]+1), ray);
			maxVals[i][j] = getMinZFromRay(ray, octree);
			ray.reset();
		}
	}

	if (BENCH_HEIGHTMAP_GEN){
		chrono::microseconds usAfterRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: " << (usAfterRayCast.count()-usBeforeRayCast.count()) << " us or: " << (usAfterRayCast.count()-usBeforeRayCast.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}
	return multiArray2DToVector(maxVals);
}

/**
 * BBX is including, meaning if value outside -> minZ if on the other.
 * This method is used to create a 2D-heightmap from an octree for a defined area
 * input:
 * 		octree	:	octree from which to create the heightmap
 * 		float	:	xStartpos: x start position of area to create heightmap
 * 						this position value is relative to min of octree / object -> 0 = begin of obj
 * 		flaot	:	yStartpos: y start position equal to xStartpos but for y
 * 		int		:	xSize: size of heightmap in x direction. Equals to number of steps taken in x direction from
 * 						xStartpos
 * 		int		:	ySize: same as xSize, but for y dimension.
 * 		float	:	stepSize: size per step in every direction (x and y equally)
 * 	return:
 * 		vector<vector<float>>	:	heightmap in 2d-vector format (height is in absolute coordinate)
 */
vector<vector<float>> areaToDepthmapWithBBX(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float bbxMinX, float bbxMinY, float bbxMaxX, float bbxMaxY){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
//	float xoffs = min[0] + xStartpos;
//	float yoffs = min[1] + yStartpos;
	// calc max vas for num of rays
//	float* maxVals = new float[sl1][sl2];
	typedef	boost::multi_array<float, 2> arr2D;
	typedef arr2D::index index;
	arr2D maxVals(boost::extents[xSize][ySize]);

	chrono::microseconds usBeforeRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	// calc zvals for rays in parallel (using maxVals array for parallel op)
#pragma omp parallel for
	for (index i = 0; i < xSize; i++){
#pragma omp parallel for
		for (index j = 0; j < ySize; j++){
			float xPos = stepSize*i+xStartpos;
			float yPos = stepSize*j+yStartpos;
			// if point is in bbx calculate, else set value of cell to minZ
			if ((xPos >= bbxMinX) && (xPos <= bbxMaxX) && (yPos >= bbxMinY) && (yPos <= bbxMaxY)){
				KeyRay ray;
				octree->computeRayKeys(point3d(xPos,yPos,min[2]-1),point3d(xPos,yPos,max[2]+1), ray);
				maxVals[i][j] = getMinZFromRay(ray, octree);
				ray.reset();
			} else {
				maxVals[i][j] = minZ;
			}
		}
	}

	if (BENCH_HEIGHTMAP_GEN){
		chrono::microseconds usAfterRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: " << (usAfterRayCast.count()-usBeforeRayCast.count()) << " us or: " << (usAfterRayCast.count()-usBeforeRayCast.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}
	return multiArray2DToVector(maxVals);
}

/**
 * BBX is including, meaning if value outside -> minZ if on the other.
 * This method is used to create a 2D-heightmap from an octree for a defined area
 * input:
 * 		octree	:	octree from which to create the heightmap
 * 		float	:	xStartpos: x start position of area to create heightmap
 * 						this position value is relative to min of octree / object -> 0 = begin of obj
 * 		flaot	:	yStartpos: y start position equal to xStartpos but for y
 * 		int		:	xSize: size of heightmap in x direction. Equals to number of steps taken in x direction from
 * 						xStartpos
 * 		int		:	ySize: same as xSize, but for y dimension.
 * 		float	:	stepSize: size per step in every direction (x and y equally)
 * 	return:
 * 		vector<vector<float>>	:	heightmap in 2d-vector format (height is in absolute coordinate)
 */
vector<vector<float>> areaToDepthmapWithBBXOptimised(OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float inBBXZ, bbx2D bbx){
	// metic min and max
	double min[3];
	double max[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	octree->getMetricMax(max[0], max[1], max[2]);
	// make lowest x and y (0/0) -> transform in object coordinates
//	float xoffs = min[0] + xStartpos;
//	float yoffs = min[1] + yStartpos;
	// calc max vas for num of rays
//	float* maxVals = new float[sl1][sl2];
//	typedef	boost::multi_array<float, 2> arr2D;
//	typedef arr2D::index index;
//	typedef fMA2DIndex index;
	fMultiArr2D maxVals(boost::extents[xSize][ySize]);
	typedef fMultiArr2D::index index;

	chrono::microseconds usBeforeRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());

	// calc zvals for rays in parallel (using maxVals array for parallel op)
#pragma omp parallel for
	for (index i = 0; i < xSize; i++){
#pragma omp parallel for
		for (index j = 0; j < ySize; j++){
			float xPos = stepSize*i+xStartpos;
			float yPos = stepSize*j+yStartpos;

			// if point is in octree (f.e. machined ot) calculate minZ with ray
			// else set to defined value
			if ((xPos <= max[0]) && (xPos >= min[0]) && (yPos <= max[1]) && (yPos >= min[1])){

				KeyRay ray;
				octree->computeRayKeys(point3d(xPos,yPos,min[2]-1),point3d(xPos,yPos,max[2]+1), ray);
				maxVals[i][j] = getMinZFromRayOptimised(ray, octree, minZ, inBBXZ, bbx);
				ray.reset();

			} else {

				// if point is in bbx set to z val for in bbx
				if ((xPos >= bbx.xMin) && (xPos <= bbx.xMax) && (yPos >= bbx.yMin) && (yPos <= bbx.yMax)){
					maxVals[i][j] = inBBXZ;
				} else {
					// if point is outside of bbx set to min z
					maxVals[i][j] = minZ;
				}

			}
		}
	}

	if (BENCH_HEIGHTMAP_GEN){
		chrono::microseconds usAfterRayCast = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: " << (usAfterRayCast.count()-usBeforeRayCast.count()) << " us or: " << (usAfterRayCast.count()-usBeforeRayCast.count())/1000000 << " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}
	return multiArray2DToVector(maxVals);
}


vector<vector<float>> areaToHeightmapWithLeafIterator(octomap::OcTree* octree,
		float xStartpos, float yStartpos, int xSize, int ySize, float stepSize,
		float minZ) {

	// benchmark
	chrono::microseconds usBeforeStateGen = chrono::duration_cast<
			chrono::microseconds>(
			chrono::system_clock::now().time_since_epoch());

	fMultiArr2D maxVals(boost::extents[xSize][ySize]);
	fill(maxVals.data(), maxVals.data() + maxVals.num_elements(), minZ);

	long numNodes = 0;
	long xIdx = 0;
	long yIdx = 0;

	for (auto iter = octree->begin_leafs(), end = octree->end_leafs();
			iter != end; iter++) {

		xIdx = (long) ((iter.getCoordinate().x() - xStartpos
				+ iter.getSize() / 2) / stepSize);
		yIdx = (long) ((iter.getCoordinate().y() - yStartpos
				+ iter.getSize() / 2) / stepSize);

		auto xMinIdx = (long) floor(xIdx - (iter.getSize() / 2 / stepSize));
		auto yMinIdx = (long) floor(yIdx - (iter.getSize() / 2 / stepSize));
		auto xMaxIdx = (long) ceil(xIdx + (iter.getSize() / 2 / stepSize));
		auto yMaxIdx = (long) ceil(yIdx + (iter.getSize() / 2 / stepSize));
		if (xMinIdx < 0) {
			xMinIdx = 0;
		}
		if (yMinIdx < 0) {
			yMinIdx = 0;
		}
		if (xMaxIdx >= xSize) {
			xMaxIdx = xSize - 1;
		}
		if (yMaxIdx >= ySize) {
			yMaxIdx = ySize - 1;
		}

		xIdx = xMinIdx;

		for (; yMinIdx <= yMaxIdx; yMinIdx++) {
			xMinIdx = xIdx;
			for (; xMinIdx <= xMaxIdx; xMinIdx++) {
				if (maxVals[xMinIdx][yMinIdx]
						< iter.getCoordinate().z() + iter.getSize() / 2) {
					maxVals[xMinIdx][yMinIdx] = iter.getCoordinate().z()
							+ iter.getSize() / 2;
				}
			}
		}
	}

	// benchmark
	chrono::microseconds usAfterStateGen = chrono::duration_cast<
			chrono::microseconds>(
			chrono::system_clock::now().time_since_epoch());
	cout << "Raycast took: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count())
			<< " us or: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000
			<< " ms or: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000000
			<< " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")

	return multiArray2DToVector(maxVals);
}


/**
 * Method to retrieve a heightmap using the leaf iterator method.
 * For small trees this method can be significantly faster the
 * method without leaf iterator, but raytracing.
 * In ONE test this leaf iterator method was approximately
 * 50 times faster.
 * input:
 * 		OcTree*		: octree: octree(pointer(&OcTree)) from which to retrieve the z-values if needed.
 * 		float		: xStartpos: metric startposition(absolut pos) of the generated heightmap.
 * 		float		: yStartpos: metric startposition(absolut pos) of the generated heightmap.
 * 		int			: xSize: number of steps in x-direction / size of heightmap in x-dimension.
 * 		int			: ySize: number of steps in y-direction / size of heightmap in y-dimension.
 * 		float		: stepSize: metric distance between one point in the heightmap and an other.
 * 						Distance made per step.
 * 		float		: minZ: minimum z value for areas which are neither in the bbx nor in the object.
 * 		float		: inBBXZ: z-value for points which are in the bbx, but not in the object.
 * 		bbx2D		: bbx: struct with min and max coordinates of the bbx. Values are metric and absolute.
 * return:
 * 		vector<vector<float>> : areaToHeightmapWithLeafIteratorAndBBX: 2D-vector with heightmap / z-values.
 */
vector<vector<float>> areaToHeightmapWithLeafIteratorAndBBX(
		octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize,
		int ySize, float stepSize, float minZ, float inBBXZ, bbx2D bbx) {

	double min[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	double max[3];
	octree->getMetricMax(max[0], max[1], max[2]);
	// only calculate values for nodes in bbx from startpos to endpos
	// / startpos + stepSize*$numSteps

	float xEndpos = xStartpos + xSize * stepSize;
	float yEndpos = yStartpos + ySize * stepSize;

	chrono::microseconds usBeforeStateGen = chrono::duration_cast<
			chrono::microseconds>(
			chrono::system_clock::now().time_since_epoch());

	// create boost multiarray and define index
	fMultiArr2D maxVals(boost::extents[xSize][ySize]);
	typedef fMultiArr2D::index index;


	point3d bbxMin(xStartpos, yStartpos, min[2]);
	point3d bbxMax(xEndpos, yEndpos, max[2]);

	if (bbxMax.x() > bbx.xMax) {
		bbxMax.x() = bbx.xMax;
	}
	if (bbxMin.x() < bbx.xMin) {
		bbxMin.x() = bbx.xMin;
	}
	if (bbxMax.y() > bbx.yMax) {
		bbxMax.y() = bbx.yMax;
	}
	if (bbxMin.y() < bbx.yMin) {
		bbxMin.y() = bbx.yMin;
	}

	long idxOfBBXxMin = (long) (floor((bbxMin.x() - xStartpos) / stepSize));
	long idxOfBBXxMax = (long) (ceil((bbxMax.x() - xStartpos) / stepSize));
	long idxOfBBXyMin = (long) floor(((bbxMin.y() - yStartpos) / stepSize));
	long idxOfBBXyMax = (long) ceil(((bbxMax.y() - yStartpos) / stepSize));
	if (idxOfBBXxMin < 0) {
		idxOfBBXxMin = 0;
	}
	if (idxOfBBXxMax > xSize - 1) {
		idxOfBBXxMax = xSize - 1;
	}
	if (idxOfBBXyMin < 0) {
		idxOfBBXyMin = 0;
	}
	if (idxOfBBXyMax > ySize - 1) {
		idxOfBBXyMax = ySize - 1;
	}

	for (index i = 0; i < maxVals.shape()[0]; i++) {
		for (index j = 0; j < maxVals.shape()[1]; j++) {
			// if point is in bbx set to in BBXZ, else set to zMin
			if ((i >= idxOfBBXxMin) && (i <= idxOfBBXxMax)
					&& (j >= idxOfBBXyMin) && (j <= idxOfBBXyMax)) {
				maxVals[i][j] = inBBXZ;
			} else {
				maxVals[i][j] = minZ;
			}
		}
	}

	// if start or endpoints are out of bbx return
	if (((xEndpos <= bbxMin.x()) && (yEndpos <= bbxMin.y()))
			|| ((xStartpos >= bbxMax.x()) && (yStartpos >= bbxMax.y()))) {
		cout << "no leaf iter needed: return" << endl;
		return multiArray2DToVector(maxVals);
	}

	long xIdx = 0;
	long yIdx = 0;

	// loop over all leafs in the optimised and adjusted bbx
	for (auto iter = octree->begin_leafs_bbx(bbxMin, bbxMax), end =
			octree->end_leafs_bbx(); iter != end; iter++) {

		if (octree->isNodeOccupied(*iter)) {

			xIdx = (long) ((iter.getCoordinate().x() - xStartpos
					+ iter.getSize() / 2) / stepSize);
			yIdx = (long) ((iter.getCoordinate().y() - yStartpos
					+ iter.getSize() / 2) / stepSize);

			auto xMinIdx = (long) floor(xIdx - (iter.getSize() / 2 / stepSize));
			auto yMinIdx = (long) floor(yIdx - (iter.getSize() / 2 / stepSize));
			auto xMaxIdx = (long) ceil(xIdx + (iter.getSize() / 2 / stepSize));
			auto yMaxIdx = (long) ceil(yIdx + (iter.getSize() / 2 / stepSize));
			if (xMinIdx < idxOfBBXxMin) {
				xMinIdx = idxOfBBXxMin;
			}
			if (yMinIdx < idxOfBBXyMin) {
				yMinIdx = idxOfBBXyMin;
			}
			if (xMaxIdx >= idxOfBBXxMax) {
				xMaxIdx = idxOfBBXxMax;
			}
			if (yMaxIdx >= idxOfBBXyMax) {
				yMaxIdx = idxOfBBXyMax;
			}

			xIdx = xMinIdx;

			for (; yMinIdx <= yMaxIdx; yMinIdx++) {
				xMinIdx = xIdx;
				for (; xMinIdx <= xMaxIdx; xMinIdx++) {
					// if new z value is bigger than the one in maxVals replace value in maxVals
					if (maxVals[xMinIdx][yMinIdx]
							< iter.getCoordinate().z() + iter.getSize() / 2) {
						maxVals[xMinIdx][yMinIdx] = iter.getCoordinate().z()
								+ iter.getSize() / 2;
					}
				}
			}
		}
	}

	if (BENCH_DEPTHMAP_GEN){
	// benchmark
	chrono::microseconds usAfterStateGen = chrono::duration_cast<
			chrono::microseconds>(
			chrono::system_clock::now().time_since_epoch());
	cout << "Raycast took: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count())
			<< " us or: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000
			<< " ms or: "
			<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000000
			<< " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

	return multiArray2DToVector(maxVals);
}


/**
 * Method to retrieve a heightmap using the leaf iterator method.
 * For small trees this method can be significantly faster the
 * method without leaf iterator, but raytracing.
 * In ONE test this leaf iterator method was approximately
 * 50 times faster.
 * input:
 * 		OcTree*		: octree: octree(pointer(&OcTree)) from which to retrieve the z-values if needed.
 * 		float		: xStartpos: metric startposition(absolut pos) of the generated heightmap.
 * 		float		: yStartpos: metric startposition(absolut pos) of the generated heightmap.
 * 		int			: xSize: number of steps in x-direction / size of heightmap in x-dimension.
 * 		int			: ySize: number of steps in y-direction / size of heightmap in y-dimension.
 * 		float		: stepSize: metric distance between one point in the heightmap and an other.
 * 						Distance made per step.
 * 		float		: minZ: minimum z value for areas which are neither in the bbx nor in the object.
 * 		float		: inBBXZ: z-value for points which are in the bbx, but not in the object.
 * 		bbx2D		: bbx: struct with min and max coordinates of the bbx. Values are metric and absolute.
 * return:
 * 		vector<vector<float>> : areaToHeightmapWithLeafIteratorAndBBX: 2D-vector with heightmap / z-values.
 */
vector<vector<float>> areaToDepthmapWithLeafIteratorAndBBX(
		octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize,
		int ySize, float stepSize, float minZ, float inBBXZ, bbx2D bbx) {

	double min[3];
	octree->getMetricMin(min[0], min[1], min[2]);
	double max[3];
	octree->getMetricMax(max[0], max[1], max[2]);
	// only calculate values for nodes in bbx from startpos to endpos
	// / startpos + stepSize*$numSteps

	float xEndpos = xStartpos + xSize * stepSize;
	float yEndpos = yStartpos + ySize * stepSize;

	chrono::microseconds usBeforeStateGen = chrono::duration_cast<
			chrono::microseconds>(
			chrono::system_clock::now().time_since_epoch());

	// create boost multiarray and define index
	fMultiArr2D maxVals(boost::extents[xSize][ySize]);
	typedef fMultiArr2D::index index;

	point3d bbxMin(xStartpos, yStartpos, min[2]);
	point3d bbxMax(xEndpos, yEndpos, max[2]);

	if (bbxMax.x() > bbx.xMax) {
		bbxMax.x() = bbx.xMax;
	}
	if (bbxMin.x() < bbx.xMin) {
		bbxMin.x() = bbx.xMin;
	}
	if (bbxMax.y() > bbx.yMax) {
		bbxMax.y() = bbx.yMax;
	}
	if (bbxMin.y() < bbx.yMin) {
		bbxMin.y() = bbx.yMin;
	}

	long idxOfBBXxMin = (long) (floor((bbxMin.x() - xStartpos) / stepSize));
	long idxOfBBXxMax = (long) (ceil((bbxMax.x() - xStartpos) / stepSize));
	long idxOfBBXyMin = (long) (floor((bbxMin.y() - yStartpos) / stepSize));
	long idxOfBBXyMax = (long) (ceil((bbxMax.y() - yStartpos) / stepSize));
	if (idxOfBBXxMin < 0) {
		idxOfBBXxMin = 0;
	}
	if (idxOfBBXxMax > xSize - 1) {
		idxOfBBXxMax = xSize - 1;
	}
	if (idxOfBBXyMin < 0) {
		idxOfBBXyMin = 0;
	}
	if (idxOfBBXyMax > ySize - 1) {
		idxOfBBXyMax = ySize - 1;
	}

	for (index i = 0; i < maxVals.shape()[0]; i++) {
		for (index j = 0; j < maxVals.shape()[1]; j++) {
			// if point is in bbx set to in BBXZ, else set to zMin
			if ((i >= idxOfBBXxMin) && (i <= idxOfBBXxMax)
					&& (j >= idxOfBBXyMin) && (j <= idxOfBBXyMax)) {
				maxVals[i][j] = inBBXZ;
			} else {
				maxVals[i][j] = minZ;
			}
		}
	}

	// if start or endpoints are out of bbx return
	if (((xEndpos <= bbxMin.x()) && (yEndpos <= bbxMin.y()))
			|| ((xStartpos >= bbxMax.x()) && (yStartpos >= bbxMax.y()))) {
		cout << "no leaf iter needed: return" << endl;
		return multiArray2DToVector(maxVals);
	}

	long xIdx = 0;
	long yIdx = 0;

	// loop over all leafs in the optimised and adjusted bbx
	for (auto iter = octree->begin_leafs_bbx(bbxMin, bbxMax), end =
			octree->end_leafs_bbx(); iter != end; iter++) {

		if (octree->isNodeOccupied(*iter)) {

			xIdx = (long) ((iter.getCoordinate().x() - xStartpos
					+ iter.getSize() / 2) / stepSize);
			yIdx = (long) ((iter.getCoordinate().y() - yStartpos
					+ iter.getSize() / 2) / stepSize);

			auto xMinIdx = (long) floor(xIdx - (iter.getSize() / 2 / stepSize));
			auto yMinIdx = (long) floor(yIdx - (iter.getSize() / 2 / stepSize));
			auto xMaxIdx = (long) ceil(xIdx + (iter.getSize() / 2 / stepSize));
			auto yMaxIdx = (long) ceil(yIdx + (iter.getSize() / 2 / stepSize));
			if (xMinIdx < idxOfBBXxMin) {
				xMinIdx = idxOfBBXxMin;
			}
			if (yMinIdx < idxOfBBXyMin) {
				yMinIdx = idxOfBBXyMin;
			}
			if (xMaxIdx >= idxOfBBXxMax) {
				xMaxIdx = idxOfBBXxMax;
			}
			if (yMaxIdx >= idxOfBBXyMax) {
				yMaxIdx = idxOfBBXyMax;
			}

			xIdx = xMinIdx;

			for (; yMinIdx <= yMaxIdx; yMinIdx++) {
				xMinIdx = xIdx;
				for (; xMinIdx <= xMaxIdx; xMinIdx++) {
					// if new z value is smaller than the one in maxVals replace value in maxVals
					if (maxVals[xMinIdx][yMinIdx]
							> iter.getCoordinate().z() - iter.getSize() / 2) {
						maxVals[xMinIdx][yMinIdx] = iter.getCoordinate().z()
								- iter.getSize() / 2;
					}
				}
			}
		}
	}

	if (BENCH_DEPTHMAP_GEN){
		// benchmark
		chrono::microseconds usAfterStateGen = chrono::duration_cast<
				chrono::microseconds>(
				chrono::system_clock::now().time_since_epoch());
		cout << "Raycast took: "
				<< (usAfterStateGen.count() - usBeforeStateGen.count())
				<< " us or: "
				<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000
				<< " ms or: "
				<< (usAfterStateGen.count() - usBeforeStateGen.count()) / 1000000
				<< " s" << endl; // @suppress("Method cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Invalid overload")
	}

	return multiArray2DToVector(maxVals);
}


/**
 * decisions for dimension are made when creating the original dm
 * totally not optimized code yet :D
 * This method is used to create a depth map for a part area of another depth map.
 * For points outside of the dm the values of the returned depth map are minZ.
 * For points inside the area of the original depthmap the values of the original dm are returned
 */
vector<vector<float>> depthMapToPartlyDepthMap(vector<vector<float>> &dm, float dmResolution, float originalXStartpos, float originalYStartpos, float xStartpos, float yStartpos, int xSize, int ySize,
		float minZ){

	vector<vector<float>> partlyDm;
	for (int i = 0; i < xSize; i++){
		partlyDm.push_back(vector<float>());
		for (int j = 0; j < ySize; j++){
			partlyDm.back().push_back(minZ);
		}
	}

	int originalXSize = dm.size();
	int originalYSize = dm[0].size();

	// indices for the original map
	int xStartIndex = (xStartpos-originalXStartpos)/dmResolution;
	int yStartIndex = (yStartpos-originalYStartpos)/dmResolution;
	// for all indices in partlydm and original dm copy values from original dm to partly dm
	for (int i = max(xStartIndex,0), end=min(originalXSize-xStartIndex, xSize)+xStartIndex; i < end; i++){
		for (int j = max(yStartIndex,0), end=min(originalYSize-yStartIndex, ySize)+yStartIndex; j < end; j++){
			partlyDm[i-xStartIndex][j-yStartIndex] = dm[i][j];
		}
	}

	return partlyDm;
}


/**
 * This method is bs for now * maybe sth wrong cause darn slow
 */
vector<vector<float>> getZFromBbx(point3d bbxmin, point3d bbxmax, OcTree* octree, float stepSize){
	vector<vector<float>> retVec;
	double min[3];
	octree->getMetricMin(min[0], min[1], min[2]);

	for (float i = bbxmin.x(); i < bbxmax.x(); i+= stepSize){
		retVec.push_back(vector<float>());
	}

	for (float i = bbxmin.x(); i < bbxmax.x(); i+= stepSize){
		for (float j = bbxmin.y(); i < bbxmax.y(); j+= stepSize){
			retVec[i].push_back(min[0]);
		}
	}
	cout << "bbxretvec created " << endl;
	//would need to rewrite to get omp speedup
	for (OcTree::leaf_bbx_iterator iter = octree->begin_leafs_bbx(bbxmin, bbxmax), end = octree->end_leafs_bbx();
			iter != end; iter++){
			int xIndex = int((iter.getX()-bbxmin.x())/stepSize);
			int yIndex = int((iter.getY()-bbxmin.y())/stepSize);
			float zHeight = (iter.getZ());
			if(zHeight > retVec[xIndex][yIndex]){
				retVec[xIndex][yIndex] = zHeight;
		}
	}
	cout << "leaf iteration finished" << endl;
	return retVec;
}

//----------------------------------------------------------------------------------------
// print functions
//----------------------------------------------------------------------------------------

/**
 * decprecated will be removed
 */
int printHello(void){
	cout << "hi" << endl;
	return 0;
}

//----------------------------------------------------------------------------------------
// open functions
//----------------------------------------------------------------------------------------

/**
 * This method is used to open the file dialogue and get the filepath of the selected file
 * input:
 * 		void
 * return:
 * 		type	: string with complete filepath
 * 		usage	: path to file
 */
string openFile(void){
	char tempFile[1024];
	FILE *f = popen("cd ~/Documents/Models && zenity --file-selection", "r");
	fgets(tempFile, 1024, f);
	string file = string(tempFile);
	return file;
}

/**
 * This method is used to open the file dialogue and retrieve the filepath, with selected file,
 * but without the file extension.
 * input:
 * 		void
 * return:
 * 		type	: string with filepath but without extension
 * 		usage	: get filepath, add extension manual
 */
string openFileRawName(void){
	string filePath = openFile();
	cout << "filpat: " << filePath << " len: " << filePath.length()  << endl;
	// remove fileextension
	string tst(filePath);
	size_t lastindex = tst.find_last_of(".");
	string rawname = tst.substr(0, lastindex);
	return rawname;
}

/**
 * This method is used to open bt files with octovis
 * ! there is no file check !
 * therefore you are self dependent, that it works =)
 * input:
 * 		path: path to file
 * return:
 * 		Null
 */
void openFileInOctovis(string path){
	stringstream ss;
	ss << "octovis " << path;
	string sysCmd = ss.str();
	cout << sysCmd << endl;
	string octovisReturn = exec(sysCmd.c_str());
	cout << octovisReturn << endl;
	sysCmd.clear();
	ss.clear();
}

//----------------------------------------------------------------------------------------
// file preprocessing
//----------------------------------------------------------------------------------------

/**
 * This method is used to create a binvox file from a 3D-Model.
 * This method uses binvox by Patrick Min.
 * Link: http://www.patrickmin.com/binvox
 * The binvox script must lie in /usr/local/bin or other path from $path
 * input:
 * 		path		: string with full path to file, including file extension
 * 		resolution	: resolution of binvox file, up to 1024 or 4096 (binvox) (in voxel per dim)
 * 		exact		: 1 if binvox -e should be enabled
 * 					  2 if not. (-e exact calculation, no GPU)
 * return:
 * 		string		: new file name
 */
string voxelizeFile(string path, int resolution, int exact){
	stringstream ss;

	if (exact == 1){
		ss << "binvox " << "-e " << "-d " << resolution << " " << path;
	} else {
		ss << "binvox " << "-d " << resolution << " " << path;
	}
	string sysCmd = ss.str();
//	system(sysCmd.c_str());
	string binvoxReturn = exec(sysCmd.c_str());
	cout << binvoxReturn << endl;
	sysCmd.clear();
	ss.clear();

	// extract new filename
	size_t lastindex = binvoxReturn.find("write_file(");
	string secondPart = binvoxReturn.substr(lastindex+11);
	lastindex = secondPart.find(")");
	string newFileName = secondPart.substr(0,lastindex);

	return newFileName;
}

/**
 * This method is used to create a bonsai tree(octree) from a binvox file.
 * This method uses binvox2bt from octomap.
 * Octomap: From Armin Hornung and Kai M. Wurm and Maren Bennewitz and Cyrill  Stachniss and Wolfram Burgard.
 * Link: https://octomap.github.io/
 * input:
 * 		path	: string with complete filepath
 * return:
 * 		string	: new file name
 */
string createBtFromBinvoxfile(string path){
	stringstream ss;
	ss << "binvox2bt " << path;
	string sysCmd = ss.str();
	string binvox2btReturn = exec(sysCmd.c_str());
	cout << binvox2btReturn;
	sysCmd.clear();
	ss.clear();

	//extract new filename
	size_t lastindex = binvox2btReturn.find("Writing octree to ");
	string secondPart = binvox2btReturn.substr(lastindex+18);
	lastindex = secondPart.find(".bt");
	string newFilePath = secondPart.substr(0,lastindex+3);
	lastindex = newFilePath.find_last_of("/");
	string newFileName = newFilePath.substr(lastindex+1);
	cout << "\t new path: " << newFilePath << " new file name " << newFileName << "\n" << endl;

	return newFileName;
}

/**
 * This method is used to retrive a path to save a file via the gui
 * return:
 * 		string		:	filepath without fileextension, but name
 */
string getOutputFilePath(void){
	char tempFile[1024];
	FILE *f = popen("cd ~/Documents/Models && zenity --file-selection --save", "r");
	fgets(tempFile, 1024, f);
	string file = string(tempFile);
	return file;
}

/**
 * !!! deprecated not necessary key ray is vector -> deep copy exists much easier!!!
 * This method is used to create a deep copy of a KeyRay object
 * input:
 * 		KeyRay		:	source object
 * return:
 * 		KeyRay		:	copie of input file
 */
KeyRay deepCopyKeyRay(KeyRay rayToCp){
	KeyRay kr;
	KeyRay::iterator key = rayToCp.begin();
	for (int i = 0; i < rayToCp.size(); i++){
		kr.addKey(*(key));
		key ++;
	}
	return kr;
}


/**
 * https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
 * Run system command and get return
 * input:
 * 		const char*	:	command to execute
 * return:
 * 		string		:	system output
 */
std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

/*
 * https://stackoverflow.com/questions/4839626/element-count-of-an-array-in-c
 * This method is used to create a human readable string from a double string of arbitrary length
 */
template<size_t N>
string doubleArrToString(double (&d)[N]){
	int arrLength = N;
	stringstream result;
	result << "[ ";
	for (int i = 0; i < arrLength; i++){
		result << d[i] << " ";
	}
	result << " ]";
	return result.str();
}

//template<size_t K>
template<size_t N, size_t M>
vector<vector<float>> array2DToVector(float (&d)[N][M]){
	cout << M << endl;
	vector<vector<float>> retVec;
	for (long i = 0; i < N; i++){
		retVec.push_back(vector<float>());
	}
	for (long i = 0; i < N; i++){
		for (long j = 0; j < M; j++){
			retVec[i].push_back(d[i][j]);
		}
	}
	return retVec;
}


vector<vector<float>> multiArray2DToVector(boost::multi_array<float, 2> d){
	//boost::multi_array<float, 2>::size_type M = d.shape()[0];
	//boost::multi_array<float, 2>::size_type N = d.shape()[1];
	auto N = d.shape()[0];
	auto M = d.shape()[1];


	//typedef fMultiArr2D::index index;

	vector<vector<float>> retVec;

	for (size_t i = 0; i < N; i++){
		retVec.push_back(vector<float>());
	}
	for (size_t i = 0; i < N; i++){
		for (size_t j = 0; j < M; j++){
			retVec[i].push_back(d[i][j]);
		}
	}

	return retVec;
}

void printOtSpecs(octomap::OcTree* ot){
	double otMin[3];
	double otMax[3];
	ot->getMetricMin(otMin[0],otMin[1],otMin[2]);
	ot->getMetricMax(otMax[0],otMax[1],otMax[2]);
	cout << "maintree is size : " << ot->calcNumNodes() << " max depth: " << ot->getTreeDepth() << endl;
	cout << "maintree resolution is: "  << ot->getResolution() << endl;
	cout << "metric dims are: [ " << otMin[0] << " , " << otMin[1] << " , " << otMin[2] << " ] :: [ "
			<< otMax[0] << " , " << otMax[1] << " , " << otMax[2] << " ]" << endl;
}

/**
 * This Method prints oc specs
 */
void compareOtSpecs(octomap::OcTree *mainOt, octomap::OcTree *machinedOt){
	// compare ot outputs
	cout << "maintree is size : " << mainOt->calcNumNodes() << " max depth: " << mainOt->getTreeDepth() << endl;
	cout << "maintree resolution is: "  << mainOt->getResolution() << endl;
	cout << "machined octree size : "   << machinedOt->calcNumNodes() << "  max depth: " << machinedOt->getTreeDepth() << endl;
	cout << "machined octree res is: "  << machinedOt->getResolution() << endl;
}


void saveOcTree(OcTree* ot, bool saveAsBinary){
	string otPath = getOutputFilePath();
	otPath = otPath.substr(0,otPath.size()-1);
	saveOcTree(ot, saveAsBinary, otPath);
}


void saveOcTree(OcTree* ot, bool saveAsBinary, string path){
	ot->updateInnerOccupancy();
	string otPath = path;
	// only save if a filepath is entered
	if (otPath == "") return;
	if (saveAsBinary) {
		// all files shall end with bt / ot
		if (0 != otPath.compare(otPath.length() - 3, 3, ".bt")){
			otPath += ".bt";
		}
		ot->writeBinary(otPath);
	} else {
		// all files shall end with bt / ot
		if (0 != otPath.compare(otPath.length() - 3, 3, ".ot")){
			otPath += ".ot";
		}
		ot->write(otPath);
	}
	cout << "otm written to: " << otPath << endl;
}


ifstream readGCodeFile(string path){
	ifstream myFile(path);
//	string line;
//	int lineNumber = 0;
	if (myFile.is_open()){
//		while (getline (myFile, line) ){
//			cout << "l[" << lineNumber << "]" << line << '\n';
//			lineNumber ++;
//		}
		return myFile;
	}
	throw new exception;
}

point3d getNextPointFromScaledGCode(ifstream &file){
	string line;
	if (getline (file, line) ){
		return getPointFromScaledGCodeLine(line);
	} else {
//		file.close();
//		throw new exception;
		return point3d();
	}
}

//TODO insert same code as in python.
// create second method, for scaled gcode, where points are "directly" read -> no checks needed
point3d getPointFromGCodeLine(string line, point3d oldPoint){
	point3d pt;
	//https://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
	transform(line.begin(), line.end(), line.begin(), ::tolower);
	int xpos = line.find("x");
	int ypos = line.find("y");
	int zpos = line.find("z");
	return pt;
}

point3d getPointFromScaledGCodeLine(string line){
	int xpos = line.find("x");
	int ypos = line.find("y");
	int zpos = line.find("z");
	point3d pt ( stof(line.substr(xpos+1,ypos-1)), stof(line.substr(ypos+1,zpos-1)), stof(line.substr(zpos+1)));
	return pt;
}


/**
 * This method is used to calculate the euclidean distance between 2 points.
 */
float euclideanDistance(point3d pt0, point3d pt1){
	return (float)sqrt(pow((pt0.x()-pt1.x()),2) + pow((pt0.y()-pt1.y()),2) + pow((pt0.z()-pt1.z()),2));
}


/**
 * This method is used to retrieve a list of all files contained in a directory.
 * input:
 * 		string			: path: directory from which to retrieve the file entries
 * output:
 * 		vector<string>	: vector containing absolute file paths in string format of all files
 * 							in path
 */
vector<string> retrieveFilesListFromPath(string path){
	vector<string> filesList;
	for (auto entry : filesystem::directory_iterator(path)){
		if (PRINT_FOUND_FILES)
		cout << entry.path() << endl;
		filesList.push_back(entry.path().u8string());
	}
	return filesList;
}


/**
 * This method is used to retrieve a list of all files contained in a directory.
 * input:
 * 		string			: path: directory from which to retrieve the file entries
 * output:
 * 		vector<string>	: vector containing absolute file paths in string format of all files
 * 							in path with file_extension==fileExtension
 */
vector<string> retrieveFilesListFromPath(string path, string fileExtension){
	vector<string> filesList;
	for (auto entry : filesystem::directory_iterator(path)){
		string fileString = entry.path().u8string();
		int i = fileString.find_last_of(".");
		if ( (i <= fileString.length()) && (i > 0) ){
			if (fileString.substr(fileString.find_last_of("."))==fileExtension){
				if (PRINT_FOUND_FILES)
				cout << entry.path() << endl;
				filesList.push_back(entry.path().u8string());
			}
		}
	}
	return filesList;
}


string getFileNameFromPath(string path){
	return path.substr(path.find_last_of("/")+1);
}


string getFileNameFromPathWithoutFileExtension(string path){
	int substrBegin = path.find_last_of("/")+1;
	return path.substr(substrBegin, path.find_last_of(".")-substrBegin);
}


string getFileNameFromPathWithoutFileExtension(string path, int numberOfExtensions){
	string filename = path;

	try {
		for (int i = 0; i < numberOfExtensions; i++){
			filename = getFileNameFromPathWithoutFileExtension(filename);
		}
		return filename;
	} catch (const exception& e){
		return filename;
	}
}


/**
 * Method to write states / depth maps to files.
 * Values are separated 							by SEPARATOR .
 * Lines are separated 								by LINE_SEPARATOR .
 * States / depth maps / 2D-arrays are sepparated 	by STATE_SEPARATOR .
 * inputs:
 * 		vector<vector<float>>&	: array2D : depth map or other 2D-matrix with
 * 											vector<vector<float>>& format
 * 		ofstream&				: file : file which to write the values from array2D
 */
template <class T>
void write2DVectorToFile(vector<vector<T>>& array2D, ofstream& file){
	for (vector<T> vec : array2D){
		for (int i = 0; i < vec.size() - 1; i ++){
			// appen state vars and (value) separator
			file << to_string(vec[i]) << SEPARATOR;
		}
		file << to_string(vec.back()) << LINE_SEPARATOR;
		// last var of line gets a line separator
	}
	file << STATE_SEPARATOR;
}

void plot2DVectorWithGnuPlot(vector<vector<float>> map){
	//TODO find scaling factor
	Gnuplot gp;
	gp << "unset key\n";
	gp << "set pm3d\n";
	gp << "set hidden3d\n";
	gp << "set view map\n";
//	gp << "set xrange [ -20 : 270 ] \n";
//	gp << "set yrange [ -20 : 270 ] \n";
	gp << "set xrange [ -10 : " + to_string(map.size()+10) + " ] \n";
	gp << "set yrange [ -10 : " + to_string(map[0].size()+10) + " ] \n";
//	gp << "set zrange [ 14 : 17 ] \n";
	//gp << "set zrange [  -10 : -9.8  ] \n";
	gp << "splot '-'\n";
	gp.send2d(map);
	gp.flush();
	cout << "gnuplot ended?" << endl;
}


//http://www.codebind.com/cpp-tutorial/c-get-current-directory-linuxwindows/
std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  getcwd( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

//TODO es ex compare_octrees funktion, vll diese für bearbeitung fertig verwenden


