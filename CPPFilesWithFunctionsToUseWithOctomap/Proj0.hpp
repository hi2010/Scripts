/*
 * Proj0.hpp
 *
 *  Created on: Dec 11, 2018
 *      Author: louis
 */

#ifndef PROJ0_HPP_
#define PROJ0_HPP_


#include <string>
#include "octomap/octomap.h"
#include "boost/multi_array.hpp"

using namespace std;

// 2d boundingbox container
struct bbx2D{
	float xMin;
	float yMin;
	float xMax;
	float yMax;
};

void testWorkflow(void);
void testWorkflow2(void);
void testWorkflow3(void);
void testWorkflow4(void);
void testWorkflow5(void);
void testWorkflow6(void);
void testWorkflow7(void);
void testWorkflow8(void);
void modelToOctovis(void);

void createTrainingdata(string gCodePath, string btPath, string machinedBtPath,
		string stateFilePath, string targetStateFilePath);

string stlToOt(void);
string stlToOt(string path);
string stlToBt(void);
string stlToBt(string path);
string btToOt(string path);
int stlToOctovis(void);
int stlToOctovis(string path);

//octomap::OcTree openBtFileWithOctomap(string filePath);
int openOtFileWithOctomap(void); //todo return OcTree*

void addKeysFromRay(octomap::KeyRay, octomap::Pointcloud&, octomap::OcTree*);
octomap::OcTree vector2dToOctree(vector<vector<float>> heightMap, float resolution, octomap::point3d org);
void addVector2dToOctree(vector<vector<float>> heightMap, octomap::OcTree* ot, octomap::point3d org);

octomap::Pointcloud createCircularPointCloud(float radius, float resolution);
octomap::OcTree createCircularOctree(float radius, float resolution,octomap::point3d center);

void insertTool(octomap::OcTree*, octomap::Pointcloud, octomap::pose6d);
void insertToolInLine(octomap::OcTree* ot, octomap::Pointcloud pc, octomap::point3d startPosition, octomap::point3d endPosition);
void insertToolInLineRoundUp(octomap::OcTree* ot, octomap::Pointcloud pc, octomap::point3d startPosition, octomap::point3d endPosition);
void insertToolInLineRoundUpMiniSteps(octomap::OcTree* ot, octomap::Pointcloud pc, octomap::point3d startPosition, octomap::point3d endPosition);
void insertToolInLineUsingOctomapComputeRayKeys(octomap::OcTree* ot, octomap::Pointcloud pc, octomap::point3d startPosition, octomap::point3d endPosition);

float getMaxZFromRay(octomap::KeyRay, octomap::OcTree*);
float getMaxZFromRayWithAlternateMin(octomap::KeyRay, octomap::OcTree*, double);

float getMinZFromRay(octomap::KeyRay, octomap::OcTree*);
float getMinZFromRayOptimised(octomap::KeyRay, octomap::OcTree*);

vector<vector<float>> areaToHeightmap(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize);
vector<vector<float>> areaToHeightmapAlternateMinZ(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize, double minZ);
vector<vector<float>> areaToHeightmapRelativeToBody(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize);
// searches for min z val; height map for max
vector<vector<float>> areaToDepthmap(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize);
vector<vector<float>> areaToDepthmapWithBBX(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float bbxMinX, float bbxMinY, float bbxMaxX, float bbxMaxY);
vector<vector<float>> areaToDepthmapWithBBXOptimised(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float inBBXZ, bbx2D bbx);

vector<vector<float>> areaToHeightmapWithLeafIterator(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ);
vector<vector<float>> areaToHeightmapWithLeafIteratorAndBBX(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float inBBXZ, bbx2D bbx);
vector<vector<float>> areaToDepthmapWithLeafIteratorAndBBX(octomap::OcTree* octree, float xStartpos, float yStartpos, int xSize, int ySize, float stepSize
		,float minZ, float inBBXZ, bbx2D bbx);

// decisions for dimension are made when creating the original dm
vector<vector<float>> depthMapToPartlyDepthMap(vector<vector<float>> &dm, float dmResolution, float originalxStartpos, float originalyStartpos, float xStartpos, float yStartpos, int xSize, int ySize,
		float minZ);

vector<vector<float>> getZFromBbx(octomap::point3d bbxmin, octomap::point3d bbxmax, octomap::OcTree* octree, float stepSize);

int printHello(void);

string openFile(void);
string openFileRawName(void);
void openFileInOctovis(string path);

string voxelizeFile(string path, int resolution, int exact);
string createBtFromBinvoxfile(string path);

string getOutputFilePath(void);

std::string exec(const char* cmd);
template<size_t N>
string doubleArrToString(double (&d)[N]);

template<size_t N, size_t M>
vector<vector<float>> array2DToVector(float (&d)[N][M]);
vector<vector<float>> multiArray2DToVector(boost::multi_array<float, 2> d);

void printOtSpecs(octomap::OcTree*);
void compareOtSpecs(octomap::OcTree*, octomap::OcTree*);

void saveOcTree(octomap::OcTree*, bool saveAsBinary);
void saveOcTree(octomap::OcTree*, bool saveAsBinary, string path);

// TODO lossy compression with 2d dm

ifstream readGCodeFile(string);

octomap::point3d getNextPointFromScaledGCode(ifstream&);

octomap::point3d getPointFromGCodeLine(string,  octomap::point3d);
octomap::point3d getPointFromScaledGCodeLine(string);

float euclideanDistance(octomap::point3d, octomap::point3d);

vector<string> retrieveFilesListFromPath(string path);
vector<string> retrieveFilesListFromPath(string path, string fileExtension);

string getFileNameFromPath(string path);
string getFileNameFromPathWithoutFileExtension(string path);
string getFileNameFromPathWithoutFileExtension(string path, int numberOfExtensions);

template <class T>
void write2DVectorToFile(vector<vector<T>>& array2D, ofstream& file);

void plot2DVectorWithGnuPlot(vector<vector<float>>);

string GetCurrentWorkingDir( void );

#endif /* PROJ0_HPP_ */
