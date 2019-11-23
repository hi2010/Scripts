#creates scad objects by cutting a box with a cylinder
import os.path
from pathlib import Path
import random as rndm

numberOfFilesToCreate = 100

number = 0

for fileCount in range(numberOfFilesToCreate):
    #definition of value ranges for new objects
    cubeDimensionRange = [10,10,10,40,40,40] #min x-z, max x-z
    #cutDepthRange = [1,30] # later create min of this val and possible
    minWallThickness = 2 # min thickness of bottom
    # for 0's and insanely high values it is equal to not use cut Location
    cutLocRange = [0,0,0,400,400,400] #min x-z, max x-z
    cutRadiusRange = [1,20]
    minCutOverlap = 1

    cubeDimDifferences = [cubeDimensionRange[3]-cubeDimensionRange[0],\
                          cubeDimensionRange[4]-cubeDimensionRange[1],\
                          cubeDimensionRange[5]-cubeDimensionRange[2]]

    #calculate random values in range, make overlap equal on both sides
    # ->prbly will create explanatory graphic later on
#    cubeX = min(max((rndm.random()*(cubeDimensionRange[3]+.5*cubeDimensionRange[0])),\
#                cubeDimensionRange[0]),cubeDimensionRange[3])
#    cubeY = min(max((rndm.random()*(cubeDimensionRange[4]+.5*cubeDimensionRange[1])),\
#                cubeDimensionRange[1]),cubeDimensionRange[4])
#    cubeZ = min(max((rndm.random()*(cubeDimensionRange[5]+.5*cubeDimensionRange[2])),\
#                cubeDimensionRange[2]),cubeDimensionRange[5])

    #create cube dimensions
    cubeX = cubeDimensionRange[0] + rndm.random() * cubeDimDifferences[0]
    cubeY = cubeDimensionRange[1] + rndm.random() * cubeDimDifferences[1]
    cubeZ = cubeDimensionRange[2] + rndm.random() * cubeDimDifferences[2]
    cubeDimensions = [cubeX,cubeY,cubeZ]


    cutRadius = min(max(rndm.random()*(cutRadiusRange[1]+(cutRadiusRange[0]*.5)),cutRadiusRange[0]),\
                    cutRadiusRange[1])
    #cut radius should be smaler than min of length and with of cube else cube just might get shorter
    cutRadius = min(min(cutRadius,cubeDimensions[0]/2),cubeDimensions[1]/2)


    #TODO
    cutLocXMin  = - cutRadius + minCutOverlap
    cutLocXMax  =   cutRadius - minCutOverlap + cubeDimensionRange[0]
    cutLocXDiff = cutLocXMax - cutLocXMin
    cutLocX = cutLocXMin + rndm.random() * cutLocXDiff
    cutLocX = min(max(cutLocRange[0],cutLocX),cutLocRange[3])

    cutLocYMin  = - cutRadius + minCutOverlap
    cutLocYMax  =   cutRadius - minCutOverlap + cubeDimensionRange[1]
    cutLocYDiff = cutLocYMax - cutLocYMin
    cutLocY = cutLocYMin + rndm.random() * cutLocYDiff
    cutLocY = min(max(cutLocRange[1],cutLocY),cutLocRange[4])
    
    #cutLocX = min(max((rndm.random()*(cutLocationRange[3]+.5*cutLocationRange[0])),\
     #                 cutLocationRange[0]),cutLocationRange[3])

    cutLocZMin  = max(minWallThickness, cubeDimensions[2] - cutLocRange[5])
    cutLocZMax  = min(cubeDimensions[2] - cutLocRange[2], cubeDimensions[2]-minCutOverlap)
    cutLocZDiff = cutLocZMax - cutLocZMin
    cutHeight  = cutLocZMin + rndm.random() * cutLocZDiff
    
    #cutDepth = min(max(rndm.random()*(cutDepthRange[1]+(cutDepthRange[0]*.5)),cutDepthRange[0]),\
    #               cutDepthRange[1])
    #cutDepth = min(cutDepth,cubeZ-minWallThickness)
    #cutHeight = cubeDimensions[2]-cutDepth
    
    #TODO
    cutLocation = [cutLocX, cutLocY, cutHeight]
    cutDepth = cubeDimensions[2]-cutLocation[2]+.1



    print("cube dims : ", cubeDimensions)
    print("cut depth : ", cutDepth)
    print("cut height: ", cutHeight)
    print("cut loc   : ", cutLocation)
    print("cut radius: ", cutRadius)


    fileName = "testFile"



    #openscad datei
    output = "cubeDimensions = " + str(cubeDimensions) + \
             ";\ncutDepth = " + str(cutDepth) + \
             ";\ncutHeight = cubeDimensions[2]-cutDepth;" + \
             "\ncutLocation = " + str(cutLocation) + \
             ";\ncutRadius = " + str(cutRadius) + \
             ";\n\ndifference(){" + \
             "\n    cube(cubeDimensions);\n" + \
             "    translate(cutLocation){\n" + \
             "        cylinder(h=cutDepth+.1,r=cutRadius);\n" + \
             "    }" + \
             "}"


    myFile = Path(fileName + str(number) + ".scad")
    while myFile.is_file():
        if myFile.is_file():
            print("file exists, filename", str(myFile.absolute()))
        else:
            print("create new file, filename", str(myFile.absolute()))
            break
        number += 1
        myFile = Path(fileName + str(number) + ".scad")

    f = open(myFile.absolute(),'w')
    f.writelines(output)
    print("creaetd:  " + str(myFile.absolute()), " (file {:d} of {:d})".format(fileCount+1,numberOfFilesToCreate), "\n")
    f.flush()
    f.close()
