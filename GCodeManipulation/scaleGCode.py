# this script is used to create scaled gcode
# used gcode may only use g91 and g01

# philosophy: print es less as possible for speed (XD but using python :D)

import numpy as np
from os import listdir
from os.path import isfile, join


#-------------------------------------------------------------------------------------------------------------------------------------
#   params
#-------------------------------------------------------------------------------------------------------------------------------------

# path to gcode file
#startFilePath = "C:/Users/*/Desktop/AlienFace.nc"
startFilePath = "C:/path_to_file"
# path for target file (scaled gcode)
modifiedFilePath = "C:/path_to_file/scaledAlienFace.nc"

# current position 3 dimensions
curPos = np.zeros(3)

# max difference from gcode point before creating a new point.
# this value should be bigger than 
# / maximum distance from current position to next gcode point for which points in between get added.
tolerance = .001

# max stepwidth of ann or in this case of the generated gcode
maxStepwidth = .1

# number of positions to create before writing to the target file
posBufferSize = 100


#-------------------------------------------------------------------------------------------------------------------------------------
#   functions
#-------------------------------------------------------------------------------------------------------------------------------------

# opens file from path and creates file from path if no file with name and
#   path exists.
#   returns [startFile, targetFile]
def openAndCreateFile(startFilePath, modifiedFilePath):
    try:
        # open gcode file for reading
        startFile = open(startFilePath,'r')
        # create target file, if already exists -> failing
        #targetFile = open(modifiedFilePath,'x')
        targetFile = open(modifiedFilePath,'w')#TODO
        return [startFile, targetFile]
    except:
        print("file creation of: " + modifiedFilePath + " \nfor: " + startFilePath + " failed!")

# https://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float
def is_number_tryexcept(s):
    """ Returns True is string is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False

def createTargetPosFromGCodeLine(line, curPos):
    line = line.lower() # change everything to lowercase
    line = line.replace('\n','') # remove linebreak
    targetPos = curPos.copy()

    #print("line: ", line, " curPos: ", curPos)
    
    # if x exists in line extract x value
    if line.find('x') > -1:
        #found x value
        tempLine = line[line.find('x')+1:]
        #print("tempLine: ", tempLine)
        spaceIdx = tempLine.find(' ')
        yIdx = tempLine.find('y')
        zIdx = tempLine.find('z')
        # if not found set to complete length
        if spaceIdx < 0:
            spaceIdx = len(tempLine)
        if yIdx < 0:
            yIdx = len(tempLine)
        if zIdx < 0:
            zIdx = len(tempLine)
        xCoord = tempLine[:min(spaceIdx,yIdx,zIdx)]
      #  print("xCoord: ", xCoord, "space,y,z-Idx: ", spaceIdx, yIdx, zIdx)
        if is_number_tryexcept(xCoord):
            targetPos[0] = float(xCoord)

    # if y exists in line extract y value
    if line.find('y') > -1:
        #found y value
        tempLine = line[line.find('y')+1:]
        spaceIdx = tempLine.find(' ')
        xIdx = tempLine.find('x')
        zIdx = tempLine.find('z')
        # if not found set to complete length
        if spaceIdx < 0:
            spaceIdx = len(tempLine)
        if xIdx < 0:
            xIdx = len(tempLine)
        if zIdx < 0:
            zIdx = len(tempLine)
        yCoord = tempLine[:min(spaceIdx,xIdx,zIdx)]
        if is_number_tryexcept(yCoord):
            targetPos[1] = float(yCoord)

    # if z exists in line extract z value
    if line.find('z') > -1:
        #found z value
        tempLine = line[line.find('z')+1:]
        spaceIdx = tempLine.find(' ')
        yIdx = tempLine.find('y')
        xIdx = tempLine.find('x')
        # if not found set to complete length
        if spaceIdx < 0:
            spaceIdx = len(tempLine)
        if yIdx < 0:
            yIdx = len(tempLine)
        if xIdx < 0:
            xIdx = len(tempLine)
        zCoord = tempLine[:min(spaceIdx,yIdx,xIdx)]
        #print("zCoord: ", zCoord, " isNum: ", is_number_tryexcept(zCoord))
        if is_number_tryexcept(zCoord):
            targetPos[2] = float(zCoord)

    #print("createPosRes: ", targetPos)
    return targetPos


def scaleGcode(startFile, targetFile, curPos):
    print("gCode scaling started")
    targetPos = np.zeros(3) # [x,y,z] only here for better understanding
    # buffer before writing
    # -> theoretically double buffered this way
    pointsList = []
    pointsList.append(curPos)
    for p in startFile:
        targetPos = createTargetPosFromGCodeLine(p, curPos)
        # if point is euqal to last point in list do not run code
        if not (targetPos == pointsList[len(pointsList)-1]).min():
         #   pointsList.append(curPos)
            #TODO need to change to numpy
            delta = targetPos - curPos
            #print("targetPos ", targetPos, "curPos ", curPos, "delta: ", delta)
            while np.absolute(delta).max() > tolerance:
                
                # no scaling needed -> append new point directly
                if np.absolute(delta).max() <= maxStepwidth:
                    curPos = curPos + delta
                    pointsList.append(curPos)
                    #print("\n\tKoks\n")
                else:
                    # this value could be stored as desired KNN output, but wont be done
                    # reason is: this can be calculated quickly and i may be doing it somewhere else
                    # maybe just in less efficient method -> not_sure_yet=True
                    delta = delta * (maxStepwidth / np.absolute(delta).max())
                    curPos = curPos + delta
                    pointsList.append(curPos)
                    #print("\n\tFrauenDerNAcht\n")

                delta = targetPos - curPos
            #print("WhileyDone")
                
                    
                
            # TODO overthink if this makes sense
            # set current position to "old" target position to avoid error propagation
            curPos = targetPos.copy()
            #DUNNO here
            pointsList.append(curPos)

            if len(targetPos) > posBufferSize:
                #write to file and clear
                oldPoint = np.array((-10000,-10000,-10000))
                for pos in pointsList:
                    #print("writing to file")
                    # check if point is new point or the same, add only if new point
                    if not (oldPoint == pos).min():
                        targetFile.write("x" + str(pos[0]) + " y" + str(pos[1]) + " z" + str(pos[2]) + "\n")
                    oldPoint = pos
                pointsList.clear()

    oldPoint = np.array((-10000,-10000,-10000))
    for pos in pointsList:
        #print("writing to file")
        #print("\n\n",pos,"\n\n")
        if not (oldPoint == pos).min():
            targetFile.write("x" + str(pos[0]) + " y" + str(pos[1]) + " z" + str(pos[2]) + "\n")
    pointsList.clear()
    
    print("gcode scaling done")
        

#-------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
#   main code
#----------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

# get files list
# src: https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
mypath = startFilePath
print(mypath)
# retrieve all gcode files in path
tempFiles = [f for f in listdir(mypath) if (isfile(join(mypath, f)))]
files = []
for f in tempFiles:
    # only scale files that were not already scaled
    if not (f.find("scScaled") == 0):
        if f.lower().find(".nc") > -1:
            files.append(f)
        elif f.lower().find(".gcode") > -1:
            files.append(f)
        elif f.lower().find(".g") > -1:
            files.append(f)

print(files)


for f in files:
    [startFile, targetFile] = openAndCreateFile( (mypath + f) , (mypath + "scScaled" + f) )
    scaleGcode(startFile, targetFile, curPos)

#import os
#print("path: ", os.getcwd())

if True:
    startFile.close()
    targetFile.flush()
    targetFile.close()
