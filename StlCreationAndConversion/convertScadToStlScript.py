#this script is windows only
#currently only works for files without blank spaces

from subprocess import call
import subprocess
import os
import datetime

#read current dir (dir in whitch this file lies)
line = os.popen("dir").read()

#split at new lines
lns = line.split("\n")

#get current path
tempSplit = lns[3].split(" ")
currentDir = tempSplit[len(tempSplit)-1]
print("Dir: ", currentDir)

#create list for files with right extension
allFiles = []
#find all stl files
stlFilesRaw = []

#loop over all cmd lines
for i in lns:
    #find files with scad or txt file extension
    if len(i) > 0:
        if ((i.lower().find("scad")>len(i)-5) | (i.lower().find("txt")>len(i)-4)):
            print(i)
            allFiles.append(i)
        elif(i.lower().find("stl")>len(i)-4):
            stlFilesRaw.append(i)

print("")
filteredFiles = []
#postprocess files to extract filenames
for dat in allFiles:
    tempDat = dat.split(" ")
    filteredFiles.append(tempDat[len(tempDat)-1])

filteredStlFiles = []
#extract filenames of stl dats
for dat in stlFilesRaw:
    filteredStlFiles.append(dat.split(" ")[len(dat.split(" "))-1])

filesToConvert = []
tempDat = ""
filesToConvWithoutExt = []
for dat in filteredFiles:
    findCount = 0
    for stDat in filteredStlFiles:
        #delete fileextension
        if dat.lower().find(".txt") >= len(dat)-5:
            tempDat = dat[:len(dat)-4]
        elif dat.lower().find(".scad") >= len(dat)-6:
            tempDat = dat[:len(dat)-5]
        #print("objs(txt,stl): ", tempDat, " ", stDat)
        if tempDat == stDat[:len(stDat)-4]:
            findCount += 1
            break
    if findCount == 0:
        filesToConvert.append(dat)
        filesToConvWithoutExt.append(tempDat)

#print("filteredFiles ",filteredFiles)
#print("filesToConv ",filesToConvert)

utcnow = datetime.datetime.utcnow()
now = datetime.datetime.now()
log = []
#utc time
log.append("\nUTC-Date: " + str(utcnow.second) + ":" + str(utcnow.minute) + ":" \
           + str(utcnow.hour) + " :: " + str(utcnow.day) + "/" + str(utcnow.month) \
           + "/" + str(utcnow.year) + "  #ss:mm:hh dd/mm/yyyy\n")
#local time
log.append("Local-Date: " + str(now.second) + ":" + str(now.minute) + ":" \
           + str(now.hour) + " :: " + str(now.day) + "/" + str(now.month) \
           + "/" + str(now.year) + "  #ss:mm:hh dd/mm/yyyy\n\n")

#convert files to stl
for i in range(len(filesToConvert)):
    dat = filesToConvert[i]
    datWExt = filesToConvWithoutExt[i]
    cmdStr = "openscad -o " + datWExt + ".stl ./" + dat
    log.append(cmdStr+"\n")
    #log.append("createing: "+ datWExt+".stl")
    print("cmdStr: ", cmdStr)
    line = os.popen(cmdStr).read()
    for ln in line:
        log.append(ln)
        if ln == "]":
            print("creatd: ", datWExt+".stl")

#create log file
datFound = 0
for dat in lns:
    if dat.lower().find("logscadtostlconverter.log") > 0:
        datFound += 1
  #      print("foundDat")

if datFound > 0:
    f = open('logScadToStlConverter.log', 'a')
    print("append to log")
else:
    f = open('logScadToStlConverter.log', 'w')
    print("cerate new log")
        
#for ln in log:
f.writelines(log)
if len(log) <= 2:
    f.write("nothing done\n\n")
f.flush()
f.close()

print("python script finished")
