#this script is windows only
#currently only works for files without blank spaces

from subprocess import call
import subprocess
import os
import datetime
import glob

#read current dir (dir in whitch this file lies)
relativePath = "./test/"
filesToConvert = glob.glob(relativePath+"*.txt") + glob.glob(relativePath+"*.scad")
print("Files found: ", os.getcwd(), " : ", filesToConvert)

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
    if dat.find(".scad") > -1:
        datWExt = dat[0:-5]
    else: # ext == .txt
        datWExt = dat[0:-4]
    cmdStr = "openscad -o " + datWExt + ".stl ./" + dat
    log.append(cmdStr+"\n")
    #log.append("createing: "+ datWExt+".stl")
    print("cmdStr: ", cmdStr)
    line = os.popen(cmdStr).read()
    log.append("created: " + datWExt)
    print("creatd: " + datWExt + ".stl")

#create log file
datFound = 0
lns = glob.glob("*.log")
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
