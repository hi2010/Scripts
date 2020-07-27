fileP = r"PATH_TO_SOURCE_GCODE_FILE"

[xMin, xMax, yMin, yMax] = [0, 300, 90, 300]
safeHeight = "1"

lastWasUp = False


# line: G-Code line
# direction: char to look for f.e. "X"
# if no val is found None is returned
def getCoordinateFromLine(line, direction):
    if line.upper().find(direction.upper()) < 0:
        return None

    tempLine = line[line.upper().find(direction.upper())+1::]
    xDigit = ""
    
    for c in tempLine:
        # may not work everytime
        if not c.isalpha():
            xDigit += c
        else:
            break
        
    try:
        xVal = float(xDigit)
    except:
        xVal = 0

    return xVal

# cramer rule https://stackoverrun.com/de/q/5648151 @ rook
def findLineIntersection(pp1, pp2, pp3, pp4):
    def line(p1, p2):
        A = p1[1] - p2[1]
        B = p2[0] - p1[0]
        C = p1[0] * p2[1] - p2[0] * p1[1]
        return A, B, -C

    def intersection(L1, L2):
        D = L1[0] * L2[1] - L1[1] * L2[0] 
        Dx = L1[2] * L2[1] - L1[1] * L2[2] 
        Dy = L1[0] * L2[2] - L1[2] * L2[0] 
        if D != 0:
            x = Dx/D 
            y = Dy/D 
            return x,y 
        else: 
             return False
    #print(pp1, pp2, pp3 ,pp4)
    return intersection(line(pp1,pp2), line(pp3, pp4))

def findIntersectionWithBorder(p1, p2):
    #TODO: if both points are outside one needs to get two intersects
    
    # calc slope
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    # no slope = no intersect
    if dx == 0 and dy == 0:
        return False
    # find all intersects
    l1 = findLineIntersection(p1, p2, [xMin, yMin], [xMax, yMin]) # ->
    l2 = findLineIntersection(p1, p2, [xMax, yMin], [xMax, yMax]) # ^
    l4 = findLineIntersection(p1, p2, [xMin, yMax], [xMin, yMin]) # v
    l3 = findLineIntersection(p1, p2, [xMax, yMax], [xMin, yMax]) # <-
    intersects = []
    intersects.append(l1)
    intersects.append(l3)
    intersects.append(l2)
    intersects.append(l4)

    ffCrit = False

    tol = 0.0001
    if xMin-tol <= p1[0] <= xMax+tol and yMin-tol <= p1[1] <= yMax+tol:
        if xMin+tol > p2[0] or p2[0] > xMax-tol or yMin+tol > p2[1] or p2[1] > yMax-tol:
            print("in to out")
            print(p1,p2)
            print(intersects)
            ffCrit = True

##    if xMin+tol > p1[0] or p1[0] > xMax-tol or yMin+tol < p1[1] or p1[1] > yMax-tol:
##        if xMin-tol <= p2[0] <= xMax+tol and yMin-tol <= p2[1] <= yMax+tol:
##            print("out to in")
##            print(p1,p2)
##            print(intersects)
##
##    if xMin+tol > p1[0] or p1[0] > xMax-tol or yMin+tol < p1[1] or p1[1] > yMax-tol:
##        if xMin+tol > p2[0] or p2[0] > xMax-tol or yMin+tol > p2[1] or p2[1] > yMax-tol:
##            print("out to out")
##            print(p1,p2)
##            print(intersects)
    
    # find all intersects with same direction
    #
    # diff of its and point is multiple of diff of points (inner line)
    # so if positive multiple -> same dir, else other
    realIts = []
    for its in intersects:
        if its is False:
            continue
        if abs(dx) > 0:
            if ((its[0] - p1[0]) / dx) >= 0:
                realIts.append(its)
            continue
        # always uneq 0 if dx is 0
        if ((its[1] - p1[1]) / dy) >= 0:
            realIts.append(its)

    if len(realIts) <= 0:
        print("no real")
        print(intersects)
        print(p1,p2)
        print("\n")
        return False
    
    # 3493
    # find its with min distance and whose location is within the boundaries
    resIts = None
    minDiff = None
    tol = .0001
    for its in realIts:
        # calc its with min dist
        newDiff = abs(p1[0] - its[0]) + abs(p1[1] - its[1])
        if minDiff is None:
            if xMin-tol <= its[0] <= xMax+tol and yMin-tol <= its[1] <= yMax+tol:
                minDiff = newDiff
                resIts = its
        elif newDiff < minDiff:
            # check if point is in
            if xMin-tol <= its[0] <= xMax+tol and yMin-tol <= its[1] <= yMax+tol:
                minDiff = newDiff
                resIts = its

    if ffCrit:
        print("rit")
        print(realIts)
        print(resIts)

    if resIts is None:
        print("\n")
        print(p1,p2)
        print(realIts)
        print(intersects)
        print("\n")
        return False
    
    return resIts

[xCoord, yCoord, zCoord] = [0, 0, safeHeight]
[lx, ly, lz] = [xCoord, yCoord, zCoord]
nf = lf = newFeed = 100
lastWasInRange = False
with open(fileP, mode='r') as readFile:
    with open(fileP[0:-3]+"edit.nc", mode='w') as writeFile:
        for line in readFile:
            xInRange = False
            yInRange = False

            nx = getCoordinateFromLine(line, "X")
            ny = getCoordinateFromLine(line, "Y")
            nz = getCoordinateFromLine(line, "Z")
            if nx is not None:
                xCoord = nx
            if ny is not None:
                yCoord = ny
            if nz is not None:
                zCoord = nz

            nf = getCoordinateFromLine(line, "F")
            if nf is not None:
                newFeed = nf

            if xMin < xCoord < xMax:
                xInRange = True
            else:
                xInRange = False

            if yMin < yCoord < yMax:
                yInRange = True
            else:
                yInRange = False

            if xInRange and yInRange:
                if lastWasInRange is False:
                    lastWasInRange = True
                    if newFeed != lf:
                        writeFile.write("G01 F" + str(newFeed) + "\n")
                    intersect = findIntersectionWithBorder([lx, ly], [xCoord, yCoord])
                    if intersect is not False:
                        writeFile.write("G01 X" + str(intersect[0]) + " Y" + str(intersect[1])+ ";is1\n")
                        writeFile.write("G01 Z" + str(zCoord)+ ";zl\n")
                writeFile.write(line[:-2]+";jl\n")
                if findIntersectionWithBorder([lx, ly], [xCoord, yCoord]) is False:
                    writeFile.write("G01 Z" + str(zCoord)+ ";zl\n")
            elif lastWasInRange and xCoord != lx and yCoord != ly:
                if newFeed != lf:
                    writeFile.write("G01 F" + str(newFeed)+ "\n")
                intersect = findIntersectionWithBorder([lx, ly], [xCoord, yCoord])
                if intersect is not False:
                    if zCoord != lz:
                        writeFile.write("G01 X" + str(intersect[0]) + " Y" + str(intersect[1]) + "Z" + str(zCoord) + ";is2z\n")
                    else:
                        writeFile.write("G01 X" + str(intersect[0]) + " Y" + str(intersect[1])+ ";is2\n")
                    
                writeFile.write("G01 Z" + safeHeight + ";zs\n")
                lastWasInRange = False
                #writeFile.write("isecDun"+str(gotNo)+"\n")
                    
            lx = xCoord
            ly = yCoord
            lz = zCoord

            lf = newFeed
