import time
import math

import PyKEP as kep
import numpy as np
import matplotlib.pyplot as plt


from collections import defaultdict
from itertools import combinations
import random
from breakup import setmass, breakup, plot_orbits
from operator import itemgetter

CUBE_RES = 10e3

#"""Returns a list of list of debris which are at time t in the same volume of resolution res."""
def cube(planets, res=CUBE_RES):
    cubemap = defaultdict(list)
    for p in planets:
        try:
            key = tuple(map(lambda foo: int(math.floor(foo/CUBE_RES)), p._r))
            cubemap[key].append(p)
        except:
            pass
    res = [foo for foo in cubemap.values() if len(foo) > 1]
    return res


def collision_prob(p1, p2):
    sigma = (p1.radius + p2.radius) ** 2 * math.pi # crossectional collision area
    dU = CUBE_RES ** 3 # volume of the cube
    Vimp = np.linalg.norm(p1._v - p2._v) # rel. velocity
    return  Vimp/dU * sigma


def setradii(planets):
    for p in planets:
        try:
            a = float(satcat[p.name.strip()].radarA)
            p.radius = math.sqrt(a/math.pi)
        except Exception as e:
            # print e
            p.radius = 0.5 # TODO: sample from radii distribution / use mean


def update(planets, ep):
    for p in planets[:]:
        try:
            p._r, p._v = map(np.array, p.eph(ep))
        except:
            # satelite decayed
            planets.remove(p)

# returns list of recent planets for relaunching sequence
def launchSequence(satcat, beginningOfSequence):
    recentSatcat = dict() #contains all the planets from chosen year till present
    recentDebris = dict() #contains all planets (debris) presented as DEB in satcat
    recentPlanets = dict() #contains all planets which are not DEB
    activeRecentObjects = dict() #recent planets which are not decayed
    for i in satcat.items():
        if len(i[1][6].split("-")) != 3:
            continue
        else:
            yearOfLaunch = i[1][6].split("-")[0]
            monthOfLaunch = i[1][6].split("-")[1]
        if len(yearOfLaunch) == 4:
            yearOfLaunchNo = int(yearOfLaunch)
        if yearOfLaunchNo >= beginningOfSequence:
            recentSatcat.update({i[0]:i[1]})
        if yearOfLaunchNo >= beginningOfSequence and 'DEB' in i[1][4]:
            recentDebris.update({i[0]:i[1]})
        if yearOfLaunchNo >= beginningOfSequence and 'DEB' not in i[1][4]:
            recentPlanets.update({i[0]:i[1]})
        #there is no decay date -> still active
        if yearOfLaunchNo >= beginningOfSequence and len(i[1][8].split("-")[0]) != 4 and 'DEB' not in i[1][4]:
            activeRecentObjects.update({i[0]:i[1]})            
    return recentPlanets, activeRecentObjects

# computes the spatial densities for later visualisation
def computeSpatialDensities(planets):
    altStart = 200 # altitude from [km]
    step = 20 # size of bins [km]
    spaceDens = []
    for alt in range(altStart, 2000, step):
        noOfPlanInAlt = 0
        for i in range(1,len(debris)):
            try:
                r, v = debris[i].eph(16*365)
            except:
                continue
            el = kep.ic2par(r, v, debris[i].mu_central_body)
            altitude = (el[0]-6378000)/1000
            if altitude >= alt and altitude < (alt+step):
                noOfPlanInAlt += 1
        vol1 = (alt+6378)**3 * math.pi * 4 / 3
        vol2 = (alt+step+6378)**3 * math.pi * 4 / 3
        totVol = vol2-vol1
        spaceDens.append(noOfPlanInAlt/totVol)
    return spaceDens

# finds recent active assets and returns a list of them
# recent active assets --> launched at most 10 years ago (life span of 10 years)
                    #  --> not debris
                    #  --> in LEO (its perigee is in range 160-2000km)
def findRecentActiveAssets(planets, currentYear, horizon):
    recentActivePlan = dict()
    for i in planets.items():
        # we assume life span of 10 years
        beginning = currentYear - horizon
        # some satcat entries are missing information - e.g. orbital status code NEA - no elements available
        try:
            apogee = int(i[1][11])
            perigee = int(i[1][12])
            # print apogee, perigee
        except: 
            apogee = 0
            perigee = 0

        semimajorA = ((apogee) + (perigee))/2
        if len(i[1][6].split("-")) != 3:
            continue
        else:
            yearOfLaunchNo = int(i[0][0:4])

        if yearOfLaunchNo >= beginning and len(i[1][8].split("-")[0]) != 4 and 'DEB' not in i[1][4] and perigee > 160 and perigee < 2000:
            # we filter only those which belong to one of the agents (US, Russia, China, EU)
            if (i[1][5].strip() in {"US","CIS","PRC"}) or i[1][5].strip() in EU_list:
                recentActivePlan.update({i[0]:i[1]})

    return recentActivePlan

def decayObjects(planets, currentYear):
    for p in planets:
        yearOfLaunch = int(p.name[0:4])
        # print yearOfLaunch
        decayYear = currentYear - 10
        if yearOfLaunch > 2016 and yearOfLaunch < decayYear:
            planets.remove(p)
            # print "removing"

# def findRecentTLE(tles, satcatPlanets, beginning):
#     recentActiveTLE = []
#     for planet in satcatPlanets.items():
#         for planetTLE in tles:
#             if planet[0].strip() == planetTLE.name.strip():
#                 recentActiveTLE.append(planetTLE)
#                 break
#     return recentActiveTLE
def findRecentTLE(tles, satcatPlanets):
    recentActiveTLE = []
    # beginning = 2015-sequenceLength
    # for satcatPlanet in satcat.items():

    # for planet in tles:
    #     if int(planet.name.strip()[0:4]) > beginning:
    #         recentActiveTLE.append(planet)

    for planet in satcatPlanets.items():
        for planetTLE in tles:
            if planet[0].strip() == planetTLE.name.strip():
                recentActiveTLE.append(planetTLE)
                break
    return recentActiveTLE


def launchRandomAssets(recentActiveTLEObjects, debris, recentActiveAssets, indicesIncl, satcatObjects):
    modifIncl = float(np.random.choice(indicesIncl[1:],1,p=weightsIncl))
    # print len(recentActiveTLEObjects)
    tleObject = random.choice(recentActiveTLEObjects)
    # planetsSATCAT = launchSequence(satcat, beginningOfSequence)
    if modifIncl < 0:
        modifIncl = 0
    if modifIncl > 180:
        modifIncl = 180
    newIncl = '{:8.4f} '.format(modifIncl) # inclination (i)
    # print tleObject.name, modifIncl
    # for planet in satcatObjects.items():
    # if planet[0].strip() == tleObject.name.strip():
    planet = satcatObjects[tleObject.name.strip()]
    # print tleObject.name.strip()
    # print satcatObjects.items()[200][0]
    newLine2 = tleObject.line2[:8] +newIncl+tleObject.line2[17:]

    dayOfYear = (((ep-16*365.25)/365.25+2016) - year)*365
    if dayOfYear < 10:
        epoch = str(year)[2:4] +"00"+ '{:2.8f}'.format(dayOfYear)
    elif dayOfYear >= 10 and dayOfYear < 100:
        epoch = str(year)[2:4] +"0"+ '{:2.8f}'.format(dayOfYear)
    elif dayOfYear >= 100:
        epoch = str(year)[2:4] + '{:2.8f}'.format(dayOfYear)

    newLine1 = tleObject.line1[:18] + epoch + tleObject.line1[32:]

    newObject = kep.planet.tle(newLine1, newLine2)
    newObject.name = str(year) + newObject.name.strip()[4:]
    # print newObject.name
    debris.append(newObject)
    addingNewPlanet = 1
    
    nameOfNew = str(year) + tleObject.name.strip()[4:]
    newSatcat = planet
    # print "adding ",nameOfNew
    # new satcat entries does not have modified inclination because the simulation works with TLE objects
    recentActiveAssets.update({nameOfNew:newSatcat})


def updateCollisionTuples(tupleOfCollisions, removedObjects):
    newListOfCommonRiskPlanets = []
    for pln in tupleOfCollisions:
        if pln[0] in removedObjects:
            tupleOfCollisions.remove(pln)
        elif pln[1] in removedObjects:
            tupleOfCollisions.remove(pln)

    for couple in tupleOfCollisions:
        for planet in newListOfCommonRiskPlanets:
            if couple[0] == planet[0]:
                planet[1] += couple[4]
                break
        else:
            newListOfCommonRiskPlanets.append([couple[0],couple[4]])

        for planet in newListOfCommonRiskPlanets:
            if couple[1] == planet[0]:
                planet[1] += couple[4]
                break
        else:
            newListOfCommonRiskPlanets.append([couple[1],couple[4]])

    return newListOfCommonRiskPlanets, tupleOfCollisions


# for cycle for multiple experiments - output file is marked with number of cycle
for qq in range(1,2):
    satcat = kep.util.read_satcat('satcat.txt')
    debris = kep.util.read_tle('full.tle', with_name=False)
    # list of all EU states and ESA and its bodies
    EU_list = ["ASRA","BEL","CZCH","DEN","ESA","ESRO","EST","EUME","EUTE","FGER","FR","FRIT","GER","GREC","HUN","IT","LTU","LUXE","NETH","NOR","POL","POR","SPN","SWED","SWTZ","UK"]

#  filtering only objects in LEO
    leo_debris = []
    ep = kep.epoch(16.*365.25)
    for p in debris:
        try:
            oe = p.osculating_elements(ep)
            if oe[0] * (1-oe[1]) < 6378000.0 + 2000000:
                leo_debris.append(p)
        except:
            pass
    debris = leo_debris
    # loading tle files of recent objects for relaunching sequence from 2005 to 2014
    # TLE files are downloaded from space-track based on NORAD number range from filtered satcat files
    recentPlanetsTLE = kep.util.read_tle('recentPlanetsTLE_2005_2014.tle', with_name=False)

    # alternative sequence for relaunching objects from 1995 to 2004
    # recentPlanetsTLE = kep.util.read_tle('recentPlanetsTLE_1995_2004.tle', with_name=False)
    beginningOfSequence = 2005

    removeHorizon = 2 #how often we remove objects - 1 -> every year, 10 -> every 10 years
    # satcat entries of recent objects used for relaunching sequence
    recentPlanets, recentActivePlanets = launchSequence(satcat, beginningOfSequence)
    # recentActiveTLE = findRecentTLE(recentPlanetsTLE, recentActivePlanets, beginningOfSequence)
    # print len(recentPlanetsTLE)

    # for computing the inclination distribution and then to sample randomly new ones
    satcatPlanetsLast15years = findRecentActiveAssets(satcat, 2015, 15)
    # print len(satcatPlanetsLast20years)
    recentActiveTLE = findRecentTLE(debris, satcatPlanetsLast15years)
    # print len(recentActiveTLE)


    # inclination is randomly sampled from inclination distribution of objects from (2005-2014)
    inclinations = []
    for planetR in recentActiveTLE:
        inclinations.append(float(planetR.line2.split()[2]))
    valuesIncl,indicesIncl=np.histogram(inclinations,bins=180)
    valuesIncl=valuesIncl.astype(np.float32)
    weightsIncl=valuesIncl/np.sum(valuesIncl)

    setradii(debris)
    setmass(debris)
    setmass(recentPlanetsTLE)

    sump = 0
    decadeYearMultiplier = 1
    numberOfCatastrophicCollisions = 0
    numberOfNonCatastrophicCollisions = 0
    totalNumberOfNewDebris = 0
    totalDebrisEvol = []
    spaceDensities = []

# list of recent active assets used for making "list of important assets" for each of the agents
    recentActiveAssets = findRecentActiveAssets(satcat, 2015, 10)

    yearOfColl = []

# length of experiment in years, standard setting 200 years
    timeHorizon = 150
    # Russia
    removedObjectsPerYear1 = 0

    # USA
    removedObjectsPerYear2 = 1

    # China
    removedObjectsPerYear3 = 0

    # EU
    removedObjectsPerYear4 = 0


    # Russia
    objectsWithPossitiveProbOfCollision1 = []
    objectsProbsOfCollisions1 = []
    
    agent1Collisions = [0] * timeHorizon
    agent1Risks = []
    agent1RisksProb = []
    agent1YearRiskProb = 0
    agent1TotalImportantAssets = []

    # US
    objectsWithPossitiveProbOfCollision2 = []
    objectsProbsOfCollisions2 = []
    agent2Collisions = [0] * timeHorizon
    agent2Risks = []
    agent2RisksProb = []
    agent2YearRiskProb = 0
    agent2TotalImportantAssets = []


    # China
    objectsWithPossitiveProbOfCollision3 = []
    objectsProbsOfCollisions3 = []
    agent3Collisions = [0] * timeHorizon
    agent3Risks = []
    agent3RisksProb = []
    agent3YearRiskProb = 0
    agent3TotalImportantAssets = []

    # EU
    objectsWithPossitiveProbOfCollision4 = []
    objectsProbsOfCollisions4 = []
    agent4Collisions = [0] * timeHorizon
    agent4Risks = []
    agent4RisksProb = []
    agent4YearRiskProb = 0
    agent4TotalImportantAssets = []

    collisionRiskCouplesVirtual = []
    agent1Threats = []
    agent2Threats = []
    agent3Threats = []
    agent4Threats = []

    # list of removed objects and who removed it
    listOfRemovedObjects = []
    listOfRemovedObjectsFinal = []

    objectsWithPossitiveProbOfCollisionAll = []
    objectsProbsOfCollisionsAll = []
    commonRisksProb = []
    commonRisksProbYear = 0

    stringOfRemovals = "every"+str(removeHorizon)+"year_rem"+str(removedObjectsPerYear1) + str(removedObjectsPerYear2) + str(removedObjectsPerYear3) + str(removedObjectsPerYear4)

    virtualRunEnabled = 1
    if removedObjectsPerYear1 == 0 and removedObjectsPerYear2 == 0 and removedObjectsPerYear3 == 0 and removedObjectsPerYear4 == 0:
        virtualRunEnabled = 0

    yearToPlot = 0
    ep = int(math.floor(16*365.25)) 
    virtualRun = 0 #initial setting, if 1 the run is virtual for computing probs of collisions
    # start of removing objects
    yearPom = 2016
    virtualDebris = []
    virtualDebris = debris[:]
    dataSaving = yearPom + 1

    # main loop, one step is 5 days
    while ep < int(math.ceil(16*365.25 + timeHorizon*365.25)) or virtualRun == 1:
        ep += 5
        year = int(math.floor((ep-16*365.25)/365.25)+2016) 
        # print year
        # print (float((ep-16*365.25)/365.25+2015) - year)*365
        # if year > 2020:
        #     start1 = time.time()
        # virtual run
        if virtualRun == 1 and year%removeHorizon == 0 and year != yearPom and virtualRunEnabled:
            ep = ep - int(removeHorizon*365.25)
            virtualRun = 0
            yearPom = year
            year = int(math.floor((ep-16*365.25)/365.25)+2016)


            # we do not remove any important assets
            for asset in agent1Threats:
                if recentActiveAssets.has_key(asset[0]):
                    agent1Threats.remove(asset)
            for asset in agent2Threats:
                if recentActiveAssets.has_key(asset[0]):
                    agent2Threats.remove(asset)
            for asset in agent3Threats:
                if recentActiveAssets.has_key(asset[0]):
                    agent3Threats.remove(asset)
            for asset in agent4Threats:
                if recentActiveAssets.has_key(asset[0]):
                    agent4Threats.remove(asset)

            pomThreats1 = agent1Threats[:]
            pomThreats2 = agent2Threats[:]
            pomThreats3 = agent3Threats[:]
            pomThreats4 = agent4Threats[:]

            # remove objects until n objects are removed
            listOfRemoved = []
            iterator1 = 0
            while iterator1 < removedObjectsPerYear1 and len(agent1Threats) != 0:
                highestThreat = max(agent1Threats, key = itemgetter(1))
                agent1Threats.remove(highestThreat)

                for d in debris:
                    if highestThreat[0] == d.name:
                        print "Russia removing", d.name, highestThreat[1]
                        debris.remove(d)
                        listOfRemovedObjects.append(["CIS","o",d.name,0,0,0,0])
                        listOfRemoved.append(d.name)
                        iterator1 += 1
                        break

            iterator2 = 0
            while iterator2 < removedObjectsPerYear2 and len(agent2Threats) != 0:
                highestThreat = max(agent2Threats, key = itemgetter(1))
                agent2Threats.remove(highestThreat)
                for d in debris:
                    if highestThreat[0] == d.name:
                        print "US removing", d.name, highestThreat[1]
                        debris.remove(d)
                        listOfRemovedObjects.append(["US","o",d.name,0,0,0,0])
                        listOfRemoved.append(d.name)
                        iterator2 += 1
                        break

            iterator3 = 0
            while iterator3 < removedObjectsPerYear3 and len(agent3Threats) != 0:
                highestThreat = max(agent3Threats, key = itemgetter(1))
                agent3Threats.remove(highestThreat)
                for d in debris:
                    if highestThreat[0]== d.name:
                        print "China removing", d.name, highestThreat[1]
                        debris.remove(d)
                        listOfRemovedObjects.append(["PRC","o",d.name,0,0,0,0])
                        listOfRemoved.append(d.name)
                        iterator3 += 1
                        break

            iterator4 = 0           
            while iterator4 < removedObjectsPerYear4 and len(agent4Threats) != 0:
                highestThreat = max(agent4Threats, key = itemgetter(1))
                agent4Threats.remove(highestThreat)
                for d in debris:
                    if highestThreat[0] == d.name:
                        print "EU removing", d.name, highestThreat[1]
                        debris.remove(d)
                        listOfRemovedObjects.append(["EU","o",d.name,0,0,0,0])
                        listOfRemoved.append(d.name)
                        iterator4 += 1
                        break

            listOfCommonRiskPlanets = []

            listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, listOfRemoved)

            while (iterator1 < removedObjectsPerYear1 or iterator2 < removedObjectsPerYear2 or iterator3 < removedObjectsPerYear3 or iterator4 < removedObjectsPerYear4) and len(listOfCommonRiskPlanets) != 0:
                if iterator1 < removedObjectsPerYear1 and len(listOfCommonRiskPlanets) != 0:

                    highestThreat = max(listOfCommonRiskPlanets, key = itemgetter(1))
                    listOfCommonRiskPlanets.remove(highestThreat)
                    # we do not remove important assets
                    if recentActiveAssets.has_key(highestThreat[0]):
                        continue
                    for d in debris:
                        if highestThreat[0] == d.name:
                            print "Russia removing common", d.name, highestThreat[1]
                            debris.remove(d)
                            listOfRemovedObjects.append(["CIS","c",d.name,0,0,0,0])
                            listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, d.name)
                            iterator1 += 1
                            break

                if iterator2 < removedObjectsPerYear2 and len(listOfCommonRiskPlanets) != 0:

                    highestThreat = max(listOfCommonRiskPlanets, key = itemgetter(1))
                    listOfCommonRiskPlanets.remove(highestThreat)
                    # we do not remove important assets
                    if recentActiveAssets.has_key(highestThreat[0]):
                        continue

                    for d in debris:
                        if highestThreat[0] == d.name:
                            print "US removing common", d.name, highestThreat[1]
                            debris.remove(d)
                            listOfRemovedObjects.append(["US","c",d.name,0,0,0,0])
                            listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, d.name)
                            iterator2 += 1
                            break

                if iterator3 < removedObjectsPerYear3 and len(listOfCommonRiskPlanets) != 0:
                    highestThreat = max(listOfCommonRiskPlanets, key = itemgetter(1))
                    listOfCommonRiskPlanets.remove(highestThreat)
                    # we do not remove important assets
                    if recentActiveAssets.has_key(highestThreat[0]):
                        continue

                    for d in debris:
                        if highestThreat[0] == d.name:
                            print "China removing common", d.name, highestThreat[1]
                            debris.remove(d)
                            listOfRemovedObjects.append(["PRC","c",d.name,0,0,0,0])
                            listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, d.name)
                            iterator3 += 1
                            break

                if iterator4 < removedObjectsPerYear4 and len(listOfCommonRiskPlanets) != 0:
                    highestThreat = max(listOfCommonRiskPlanets, key = itemgetter(1))
                    listOfCommonRiskPlanets.remove(highestThreat)
                    # we do not remove important assets
                    if recentActiveAssets.has_key(highestThreat[0]):
                        continue

                    for d in debris:
                        if highestThreat[0] == d.name:
                            print "EU removing common", d.name, highestThreat[1]
                            debris.remove(d)
                            listOfRemovedObjects.append(["EU","c",d.name,0,0,0,0])
                            listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, d.name)
                            iterator4 += 1
                            break

            for remObj in listOfRemovedObjects:
                for threat1 in pomThreats1:
                    if remObj[2] == threat1[0]:
                        remObj[3] += threat1[1]
                for threat2 in pomThreats2:
                    if remObj[2] == threat2[0]:
                        remObj[4] += threat2[1]
                for threat3 in pomThreats3:
                    if remObj[2] == threat3[0]:
                        remObj[5] += threat3[1]
                for threat4 in pomThreats4:
                    if remObj[2] == threat4[0]:
                        remObj[6] += threat4[1]

            listOfRemovedObjectsFinal.extend(listOfRemovedObjects)
            listOfRemovedObjects = []
            # switching to real run 
            print "Switching to real run ",year

        if virtualRun == 0 and virtualRunEnabled and dataSaving == year:
            dataSaving += 1

            update(recentActiveTLE, ep)
            # increasing number of launched assets over the years
            for i in range(0,(year-2000)):
                launchRandomAssets(recentActiveTLE, debris, recentActiveAssets, indicesIncl, satcatPlanetsLast15years)
            setmass(debris)
            
            agent1RisksProb.append(agent1YearRiskProb)
            agent2RisksProb.append(agent2YearRiskProb)
            agent3RisksProb.append(agent3YearRiskProb)
            agent4RisksProb.append(agent4YearRiskProb)

            commonRisksProb.append(commonRisksProbYear)

            agent1YearRiskProb = 0
            agent2YearRiskProb = 0
            agent3YearRiskProb = 0
            agent4YearRiskProb = 0
            commonRisksProbYear = 0

            noAssetsAgent1 = 0
            noAssetsAgent2 = 0
            noAssetsAgent3 = 0
            noAssetsAgent4 = 0
            # statistics - file output
            for planet in recentActiveAssets.items():
                if "CIS" == planet[1][5].strip():
                    noAssetsAgent1 += 1
                if "US" == planet[1][5].strip():
                    noAssetsAgent2 += 1
                elif "PRC" == planet[1][5].strip():
                    noAssetsAgent3 += 1
                elif planet[1][5].strip() in EU_list:
                    noAssetsAgent4 += 1

            agent1TotalImportantAssets.append(noAssetsAgent1)
            agent2TotalImportantAssets.append(noAssetsAgent2)
            agent3TotalImportantAssets.append(noAssetsAgent3)
            agent4TotalImportantAssets.append(noAssetsAgent4)


        if virtualRun == 0 and year%removeHorizon == 0 and year == yearPom and year != (2016+timeHorizon) and virtualRunEnabled:
             
            virtualRun = 1
            yearPom = year
            virtualDebris = debris[:]

            print "Switching to virtual run ",year

        if not virtualRunEnabled and year == yearPom:
            yearPom += 1

            update(recentActiveTLE, ep)
            # increasing number of launched assets over the years
            for i in range(0,(year-2000)):
                launchRandomAssets(recentActiveTLE, debris, recentActiveAssets, indicesIncl, satcatPlanetsLast15years)
            setmass(debris)


            agent1RisksProb.append(agent1YearRiskProb)
            agent2RisksProb.append(agent2YearRiskProb)
            agent3RisksProb.append(agent3YearRiskProb)
            agent4RisksProb.append(agent4YearRiskProb)

            commonRisksProb.append(commonRisksProbYear)

            agent1YearRiskProb = 0
            agent2YearRiskProb = 0
            agent3YearRiskProb = 0
            agent4YearRiskProb = 0
            commonRisksProbYear = 0

            noAssetsAgent1 = 0
            noAssetsAgent2 = 0
            noAssetsAgent3 = 0
            noAssetsAgent4 = 0
            # statistics - file output
            for planet in recentActiveAssets.items():
                if "CIS" == planet[1][5].strip():
                    noAssetsAgent1 += 1
                if "US" == planet[1][5].strip():
                    noAssetsAgent2 += 1
                elif "PRC" == planet[1][5].strip():
                    noAssetsAgent3 += 1
                elif planet[1][5].strip() in EU_list:
                    noAssetsAgent4 += 1

            agent1TotalImportantAssets.append(noAssetsAgent1)
            agent2TotalImportantAssets.append(noAssetsAgent2)
            agent3TotalImportantAssets.append(noAssetsAgent3)
            agent4TotalImportantAssets.append(noAssetsAgent4)


        if yearToPlot != int(math.ceil((ep-16*365.25)/365.25)) and virtualRun == 0:
            yearToPlot = int(math.ceil((ep-16*365.25)/365.25))
            # statistics - output file
            totalDebrisEvol.append(len(debris))
            # computing space densitites every 100 years
            if yearToPlot == 1 or yearToPlot == 100:
                spaceDensities.append(computeSpatialDensities(debris))
        
        # getting month for relaunching sequences
        month = math.ceil((((ep-16*365.25)/365.25 +0.000001) - math.floor((ep-16*365.25)/365.25))*12)
        monthDouble = (((ep-16*365.25)/365.25 +0.000001) - math.floor((ep-16*365.25)/365.25))*12
        if decadeYearMultiplier != int((year - beginningOfSequence) / 10):
            # reloading recent planets after removing them for monthly relaunch sequence
            recentPlanets,recentActivePlanets = launchSequence(satcat, beginningOfSequence)

        # repeating relaunch sequence every 10 years
        decadeYearMultiplier =  int((year - beginningOfSequence) / 10)
        yearToCopyLaunches = year - 10*decadeYearMultiplier

        # handling the correct format of months
        if month < 10:
            yearAndMonth = str(int(yearToCopyLaunches)) +'-0'+str(int(month))
        else:
            yearAndMonth = str(int(yearToCopyLaunches)) +'-'+str(int(month))

        # repeating previous launching sequence (2005-2015)
        addingNewPlanet = 0;
        if virtualRun == 0:
            for planet in recentPlanets.items():
                if yearAndMonth in planet[1][6]:
                    for planetTLE in recentPlanetsTLE:
                        if planet[0].strip() == planetTLE.name.strip(): 
                            inclin = float(planetTLE.line2.split()[2])
                            # add random number to slightly modify the inclination +/- 2 degrees
                            # modifIncl = float(inclin + (4*random()-2))
                            # inclination is randomly sampled from inclination distribution of objects from (2005-2014)
                            modifIncl = float(np.random.choice(indicesIncl[1:],1,p=weightsIncl))
                            # print new_random_inclination
                            if modifIncl < 0:
                                modifIncl = 0
                            if modifIncl > 180:
                                modifIncl = 180
                            newIncl = '{:8.4f} '.format(modifIncl) # inclination (i)
                            # print newIncl
                            newLine2 = planetTLE.line2[:8] +newIncl+planetTLE.line2[17:]

                            dayOfYear = (((ep-16*365.25)/365.25+2016) - year)*365
                            if dayOfYear < 10:
                                epoch = str(year)[2:4] +"00"+ '{:2.8f}'.format(dayOfYear)
                            elif dayOfYear >= 10 and dayOfYear < 100:
                                epoch = str(year)[2:4] +"0"+ '{:2.8f}'.format(dayOfYear)
                            elif dayOfYear >= 100:
                                epoch = str(year)[2:4] + '{:2.8f}'.format(dayOfYear)

                            newLine1 = planetTLE.line1[:18] + epoch + planetTLE.line1[32:]

                            newObject = kep.planet.tle(newLine1, newLine2)
                            newObject.name = str(year) + newObject.name.strip()[4:]
                            # print newObject.name
                            debris.append(newObject)
                            addingNewPlanet = 1
                            
                            nameOfNew = str(year) + planet[0][4:]
                            newSatcat = planet[1]
                            # new satcat entries does not have modified inclination because the simulation works with TLE objects
                            recentActiveAssets.update({nameOfNew:newSatcat})

                    recentPlanets.pop(planet[0])
            # updating the mass information only if there was a new planet added
            if addingNewPlanet == 1:        
                setmass(debris)

        #getting a list of active recent planets to form list of important assets for each agent
        # a = len(recentActiveAssets) 
        recentActiveAssets = findRecentActiveAssets(recentActiveAssets, year, 10)
        # print "we lost assets", a - len(recentActiveAssets) 


        if virtualRun == 1:
            # print "virtual", len(virtualDebris)
            update(virtualDebris, ep)
            volumes = cube(virtualDebris)
        else:
            # print "real", len(debris)
            decayObjects(debris, year)
            update(debris, ep)
            volumes = cube(debris)

        maxp = 0
        for volume in volumes:
            for p1, p2 in combinations(volume, 2):
                if p1.name == p2.name and p1.mass ==  p2.mass:
                    continue
                # print p1.name,p2.name
                p1.name = p1.name.strip()
                p2.name = p2.name.strip()
                if tuple(p1._v) == tuple(p2._v):
                    pass
                else:
                    # probability of collision
                    Pij = collision_prob(p1, p2) * 2
                    # probability of collision over 5 days
                    P = Pij * 5. * kep.DAY2SEC
                    maxp = max(maxp, P)
                    sump += P
                    if P > 0 and virtualRun == 1:

                        if recentActiveAssets.has_key(p1.name.strip()) and recentActiveAssets.has_key(p2.name.strip()):
                            plnt1 = recentActiveAssets[p1.name.strip()]
                            plnt2 = recentActiveAssets[p2.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), plnt1[5].strip(), plnt2[5].strip(), P])
                        elif recentActiveAssets.has_key(p1.name.strip()) and not recentActiveAssets.has_key(p2.name.strip()):
                            plnt1 = recentActiveAssets[p1.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), plnt1[5].strip(), "ROW", P])
                        elif not recentActiveAssets.has_key(p1.name.strip()) and recentActiveAssets.has_key(p2.name.strip()):
                            plnt2 = recentActiveAssets[p2.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), "ROW", plnt2[5].strip(), P])
                        else:
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(),"ROW", "ROW", P])
                        

                        if recentActiveAssets.has_key(p1.name.strip()):
                            plnt1 = recentActiveAssets[p1.name.strip()]
                            if plnt1[5].strip() == "CIS":
                                for nameOfObject in agent1Threats:
                                    if nameOfObject[0] == p2.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent1Threats.append([p2.name.strip(), P])

                            if plnt1[5].strip() == "US":
                                for nameOfObject in agent2Threats:
                                    if nameOfObject[0] == p2.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent2Threats.append([p2.name.strip(), P])

                            if plnt1[5].strip() == "PRC":
                                for nameOfObject in agent3Threats:
                                    if nameOfObject[0] == p2.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent3Threats.append([p2.name.strip(), P])  

                            if plnt1[5].strip() in EU_list:
                                for nameOfObject in agent4Threats:
                                    if nameOfObject[0] == p2.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent4Threats.append([p2.name.strip(), P])

                        if recentActiveAssets.has_key(p2.name.strip()):
                            plnt2 = recentActiveAssets[p2.name.strip()]
                            if plnt2[5].strip() == "CIS":
                                for nameOfObject in agent1Threats:
                                    if nameOfObject[0] == p1.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent1Threats.append([p1.name.strip(), P])

                            if plnt2[5].strip() == "US":
                                for nameOfObject in agent2Threats:
                                    if nameOfObject[0] == p1.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent2Threats.append([p1.name.strip(), P])

                            if plnt2[5].strip() == "PRC":
                                for nameOfObject in agent3Threats:
                                    if nameOfObject[0] == p1.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent3Threats.append([p1.name.strip(), P])  

                            if plnt2[5].strip() in EU_list:
                                for nameOfObject in agent4Threats:
                                    if nameOfObject[0] == p1.name.strip():
                                        nameOfObject[1] += P
                                        break
                                else:
                                    agent4Threats.append([p1.name.strip(), P])
                       
                    if P > 0 and virtualRun == 0:
                        commonRisksProbYear += P
                        if recentActiveAssets.has_key(p1.name.strip()):
                            planetA = recentActiveAssets[p1.name.strip()]
                            if planetA[5].strip() == "CIS":
                                print "Russian object is threatened by ", p2.name
                                agent1YearRiskProb += P
                            elif planetA[5].strip() == "US":
                                print "American object is threatened by ", p2.name
                                agent2YearRiskProb += P
                            elif planetA[5].strip() == "PRC":
                                print "Chinese object is threatened by ", p2.name
                                agent3YearRiskProb += P
                            elif planetA[5].strip() in EU_list:
                                print "EU object is threatened by ", p2.name
                                agent4YearRiskProb += P

                        if recentActiveAssets.has_key(p2.name.strip()):
                            planetB = recentActiveAssets[p2.name.strip()] 
                            if planetB[5].strip() == "CIS":
                                print "Russian object is threatened by ", p1.name
                                agent1YearRiskProb += P
                            elif planetB[5].strip() == "US":
                                print "American object is threatened by ", p1.name
                                agent2YearRiskProb += P
                            elif planetB[5].strip() == "PRC":
                                print "Chinese object is threatened by ", p1.name
                                agent3YearRiskProb += P
                            elif planetB[5].strip() in EU_list:
                                print "EU object is threatened by ", p1.name
                                agent4YearRiskProb += P  

                    # Monte-carlo simulation - there is a collision if random number is lower than prob of collision
                    if random.random() < P and virtualRun == 0:

                        for planet in recentActiveAssets.items():
                            # if agent's assset is in risk of collision we save the object which threaten it
                            indexInArray = year-2016

                            if (planet[0].strip() == p1.name.strip() or planet[0].strip() == p2.name.strip()) and "CIS" == planet[1][5].strip():
                                agent1Collisions[indexInArray] += 1

                            if (planet[0].strip() == p1.name.strip() or planet[0].strip() == p2.name.strip()) and "US" == planet[1][5].strip():
                                agent2Collisions[indexInArray] += 1

                            if (planet[0].strip() == p1.name.strip() or planet[0].strip() == p2.name.strip()) and "PRC" == planet[1][5].strip():
                                agent3Collisions[indexInArray] += 1

                            if (planet[0].strip() == p1.name.strip() or planet[0].strip() == p2.name.strip()) and planet[1][5].strip() in EU_list:
                                agent4Collisions[indexInArray] += 1

                        yearOfColl.append(year)

                        print
                        print "!" * 25 + " BOOOOOM " + "!" * 25
                        print
                        print 'planet', p1.name, 'with mass', p1.mass
                        print 'planet', p2.name, 'with mass', p2.mass
                        dv = np.linalg.norm(p1._v - p2._v)
                        print 'collision velocity', dv
                        if p2.mass < p1.mass:
                            catastrophRatio = (p2.mass*dv**2)/(2*p1.mass*1000)
                        else:
                            catastrophRatio = (p1.mass*dv**2)/(2*p2.mass*1000)
                        print 'catastrophic ratio:', catastrophRatio
                        if catastrophRatio<40:
                            print 'NON CATASTROPHIC COLLISION - NO DEBRIS'
                            numberOfNonCatastrophicCollisions += 1
                        else:
                            debris1, debris2 = breakup(ep, p1, p2)
                            print (len(debris1) + len(debris2)), 'new pieces from collision'
                            totalNumberOfNewDebris += (len(debris1) + len(debris2))
                            debris.extend(debris1)
                            debris.extend(debris2)
                            numberOfCatastrophicCollisions += 1
                            setmass(debris)
                            # plot_orbits(ep, debris1, debris2)
        # if year > 2020:
        #     end = time.time()
        #     print "Collision computation takes:", (end-start)

        print '%.2f %d %d %d %10.8f' % (float(ep)/365.25 - 16, len(debris), len(volumes), max(map(len, volumes)) if len(volumes) else 0, maxp)
        # if year > 2020:
        #     end1 = time.time()
        #     print "One run takes:", (end1-start1)
    # print objectsWithPossitiveProbOfCollision
    # plot_orbits(ep, recentActiveAssets, recentActiveAssets)
    # print len(recentActiveAssets)
    # plot_orbits(ep, recentActiveTLE[-200:], debris[-200:])
        
    print 'There were', numberOfCatastrophicCollisions,'catastrophic collisions'
    print 'There were', numberOfNonCatastrophicCollisions,'NON-catastrophic collisions'
    print 'From collisions there were', totalNumberOfNewDebris,'of new debris'
    years = range(0,yearToPlot+1,1)

    totalDebrisEvol.append(len(debris))
    spaceDensities.append(computeSpatialDensities(debris))


    fileName = "5cmObjects_"+stringOfRemovals+"_2xProb" + str(qq) + ".txt"
    f = open(fileName,'w')
    f.write("%s\n" % years)
    f.write("%s\n" % totalDebrisEvol)
    f.write("%s\n" % numberOfCatastrophicCollisions)
    f.write("%s\n" % numberOfNonCatastrophicCollisions)
    f.write("%s\n" % totalNumberOfNewDebris)
    for dens in spaceDensities:
        f.write("%s\n" % dens)
    f.close()


    fileName = "agentsRisks"+stringOfRemovals+"obj" + str(qq) + ".txt"
    f = open(fileName,'w')
    f.write("%s\n" % agent1Collisions)
    f.write("%s\n" % agent2Collisions)
    f.write("%s\n" % agent3Collisions)
    f.write("%s\n" % agent4Collisions)
    # f.write("%s\n" % agent1Risks)
    # f.write("%s\n" % agent2Risks)
    # f.write("%s\n" % agent3Risks)
    # f.write("%s\n" % agent4Risks)
    f.write("%s\n" % agent1RisksProb)
    f.write("%s\n" % agent2RisksProb)
    f.write("%s\n" % agent3RisksProb)
    f.write("%s\n" % agent4RisksProb)

    f.write("%s\n" % commonRisksProb)
    f.write("%s\n" % agent1TotalImportantAssets)
    f.write("%s\n" % agent2TotalImportantAssets)
    f.write("%s\n" % agent3TotalImportantAssets)
    f.write("%s\n" % agent4TotalImportantAssets)
    f.close()

    fileName = "listOfRemoved"+stringOfRemovals+"obj" + str(qq) + ".txt"
    f = open(fileName,'w')
    for remObject in listOfRemovedObjectsFinal:
        f.write("%s\n" % remObject)
    f.close()
