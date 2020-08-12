from __future__ import print_function
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import re
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import copy




def IsG01(line):
  l = line.split()
  Xvalues = re.findall(r'X[-+]?\d*\.?\d+',line)
  Yvalues = re.findall(r'Y[-+]?\d*\.?\d+',line)
  if(len(Xvalues) == 0 and len(Yvalues) == 0):
    return False
  if "G1" in l:
    return True
  if "G0" in l:
    return True
  return False

def GetDistance(x1,y1,x2,y2):
  return math.sqrt(math.pow(float(x1)-float(x2),2) +math.pow(float(y1)-float(y2),2))


def GetEValue(Line):
  Evalues = re.findall(r'E[-+]?\d*\.?\d+',Line)
  if(len(Evalues) > 0):
    return float(Evalues[0][1:])
  else:
    return -1

def GetXValue(Line):
  Xvalues = re.findall(r'X[-+]?\d*\.?\d+',Line)
  if(len(Xvalues) > 0):
    GetXValue.PreviousSpeed = Xvalues[0][1:]
    return float(Xvalues[0][1:])
  else:
    return float(GetXValue.PreviousSpeed)
GetXValue.PreviousSpeed = 0

def GetZValue(Line):
  Zvalues = re.findall(r'Z[-+]?\d*\.?\d+',Line)
  if(len(Zvalues) > 0):
    GetZValue.PreviousSpeed = Zvalues[0][1:]
    return float(Zvalues[0][1:])
  else:
    return float(GetZValue.PreviousSpeed)
GetZValue.PreviousSpeed = 0

def GetFValue(Line):
  Fvalues = re.findall(r'F[-+]?\d*\.?\d+',Line)
  if(len(Fvalues) > 0):
    GetFValue.PreviousSpeed = Fvalues[0][1:]
    return float(Fvalues[0][1:])
  else:
    return float(GetFValue.PreviousSpeed)
GetFValue.PreviousSpeed = 0

def GetYValue(Line):
  Yvalues = re.findall(r'Y[-+]?\d*\.?\d+',Line)
  if(len(Yvalues) > 0):
    GetYValue.PreviousSpeed = Yvalues[0][1:]
    return float(Yvalues[0][1:])
  else:
    return float(GetYValue.PreviousSpeed)
GetYValue.PreviousSpeed = 0
"""
with open(filename,'r') as file:
  data = file.read()

data = data.split('\n')

Header = [['GCode','Xposition','Yposition','FValue','Evalue','Zvalue']]
DataPoints = [copy.deepcopy(Header)]
CurrentSet = 0
LastZ = 0

Optimal_Paths = []

def CalculateTour(CurrentSet):
  NodeLookup = [["NameNum","X","Y","Gcode"]]
  NodeWeights = []
  for Point in range(2,len(DataSet)):
    if(Point%100 == 0):
      print("Currently calculating point {}".format(Point))
    if(DataSet[Point][0] == "G1"):
      CurrentList = []
      for k in NodeLookup[1:]:
        CurrentList.append(int(GetDistance(k[1],k[2],DataSet[Point-1][1],DataSet[Point-1][1]) + 2000))
      NodeLookup.append(['a{}'.format(Point),DataSet[Point-1][1],DataSet[Point-1][2]])
      CurrentList.append(0)
      NodeWeights.append(CurrentList)

      #Now time to add the connected line segment
      CurrentList = []
      for k in NodeLookup[1:-1]:
        CurrentList.append(int(GetDistance(k[1],k[2],DataSet[Point][1],DataSet[Point][1]) + 2000))
      NodeLookup.append(['b{}'.format(Point),DataSet[Point][1],DataSet[Point][2]])
      CurrentList.append(0)
      CurrentList.append(0)
      NodeWeights.append(CurrentList)

  print("{} Nodes in NodeLookup".format(len(NodeLookup)))
  print("{} Nodes in NodeWeights".format(len(NodeWeights)))

  for i in NodeWeights[:-1]:
    if(NodeWeights.index(i)%100 == 0):
      print("Fixing array {}".format(NodeWeights.index(i)))
    while len(i) != len(NodeWeights[-1]):
      i.append(0)

  #print(NodeWeights)

  Optimization_array = np.array(NodeWeights)
  if(len(NodeWeights) < 2):
    continue


  #TIME TO SOLVE!!!!!!
  from concorde.tsp import TSPSolver

  print("Printing to TSPLIB file format")
  solver = TSPSolver.from_data_explicit(weightMatrix = Optimization_array,norm="GEO")

  print("Starting Concorde Solver")
  Optimal_Path = solver.solve(verbose=False,time_bound=10)

  SolutionGcode = ""

  currentpos = [0,0]
  for node in Optimal_Path:
    if(node%2 == 1):
      if(NodeLookup[node+1][1] != currentpos[0] or NodeLookup[node+1][2] != currentpos[1]):
        currentpos = [NodeLookup[node+1][1],NodeLookup[node+1][2]]
        SolutionGcode += "G0{} X{} Y{}".format(" 7200",)



OptimizedGcode = []
Startx,Starty = 0

for line in data:
  if len(re.findall(r"G1",line)) > 0:
    DataPoints[CurrentSet].append(['G1',GetXValue(line),GetYValue(line),GetFValue(line),GetEValue(line),GetZValue(line)])
  elif (len(re.findall(r"G0",line)) > 0):
    if (GetZValue(line) != -1):
      OptimizedGcode += CalculateTour(CurrentSet)
      DataPoints.append(copy.deepcopy(Header))
      CurrentSet += 1
      Start = []
    DataPoints[CurrentSet].append(['G0',GetXValue(line),GetYValue(line),GetFValue(line),GetEValue(line),GetZValue(line)])
  elif (len(re.findall(r';',line)) > 0 ):
    OptimizedGcode += line + '\n'
  else:
    OptimizedGcode += CalculateTour(CurrentSet)
    OptimizedGcode += line + '\n'
    CurrentSet += 1
    DataPoints.append(copy.deepcopy(Header))


with open(filename.replace(".","optimized."),"w") as file:
  file.write(OptimizedGcode)
"""
def create_data_model(dataTSP,start):
    """Stores the data for the problem."""
    data = {}
    data['distance_matrix'] = dataTSP  # yapf: disable
    data['num_vehicles'] = 1
    data['depot'] = [start]
    
    data['ends'] = [len(dataTSP)-1]
    return data


def print_solution(manager, routing, assignment):
    """Prints assignment on console."""
    index = routing.Start(0)
    plan_output = 'Route for vehicle 0:\n'
    route_distance = 0
    Fulllist = []
    while not routing.IsEnd(index):
        Fulllist.append(manager.IndexToNode(index))
        plan_output += ' {} ->'.format(manager.IndexToNode(index))
        previous_index = index
        index = assignment.Value(routing.NextVar(index))
        route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)
    plan_output += ' {}\n'.format(manager.IndexToNode(index))
    print('Route distance: {}miles\n'.format(route_distance))
    return Fulllist


def CalcTSP(dataTSP,StartDepo):
    """Entry point of the program."""
    # Instantiate the data problem.
    data = create_data_model(dataTSP,StartDepo)

    # Create the routing index manager.
    #print(data['depot'])
    #print(data['ends'])
    manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),
                                           data['num_vehicles'], data['depot'],data['ends'])

    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)


    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['distance_matrix'][from_node][to_node]

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
    search_parameters.time_limit.seconds = 100000

    # Solve the problem.
    assignment = routing.SolveWithParameters(search_parameters)

    #print(assignment)
    #input()

    # Print solution on console.
    if assignment:
        return print_solution(manager, routing, assignment)

BforeTot = 0
AfterTot = 0
def toGCodeFormat(line):
  if(line["G"] == 1):
    return "G1 X{} Y{} E{} F{}".format(line['X'][1],line['Y'][1],line['E'],line['F']) + ("" if line['Z'][0] == line['Z'][1] else " Z{}".format(line['Z'][1]))
  else:
    return "G0 X{} Y{} F{}".format(line['X'][1],line['Y'][1],line['F']) + ("" if line['Z'][0] == line['Z'][1] else " Z{}".format(line['Z'][1]))

def toGCodeEFormat(line,Eval):
  if(line["G"] == 1):
    return [Eval +line['E'],  "G1 X{} Y{} E{} F{}".format(line['X'][1],line['Y'][1],Eval + line['E'],line['F']) + ("" if line['Z'][0] == line['Z'][1] else " Z{}".format(line['Z'][1]))]
  else:
    return [Eval,"G0 X{} Y{} F{}".format(line['X'][1],line['Y'][1],line['F']) + ("" if line['Z'][0] == line['Z'][1] else " Z{}".format(line['Z'][1]))]
  

def toLineFormat(Line,StartPos,Speed):
  G = 1  if Line[1] == "1" else 0
  X = GetXValue(Line)
  Y = GetYValue(Line)
  E = 0 if GetEValue(Line) == -1 else (GetEValue(Line) - StartPos[3])
  F = Speed if GetFValue(Line) == -1 else GetFValue(Line)

  Z = StartPos[2] if GetZValue(Line) == -1 else GetZValue(Line)
  return {'G' : G, 'X' : [StartPos[0],X], 'Y' : [StartPos[1],Y], 'E' : E  , 'F' : F, 'Z' : [StartPos[2],Z]}

def flipLine(line):
  return {'G' : line['G'], 'X' : [line['X'][1],line['X'][0]], 'Y' : [line['Y'][1],line['Y'][0]], 'E' : line['E'] , 'F' : line['F'], 'Z' : line['Z']}

def solve(ListOfGCommands,StartingPos):
  global BforeTot
  global AfterTot
  Lines = []
  Pos = StartingPos
  total = 0
  speed = 10

  for command in ListOfGCommands:
    currentLine = toLineFormat(command,Pos,speed)
    if(currentLine["G"] == 0):
      total = total + GetDistance(currentLine["X"][0],currentLine["Y"][0],currentLine["X"][1],currentLine["Y"][1])
    Lines.append(currentLine)
    speed = currentLine["F"]
    Pos = [currentLine["X"][1],currentLine["Y"][1],currentLine["Z"][1],Pos[3] + currentLine["E"]]
  
  G0Speed = speed
  for i in Lines:
    if(i["G"] == 0):
      G0Speed = i['F']
      break

  ####################################################
  #                SOLVE THE F******                 #
  ####################################################

  ##Convert Line format into TSPLIB
  TSPDATA = []
  newLines = []
  startingPoint = 0
  maxDist = -1
  for i in Lines:
    if(i["G"] == 1):
      newLines.append(i)

  
  if(newLines == 0):
    return {'Gcode' : [], 'EndPoint' : [StartingPos[0],StartingPos[1],Pos[2],StartingPos[3]],"Good" : True}
  if(len(newLines) < 4):
    solution = {'Gcode' : [], 'EndPoint' : Pos,"Good" : False}
    for Line in Lines:
      solution['Gcode'].append(toGCodeFormat(Line))
    return(solution)

  ConstantFact = 10000000
  for i in newLines:
    DistanceArray = []
    DistanceArray2 = []
    for j in newLines:
      if(i == j):
        DistanceArray.append(0)
        DistanceArray.append(0)
        DistanceArray2.append(0)
        DistanceArray2.append(0)
      else:
        DistanceArray.append(GetDistance(i["X"][0],i["Y"][0],j["X"][0],j["Y"][0]) + ConstantFact)
        DistanceArray.append(GetDistance(i["X"][0],i["Y"][0],j["X"][1],j["Y"][1])+ ConstantFact)
        DistanceArray2.append(GetDistance(i["X"][1],i["Y"][1],j["X"][0],j["Y"][0])+ ConstantFact)
        DistanceArray2.append(GetDistance(i["X"][1],i["Y"][1],j["X"][1],j["Y"][1])+ ConstantFact)

    ##Find distance from starting point
    Dist1 = GetDistance(i["X"][0],i["Y"][0],StartingPos[0],StartingPos[1])
    Dist2 = GetDistance(i["X"][1],i["Y"][1],StartingPos[0],StartingPos[1])
    #print("Point1 = {}X {}Y".format(i["X"][0],i["Y"][0]))
    if(maxDist >= Dist1 or maxDist == -1):
      maxDist = Dist1
      startingPoint = newLines.index(i)*2
      #print("Found new min: Point {}".format(startingPoint))
    elif(maxDist >= Dist2):
      maxDist > Dist2
      startingPoint = newLines.index(i)*2 + 1
      #print("Found new min: Point {}".format(startingPoint))

    DistanceArray.append(0)
    DistanceArray2.append(0)
    TSPDATA.append(DistanceArray)
    TSPDATA.append(DistanceArray2)
  TSPDATA.append([0 for x in range(len(TSPDATA[0]))])
  #print(startingPoint)
  TSP = CalcTSP(TSPDATA,startingPoint)
  print("Starting distance: {}".format(total))

  solved = []
  if(TSP[0]%2 == 0):
    solved.append({'G' : 0, 'X' : [StartingPos[0],newLines[int(TSP[0]/2)]['X'][0]], 'Y' : [StartingPos[0],newLines[int(TSP[0]/2)]['Y'][0]], 'E' : -1 , 'F' : G0Speed, 'Z' : [StartingPos[2],Pos[2]]})
  else:
    solved.append({'G' : 0, 'X' : [StartingPos[0],newLines[int(TSP[0]/2)]['X'][1]], 'Y' : [StartingPos[0],newLines[int(TSP[0]/2)]['Y'][1]], 'E' : -1 , 'F' : G0Speed, 'Z' : [StartingPos[2],Pos[2]]})
 
  currentE = StartingPos[3]
  print(StartingPos)
  for i in range(0,len(TSP),2):
    #print(TSP[i]/2)
    #print(newLines[int(TSP[i]/2)])
    #print("Dist1 = {}".format(GetDistance(solved[-1]['X'][1], solved[-1]['Y'][1], newLines[int(TSP[i]/2)]['X'][0],newLines[int(TSP[i]/2)]['Y'][0])))
    #print("Dist2 = {}".format(GetDistance(solved[-1]['X'][1], solved[-1]['Y'][1], newLines[int(TSP[i]/2)]['X'][1],newLines[int(TSP[i]/2)]['Y'][1])))
    if(TSP[i]%2 == 0):
      if GetDistance(solved[-1]['X'][1], solved[-1]['Y'][1], newLines[int(TSP[i]/2)]['X'][0],newLines[int(TSP[i]/2)]['Y'][0]) > .1:
        solved.append({'G' : 0, 'X' : [solved[-1]['X'][1],newLines[int(TSP[i]/2)]['X'][0]], 'Y' : [solved[-1]['Y'][1],newLines[int(TSP[i]/2)]['Y'][0]], 'E' : -1 , 'F' : G0Speed, 'Z' : [solved[-1]['Z'][0],solved[-1]['Z'][1]]})
      solved.append(newLines[int(TSP[i]/2)])
    else:

      if GetDistance(solved[-1]['X'][1], solved[-1]['Y'][1], newLines[int(TSP[i]/2)]['X'][1],newLines[int(TSP[i]/2)]['Y'][1]) > .1:
        solved.append({'G' : 0, 'X' : [solved[-1]['X'][1],newLines[int(TSP[i]/2)]['X'][1]], 'Y' : [solved[-1]['Y'][1],newLines[int(TSP[i]/2)]['Y'][1]], 'E' : -1 , 'F' : G0Speed, 'Z' : [solved[-1]['Z'][0],solved[-1]['Z'][1]]})
      solved.append(flipLine(newLines[int(TSP[i]/2)]))

  #[print(x) for x in solved]
  BforeTot = BforeTot + total

  total = 0
  solution = {'Gcode' : [], 'EndPoint' : []}
  Line = []
  for Line in solved:
    line = toGCodeEFormat(Line,currentE)
    currentE = line[0]
    #print(currentE)
    if(Line["G"] == 0 and solved.index(Line) != 0):
      total = total + GetDistance(Line["X"][0],Line["Y"][0],Line["X"][1],Line["Y"][1])
    solution['Gcode'].append(line[1])

  solution = {'Gcode' : solution['Gcode'], 'EndPoint' : [solved[-1]['X'][1],solved[-1]['Y'][1],solved[-1]['Z'][1],Pos[3]],"Good" : True}
  print("Ending distance: {}".format(total))
  AfterTot = AfterTot + total
  return(solution)

def main(filePath):
  #Load gcode delimit by line
  FinalGcode = []
  data = ""
  with open(filePath,'r') as file:
    data = file.read()
    data = data.split('\n')
  #Next we split up the lines based on what they start with. If they are G# then group them, otherwise just put them on the stack

  #CurrentPos is of the format X,Y,Z,E
  currentPos = [0,0,0,0]
  group = []
  previousZ = 0
  for line in data:
    if(IsG01(line) and GetZValue(line) == previousZ):
      group.append(line)
      #Start grouping
    else:
      #If its a new layer
      if(len(group) != 0):
          solution = solve(group,currentPos)
          if solution["Good"] :
            FinalGcode = FinalGcode + solution["Gcode"]
          else: 
            FinalGcode = FinalGcode + group
          currentPos = solution["EndPoint"] 
          #print(currentPos)
          if(currentPos[3] < 0):
            currentPos[3] = 0
          group = []
          if(GetZValue(line) != previousZ):
            group.append(line)
          else:
            FinalGcode.append(line)
      else:
        if(GetZValue(line) != previousZ):
            group.append(line)
        else:
          FinalGcode.append(line)
    previousZ = currentPos[2]
  return '\n'.join(FinalGcode)
    





filename = sys.argv[1]
OptimizedGcode = main(filename)

with open(filename.replace(".","optimized."),"w") as file:
  file.write(OptimizedGcode)

print(BforeTot)
print(AfterTot)