#Author: Jesse Evans
#Version 3.3
#Derived from the 3.2 version of FDMTSPOptimizer this version tries to make larger subproblems in order to understand subproblem effects on the solving process.
#8/11/2020

#Imports
from __future__ import print_function
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
#from tsp_solver.greedy import solve_tsp
import re
import sys
import math
import numpy as np
import copy
import datetime

Version = 3.3

#####################################################################
                        #Section 1: GCode Parsers and Simple Functions#
#####################################################################


#Returns the Length of a distance matrix
def GetLength(Order,Matrix):
  l_length = 0
  lastNum = -1
  for i in Order:
    if lastNum == -1:
      lastNum = i
      continue
    l_length = l_length + Matrix[lastNum][i]
    lastNum = i
  return l_length

#Given a string of text, returns true if the line starts with either G0 or G1
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


#Given a list of bucket sizes, and a index for all of them, returns which bucket its in
def BucketSearch(Index, Buckets):
  count = 0
  for x in Buckets:
    Index = Index - x
    if Index < 0:
      return count
    count = count + 1

#; are important
def GetSemiColon(line):
  Colon = re.findall(r';',line)
  if(len(Colon) > 0):
    return True
  return False

#A simple Distance formula that returns a distance between point (x1,y1) and point (x2,y2)
def GetDistance(x1,y1,x2,y2):
  return math.sqrt(math.pow(float(x1)-float(x2),2) +math.pow(float(y1)-float(y2),2))

#Returns the GValue specified by a G0 or G1 line of gcode
#Test simply is here for continuity
def GetGValue(Line, Test = False):
  Gvalues = re.findall(r'G[-+]?\d*\.?\d+',Line)
  if(len(Gvalues) > 0):
    return int(Gvalues[0][1:])
  else:
    return -1


#Returns the EValue specified by a G0 or G1 line of gcode
#Test simply is here for continuity
def GetEValue(Line, Test = False):
  Evalues = re.findall(r'E[-+]?\d*\.?\d+',Line)
  if(len(Evalues) > 0):
    return float(Evalues[0][1:])
  else:
    return -1

#Returns a X value specifed by a G0 or G1 line of gcode
#If no X value is specified then it will return the last X value. Test overrides
def GetXValue(Line,Test = False):
  Xvalues = re.findall(r'X[-+]?\d*\.?\d+',Line)
  if(len(Xvalues) > 0):
    GetXValue.PreviousSpeed = Xvalues[0][1:]
    return float(Xvalues[0][1:])
  else:
    return float(GetXValue.PreviousSpeed) if not Test else -1
GetXValue.PreviousSpeed = 0


#Returns a Z value specifed by a G0 or G1 line of gcode
#If no X value is specified then it will return the last X value. Test overrides
def GetZValue(Line,Test = False):
  Zvalues = re.findall(r'Z[-+]?\d*\.?\d+',Line)
  if(len(Zvalues) > 0):
    GetZValue.PreviousSpeed = Zvalues[0][1:]
    return float(Zvalues[0][1:])
  else:
    return float(GetZValue.PreviousSpeed) if not Test else -1
GetZValue.PreviousSpeed = 0


#Returns a F value specifed by a G0 or G1 line of gcode
#If no F value is specified then it will return the last F value. Test overrides
def GetFValue(Line,Test = False):
  Fvalues = re.findall(r'F[-+]?\d*\.?\d+',Line)
  if(len(Fvalues) > 0):
    GetFValue.PreviousSpeed = Fvalues[0][1:]
    return float(Fvalues[0][1:])
  else:
    return float(GetFValue.PreviousSpeed) if not Test else -1
GetFValue.PreviousSpeed = 0


#Returns a Y value specifed by a G0 or G1 line of gcode
#If no Y value is specified then it will return the last Y value. Test overrides
def GetYValue(Line, Test = False):
  Yvalues = re.findall(r'Y[-+]?\d*\.?\d+',Line)
  if(len(Yvalues) > 0):
    GetYValue.PreviousSpeed = Yvalues[0][1:]
    return float(Yvalues[0][1:])
  else:
    return float(GetYValue.PreviousSpeed) if not Test else -1
GetYValue.PreviousSpeed = 0

#####################################################################
              #Section 2: Google Ortools Solver#
#####################################################################

def create_data_model(dataTSP,start,end):
    """Stores the data for the problem."""
    data = {}
    data['distance_matrix'] = dataTSP  # yapf: disable
    data['num_vehicles'] = 1
    data['depot'] = [start]
    
    data['ends'] = [end]
    return data


def get_solution(manager, routing, assignment):
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
    #print('Route distance: {}miles\n'.format(route_distance))
    #print(plan_output)
    return Fulllist


def CalcTSP(dataTSP,StartDepo,EndDepo):
    """Entry point of the program."""
    # Instantiate the data problem.
    data = create_data_model(dataTSP,StartDepo,EndDepo)

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
    search_parameters.time_limit.seconds = 100

    # Solve the problem.
    assignment = routing.SolveWithParameters(search_parameters)

    #print(assignment)
    #input()

    # Return assignment
    if assignment:
        return get_solution(manager, routing, assignment)

#####################################################################
                    #Section 3: Class Definitions#
#####################################################################








#Object for printed lines
class Line:
  GValue = -1
  XValue = [-1,-1]
  EValue = -1
  FValue = [-1,-1]
  ZValue = [-1,-1]
  Z0Speed = -1

  #Fills the data in the class based on given GCode
  #self - self
  #Gcode -  A string containing a line of valid GCode
  #StartPos - A starting Position to know where the line started from
  #Speed - The default speed that this print should have in case no F value is given
  def __init__(self, Gcode, StartPos, Speed):
    G = 1  if Gcode[1] == "1" else 0
    X = GetXValue(Gcode)
    Y = GetYValue(Gcode)
    E = 0 if GetEValue(Gcode) == -1 else (GetEValue(Gcode) - StartPos[3])
    F = Speed if GetFValue(Gcode) == -1 else GetFValue(Gcode)
    Z = StartPos[2] if GetZValue(Gcode) == -1 else GetZValue(Gcode)
    
    try:
      if(E/GetDistance(StartPos[0],StartPos[1],X,Y) > 10 ):
        print("WARNING")
        print("Very large E Value found")
        print("It is very likely that this is a bug")
        print("E Value = " , E)
        print("Distance was ",GetDistance(StartPos[0],StartPos[1],X,Y))
        print(Gcode)
        print("Press any key to continue")
        input()
    except(ZeroDivisionError):
      #Probably just a weird error??
      pass
    self.GValue = G
    self.XValue = [StartPos[0],X]
    self.YValue = [StartPos[1],Y]
    self.EValue = E
    self.FValue = F
    self.ZValue = [StartPos[2],Z]
  
  #Generates Gcode from this line object, this function only does the end point of the line,
  #and assumes that the printer is already at the starting point
  def ToGCode(self,Eval):
    if(self.GValue == 1):
      return "G1 X{} Y{} E{} F{}".format(self.XValue[1],self.YValue[1],Eval + self.EValue,self.FValue) + ("" if self.ZValue[0] == self.ZValue[1] else " Z{}".format(self.ZValue[1]))
    else:
      return "G0 X{} Y{} F{}".format(self.XValue[1],self.YValue[1],self.FValue) + ("" if self.ZValue[0] == self.ZValue[1] else " Z{}".format(self.ZValue[1]))

  #Debug print
  def Dump(self):
    print("LINE DUMP")
    print(self.ToGCode())

  def Swap(self):
    self.XValue = [self.XValue[1],self.XValue[0]]
    self.YValue = [self.YValue[1],self.YValue[0]]










#Sub problems are what the problem is broken into first. These subproblems are created by using code other then valid G0 or G1 being present.
#For the purposes 
class SubProblem:

  #A list of line objects 
  Lines = []
  #This is code used as a delimiter but is still very important to the overall printing process
  #It is stored here to be stiched back in later
  Preamble = ""

  #These are the starting and ending points. This is used for various calculations and the final solved product
  EndPoint = [-1,-1,-1]
  StartPoint = [-1,-1,-1]

  Result = ""

  def __init__(self,l_Lines,l_Preamble, l_Start = [-1,-1,-1],l_End = [-1,-1,-1]):
    self.Lines = l_Lines
    self.Preamble = l_Preamble
    self.StartPoint = l_Start
    self.EndPoint = l_End
  
  #Debug print
  def Dump(self):
    print()
    print("SUBPROBLEM DUMP")
    print("{} Lines".format(len(self.Lines)))
    print()
    print(self.Preamble)
    for line in self.Lines:
      pass
      #line.Dump()
  
  def GetLines(self):
    return [x for x in self.Lines]

  def GetSize(self):
    return len(self.Lines)

  #Solves the subproblem and reorders the lines
  def Solve(self, Start, End, Open = False):

    #TODO::BUG
    #When the Start and End are on the same line, and there is more then one line, we need to find the point closest to the end point, and set that as the new End.
    if(len(self.Lines) > 1):
      if(Start%2 == 0):
        if Start + 1 == End:
          Dist = 99999
          lead = 0
          for line in self.Lines:
            if line != self.Lines[int((End-1)/2)]:
              d1 = (GetDistance(self.Lines[int((End-1)/2)].XValue[1],self.Lines[int((End-1)/2)].YValue[1],line.XValue[0],line.YValue[0]))
              if d1 < Dist:
                lead = self.Lines.index(line) *2 
                Dist = d1
              d2 = (GetDistance(self.Lines[int((End-1)/2)].XValue[1],self.Lines[int((End-1)/2)].YValue[1],line.XValue[1],line.YValue[1]))
              if d2 < Dist:
                lead = self.Lines.index(line) *2 + 1 
                Dist = d2
          End = lead

      else:
        if Start - 1 == End:
          Dist = 99999
          lead = 0
          for line in self.Lines:
            if line != self.Lines[int((End)/2)]:
              d1 = (GetDistance(self.Lines[int((End)/2)].XValue[1],self.Lines[int((End)/2)].YValue[1],line.XValue[0],line.YValue[0]))
              if d1 < Dist:
                lead = self.Lines.index(line) *2 
                Dist = d1
              d2 = (GetDistance(self.Lines[int((End)/2)].XValue[1],self.Lines[int((End)/2)].YValue[1],line.XValue[1],line.YValue[1]))
              if d2 < Dist:
                lead = self.Lines.index(line) *2 + 1 
                Dist = d2
          End = lead 

    DistanceMatrix = []
    for line in range(len(self.Lines)):
      lineIndex = line
      Eval1 = []
      Eval2 = []

      for Node in DistanceMatrix:
        NodeIndex = DistanceMatrix.index(Node)
        #print("Node: {}".format(NodeIndex))
        #print("line: {}".format(lineIndex))
        if NodeIndex %2 == 0:
          d1 = GetDistance(self.Lines[int(NodeIndex/2)].XValue[0],self.Lines[int(NodeIndex/2)].YValue[0],self.Lines[int(lineIndex)].XValue[0],self.Lines[int(lineIndex)].YValue[0])
          d2 = GetDistance(self.Lines[int(NodeIndex/2)].XValue[0],self.Lines[int(NodeIndex/2)].YValue[0],self.Lines[int(lineIndex)].XValue[1],self.Lines[int(lineIndex)].YValue[1])
          Node.append(d1)
          Node.append(d2)
          Eval1.append(d1)
          Eval2.append(d2)
        else:
          d1 = GetDistance(self.Lines[int((NodeIndex - 1)/2)].XValue[1],self.Lines[int((NodeIndex - 1)/2)].YValue[1],self.Lines[int(lineIndex)].XValue[0],self.Lines[int(lineIndex)].YValue[0])
          d2 = GetDistance(self.Lines[int((NodeIndex - 1)/2)].XValue[1],self.Lines[int((NodeIndex - 1)/2)].YValue[1],self.Lines[int(lineIndex)].XValue[1],self.Lines[int(lineIndex)].YValue[1])
          Node.append(d1)
          Node.append(d2)
          Eval1.append(d1)
          Eval2.append(d2)
      Eval1.append(-1)
      Eval1.append(-1)
      Eval2.append(-1)
      Eval2.append(-1)
      DistanceMatrix.append(Eval1)
      DistanceMatrix.append(Eval2)
    
    if Open :
      DistanceMatrix.append([0 for x in DistanceMatrix])
      for x in DistanceMatrix:
        x.append(0)

    Total = sum([sum(x) + 2 for x in DistanceMatrix])

    DistanceMatrix = [[((x + Total) if (x != -1) else 0) for x in v] for v in DistanceMatrix]

    #print()
    Assignment = CalcTSP(DistanceMatrix,Start,(len(DistanceMatrix) - 1) if Open else End)

    
    print("Previous Rout Length was {}".format(GetLength([x for x in range(len(Assignment))],DistanceMatrix)))
    print("New Rout Length was {}".format(GetLength(Assignment,DistanceMatrix)))
    NewOrder = []
    for x in Assignment[::2]:
      if(x %2 != 0 ):
        self.Lines[int((x - 1) /2)].Swap()
        NewOrder.append(self.Lines[int((x-1)/2)])
      else:
        NewOrder.append(self.Lines[int(x/2)])
    self.Lines = NewOrder

    self.StartPoint = [self.Lines[0].XValue[0], self.Lines[0].YValue[0], -1]
    self.EndPoint = [self.Lines[-1].XValue[1], self.Lines[-1].YValue[1], -1]

    return 
  
  #Returns the Endpoint
  #If endpoint has not been set then it will return the last line end in this problem
  #If no lines are present, it will attempt to find a valid end point in the preamble
  #If that also fails it will return the StartPoint in a last ditch effort
  #WARNING THIS COULD RETURN [-1,-1,-1] IF NO START HAS BEEN SET EITHER
  def Travel(self):
    if(self.EndPoint != [-1,-1,-1]):
      return self.EndPoint
    elif(len(self.Lines) != 0):
      return [self.Lines[-1].XValue[1], self.Lines[-1].XValue[1], self.EndPoint[2]]

      ##TODO:: Preamble stuff I'm not sure what to do with it yet
    else:
      return self.StartPoint


  def Clean(self,FirstLayer = False):
    data = self.Preamble.split('\n')
    NewPreamble = ""
    for line in data:
      if(IsG01(line)):
        if(GetXValue(line,True) or GetYValue(line,True)):
          if not FirstLayer and not GetSemiColon(line):
            #print(line)
            continue
      NewPreamble = NewPreamble + line + '\n'
    self.Preamble = NewPreamble 

  
  def Stitch(self,Eval,CurrentPos,Height):
    
    Result = self.Preamble

    if(len(self.Lines) < 1):
      self.Result = Result
      return [Eval,CurrentPos]

    if(CurrentPos[0] != self.Lines[0].XValue[0]) and (CurrentPos[1] != self.Lines[0].YValue[0]):
      pass
    Result = Result + "G0 X{} Y{} F{} Z{}".format(self.Lines[0].XValue[0],self.Lines[0].YValue[0],self.Lines[0].FValue*4.5,Height) + '\n'
    Result = Result + self.Lines[0].ToGCode(Eval) + '\n'
    Eval = Eval + self.Lines[0].EValue
    for line in self.Lines[1:]:
      Result = Result + "G0 X{} Y{} F{} Z{}".format(line.XValue[0],line.YValue[0],line.FValue*4.5,Height) + '\n'
      Result = Result + line.ToGCode(Eval) +  '\n'
      Eval = Eval + line.EValue
    
    self.Result = Result
    return [Eval,[self.Lines[-1].XValue[1],self.Lines[-1].YValue[1]]]


    
            













class LayerProblem:
  #This class is a collection of subproblems that need to be solved in conjunction with one another. 
  #While solving the subproblems is nice, it often yeilds sub-par results
  #This is a generalized TSP where all subproblems are "States" and the lines are cities. 

  #A list of subproblems
  Subproblems = []

  #This is the high of the current layer or Z value for all Lines used
  Height = -1

  #This is the starting position, where the previous layer left off
  StartingPoint = []

  #This is the ending position, where the next layer needs to begin closest to
  EndingPoint = []

  Result = ""


  def __init__(self,l_Subproblems,l_Height):
    self.SubProblems = l_Subproblems
    self.Height = l_Height
  
  def Dump(self):
    print()
    print("LAYER DUMP")
    print("Z = ", self.Height)
    print()
    for problem in self.SubProblems:
      problem.Dump()
  
  def Solve(self, FirstLayer = False, EVal = 0):
    pointCloud = []
    start = 1
    
    
    print("Optimizing Layer Z = {}".format(self.Height))

    for Problem in self.SubProblems:
      pointCloud = pointCloud + Problem.GetLines()

    if(len(pointCloud) == 0):
      #There are no lines to optimize, this layer is already as good as it can get!
      self.EndingPoint = self.SubProblems[-1].Travel()
      return self.EndingPoint

    #If this is the first layer we just assume that the first line cannot be changed
    if (FirstLayer):
      start = 0
      self.StartingPoint = [pointCloud[1].XValue[0],pointCloud[1].YValue[0],EVal]
    else:
      #this is harder, we have to find the closest point in the entire layer!
      minimum = [999999,0]
      it = 0
      for point in pointCloud:
        if(GetDistance( point.XValue[0], point.YValue[0],self.StartingPoint[0],self.StartingPoint[1]) < minimum[0]):
          minimum = [GetDistance( point.XValue[0], point.YValue[0],self.StartingPoint[0],self.StartingPoint[1]), it*2]
          #print(minimum)
        if(GetDistance( point.XValue[1], point.YValue[1],self.StartingPoint[0],self.StartingPoint[1]) < minimum[0]):
          minimum = [GetDistance( point.XValue[1], point.YValue[1],self.StartingPoint[0],self.StartingPoint[1]), it*2 + 1]
          #print(minimum)
        it = it + 1
      it = minimum[1]
      it = int(it)
      #print(it)
      #print(pointCloud)
      self.StartingPoint = [pointCloud[int(it/2)].XValue[0],pointCloud[int(it/2)].YValue[0],EVal] if it % 2 == 0 else [pointCloud[int((it-1 )/2)].XValue[1],pointCloud[int((it - 1)/2)].YValue[1],EVal]
      start = it
    
    #Now to make the matrix to actually solve this
    print("Calculating Distance Matrix")
    DistanceMatrix = []
    count = 0
    for Problem in self.SubProblems:
      Size = Problem.GetSize()
      
      for i in range(Size):
        if (i%10 + 1) == 1:
          print("Element {} in problemset {} added".format(i + 1,self.SubProblems.index(Problem)))
        Elem = []
        Elem2 = []
        for n in range(len(DistanceMatrix)):
          if n < count*2:
            if n % 2 == 0:
              DistanceMatrix[n].append(GetDistance( pointCloud[int(count + i)].XValue[0],pointCloud[int(count + i)].YValue[0], pointCloud[int(n/2)].XValue[0],pointCloud[int(n/2)].YValue[0]))
              Elem.append(DistanceMatrix[n][-1])
              DistanceMatrix[n].append(GetDistance( pointCloud[count + i].XValue[1],pointCloud[count + i].YValue[1], pointCloud[int(n/2)].XValue[0],pointCloud[int(n/2)].YValue[0]))
              Elem2.append(DistanceMatrix[n][-1])
            else:
              DistanceMatrix[n].append(GetDistance( pointCloud[count + i].XValue[0],pointCloud[count + i].YValue[0], pointCloud[int((n - 1)/2)].XValue[1],pointCloud[int((n - 1)/2)].YValue[1]))
              Elem.append(DistanceMatrix[n][-1])
              DistanceMatrix[n].append(GetDistance( pointCloud[count + i].XValue[1],pointCloud[count + i].YValue[1], pointCloud[int((n - 1)/2)].XValue[1],pointCloud[int((n - 1)/2)].YValue[1]))
              Elem2.append(DistanceMatrix[n][-1])
          else:
            DistanceMatrix[n].append(0)
            DistanceMatrix[n].append(0)
            Elem.append(0)
            Elem2.append(0)
        Elem.append(0)
        Elem.append(0)
        Elem2.append(0)
        Elem2.append(0)
        DistanceMatrix.append(Elem)
        DistanceMatrix.append(Elem2)
      count = count + Size
    

    print("Cleaning Distance Matrix")

    #These next steps are for cleaning up the Distance matrix
    #First we have to know all the weights of every possible edge in the matrix
    Total = sum([sum(x) for x in DistanceMatrix])

    #Now we setup buckets for the search algorithm
    Buckets = [x.GetSize()*2 for x in self.SubProblems]

    #Next we must add the total to all non free moves
    for i in range(len(DistanceMatrix)):
      for v in range(len(DistanceMatrix)):
        if BucketSearch(v,Buckets) != BucketSearch(i,Buckets):
          DistanceMatrix[i][v] = DistanceMatrix[i][v] + Total 

    #This adds a free node at the end, the point is that we can end there, but just ignore it afterwards
    DistanceMatrix.append([0 for x in DistanceMatrix])
    for x in DistanceMatrix:
      x.append(0)

    print("Solving Problem")
    #[print(x) for x in DistanceMatrix]

    Assignment = CalcTSP(DistanceMatrix,start,len(DistanceMatrix) - 1)
    #Assignment = solve_tsp(DistanceMatrix, endpoints = (start,len(DistanceMatrix) - 1))

    lastNode = -1
    SubProblemOrder = []
    NodeStarts = []
    NodeEnds = []
    for Node in Assignment:
      if BucketSearch(Node,Buckets) != lastNode:
        lastNode = BucketSearch(Node,Buckets)
        NodeStarts.append(Node)
        if(Assignment.index(Node) != 0):
          NodeEnds.append(Assignment[Assignment.index(Node) - 1])
        SubProblemOrder.append(lastNode) 
    NodeEnds.append(Assignment[-1])

    Offset = []
    for x in SubProblemOrder:
      Offset.append(sum(Buckets[0:x]))

    for problem in range(len(SubProblemOrder) - 1):
      self.SubProblems[SubProblemOrder[problem]].Solve(NodeStarts[problem] - Offset[problem],NodeEnds[problem] - Offset[problem])

    self.SubProblems[SubProblemOrder[-1]].Solve(NodeStarts[-1] - Offset[-1],NodeEnds[-1] - Offset[-1], True)

    for x in self.SubProblems:
      if(x.GetSize() < 1):
        SubProblemOrder.insert(self.SubProblems.index(x),self.SubProblems.index(x))

    newList = []
    for x in SubProblemOrder:
      newList.append(self.SubProblems[x])
    
    for i in SubProblemOrder:
      if SubProblemOrder.count(i) > 1:
        print("ERROR FOUND!")
        input()

    for i in newList:
      if newList.count(i) > 1:
        print("ERROR FOUND!")
        input()
    self.SubProblems = newList

    #print(SubProblemOrder)
    self.EndingPoint = self.SubProblems[-1].Travel()
  
  def Clean(self,FirstLayer = False):
    #print("Z = {}".format(self.Height))
    for problem in self.SubProblems:
      problem.Clean(FirstLayer)

  def Stitch(self, Eval, CurrentPos):
    Res = [Eval,CurrentPos]
    for problem in self.SubProblems:
      Res = problem.Stitch(Res[0],Res[1],self.Height)
      self.Result = self.Result + problem.Result

    return [Res[0],Res[1]]
    

       



      








      

    


class FDMGCodeModel:
  #This class just stores all the Layer Problems and attempts to stitch them all back together
  Layers = []

  Result = ""
  
  def __init__(self, l_Layers):
    self.Layers = l_Layers
  
  
  def Dump(self):
    print("Starting whole model dump")
    print("MODEL DUMP")
    for layer in self.Layers:
      layer.Dump()

  def Solve(self):
    for layer in self.Layers:
      if self.Layers.index(layer) == 0:
        layer.Solve(FirstLayer = True)
      else:
        layer.StartingPoint = self.Layers[self.Layers.index(layer) - 1].EndingPoint  
        layer.Solve(FirstLayer = False)
    #Make sure themodel is clean
    for layer in self.Layers:
      if layer.Height == 0:
        layer.Clean(True)
      else:
        layer.Clean()


  def Stitch(self):
    #print("Layer")
    CurrentEValue = 0
    CurrentPos = self.Layers[0].StartingPoint
    #print(CurrentPos)
    for layer in self.Layers:
      if (layer.StartingPoint != []):
        CurrentPos = [layer.StartingPoint[0], layer.StartingPoint[1]]
        #print(CurrentPos)
        break

    Res = [CurrentEValue,CurrentPos]

    self.Result = ";This GCode has been Optimized with a TSP Optimization script"
    self.Result = self.Result + "Date Optimized: {}".format(datetime.time.strftime("%m/%d/%Y"))
    self.Result = self.Result + "Script Version: {}".format(Version)
    for layer in self.Layers:
      Res = layer.Stitch(Res[0],Res[1])
      self.Result = self.Result + layer.Result
    

    #print(self.Result)
    return self.Result




#####################################################################
                    #Section 4: Program#
#####################################################################

#Parser for scanning the script, creates a valid GCode Model from Gcode
#Code - The entire GCode as a list, delimited by newlines
def GCodeParser(Code):

  #Note: For our purposes only lines of code that contain a G1 X Y and E can be counted as
  #Starting code for a subproblem
  
  #Temp lists for storing information
  Lines = []
  SubProblems = []
  Layers = []

  #Marks if there is a problem
  NoCurrentProblem = True
  Preamble = ""

  count = 0
  for line in Code:
    count = count + 1
    if count % 100 == 0:
      print("{} out of {} lines processed".format(count,len(Code))) 


    if(line == "G0 F4285.7 X122.547 Y116.014"):
      pass
      #print(Line)

    if(GetZValue(line,True) != -1):
      #A new layer has begun!
      
      Start = [0,0,0]
      if(len(SubProblems) != 0):
        Start = SubProblems[-1].Travel()
      problem = SubProblem(Lines,Preamble,l_Start = Start,l_End = [GCodeParser.Position[x] for x in range(3)])
      
      Preamble = ""
      SubProblems.append(problem)
      Lines = []
      newLayer = LayerProblem(SubProblems,GCodeParser.Position[2])
      SubProblems = []
      Layers.append(newLayer)
      NoCurrentProblem = True


    #If there is no current subproblem, find one
    if(NoCurrentProblem):
      if(GetGValue(line) == 1 and
         GetXValue(line,True) != -1 and
         GetYValue(line,True) != -1 and
         GetEValue(line,True) != -1):
         NoCurrentProblem = False

         #Travel Line and add it
         Currentline = Line(line,GCodeParser.Position,10)
         Lines.append(Currentline)
      else:

        #otherwise add it to the preamble for this subproblem and travel down it
        Preamble = Preamble + line + '\n'

    else:
      #Make sure that the next line is valid to add to the subproblem
      if(GetXValue(line,True) != -1 and
         GetYValue(line,True) != -1 and
         IsG01(line) == True):
         #Line is valid add if needed and travel
         Currentline = Line(line,GCodeParser.Position,10)
         if(GetGValue(line) == 1):
           Lines.append(Currentline)
      else:
        #Line is invalid start new subproblem and cap off the last one
         Start = [0,0,0]
         if(len(SubProblems) != 0):
          Start = SubProblems[-1].Travel()
         problem = SubProblem(Lines,Preamble,l_Start=Start,l_End = [GCodeParser.Position[x] for x in range(3)])

         Preamble = line + '\n'
         Lines = []
         NoCurrentProblem = True
         SubProblems.append(problem)

    if(IsG01(line)):
      Currentline = Line(line,GCodeParser.Position,10)
      GCodeParser.Position = [Currentline.XValue[1],Currentline.YValue[1],Currentline.ZValue[1],GCodeParser.Position[3] +  (  Currentline.EValue if GetGValue(line) == 1 else 0)]
      print(GCodeParser.Position)
  Start = [0,0,0]
  if(len(SubProblems) != 0):
    Start = SubProblems[-1].Travel()
  problem = SubProblem(Lines,Preamble,l_Start=Start,l_End = [GCodeParser.Position[x] for x in range(3)])
  SubProblems.append(problem)
  newLayer = LayerProblem(SubProblems,GCodeParser.Position[2])
  Layers.append(newLayer)
  newModel = FDMGCodeModel(Layers)
  return newModel
GCodeParser.Position = [0,0,0,0] #Of the format X,Y,Z,E

#The main program, all the magic happens here. 
def main(filePath):
  data = ""
  with open(filePath,'r') as file:
    data = file.read()
    data = data.split('\n')

  print("Creating The Model")
  print()
  Model = GCodeParser(data)
  print("Model Created")
  print()
  print("Solving Model")
  Model.Solve()
  print("Cleaning")
  print("Stitching Model")
  return Model.Stitch()
  
#For neatness
if __name__ ==  "__main__":
  if(len(sys.argv) > 1):
    filename = sys.argv[1]
  else:
    filename = "UM3E_20mmTestCube_repaired.gcode"
  OptimizedGcode = main(filename)

  with open(filename.replace(".g","optimized.g"),"w") as file:
    file.write(OptimizedGcode)
