
from PySteppables import *
import CompuCell
import sys
import random
import os.path

INIT_MAX_DIVISIONS = 4.0
KILL_PERIOD = 2000
GLU_ABSORBED_FRAC = 0.3
GROWTH_THRESHOLD = 30.
DEATH_THRESHOLD = 1.
INITIAL_CONCENTRATION = 5.

# define a data structure to hold a tree of 
# stem cells and thier progeny
class CellNode: 
    
    def __init__(self, cell, attr, mcs):
        """
        Creates a new CellNode structure, and
        initializes the fields from a given CC3D cell.
        """
        
        # copy cell id
        self.id = cell.id
        self.creationTime = mcs
        
        # no death time
        self.deathTime = None
        
        # this dictionary has the same keys as the main cell
        # dict, but has values of another dictionary with 
        # keys of mcs, and values of values
        self.attr = {}
                
        for k,v in attr.items():

            self.attr[k] = {mcs : v}
        
        # only set the parent id if this cell has a parent
        if attr["Parent_ID"] != cell.id:
            self.parentId = attr["Parent_ID"]
        else:
            self.parentId = None
        
        # no chidren yet.
        self.childIds = []


# a single global dictionary to hold all of the nodes
# this is indexed by cell id.
nodes = {}


# A function to append newly created nodes to our
# tree
def AppendTree(steppable, cell, mcs):
    
    # grab the dict of attributes
    attr =  steppable.getDictionaryAttribute(cell)
    
    # make a new cell that we will add
    c = CellNode(cell, attr, mcs)
   
    # add to the global list of nodes
    nodes[c.id] = c
    
    # if the node has a parent, add this new node to that 
    # parent nodes list of child nodes. 
    # only cells created through mitosis have parent cells, 
    # parent cells have different ids than the current cell
    if c.parentId != None:
        parent = nodes[c.parentId]
        parent.childIds.append(c.id)

   
    
def KillNode(cell, mcs):
    """
    Find a node, and set its death time
    """
    if nodes.has_key(cell.id):
        node = nodes[cell.id]
        node.deathTime = mcs
    
    
    


# set the path of the file where you want to write the output to, 
# the '~' character gets expanded to your home directory with
# the os.path.expanduser function, this works on both 
# unix and windows
OUTPUT_FILE_NAME = os.path.expanduser("~/my_outputfile.txt")

# handle to an global output file, initilize in the 
# initializer steppable
outputFile = None

def initOutputFile(fname = OUTPUT_FILE_NAME):
    """
    creates and initialized an output file, and writes a header.
    """
    print("creating output, name=" + fname)
    
    global outputFile
    outputFile = open(fname, "w")
    
    s = "MCS"
    s += ", "
    s += "EVENT_ID"
    s += ", "
    s += "CELL_TYPE"
    s += ","
    s += "CELL_ID"
    s += ", "
    s += "Parent_ID"
    s += ", "
    s += "Last_Stem_Cell_ID"
    s += ", "
    s += "P_STEM"
    s += ", "
    s += "MAX_DIVISIONS"
    s += ", "
    s += "CURRENT_DIVISIONS"
    s += "\n"
    
    outputFile.write(s)
    outputFile.flush()
    
    

def writeCell(steppable, eventId, mcs, cell):
    """ 
    Writes the state of a cell as a record in an output file.
    """
    global outputFile
    
    dict =  steppable.getDictionaryAttribute(cell)
   
    s = str(mcs)
    s += ", "
    s += str(eventId)
    s += ", "
    s += str(cell.type)
    s += ", "
    s += str(cell.id)
    s += ","
    s += str(dict["Parent_ID"])
    s += ", "
    s += str(dict["Last_Stem_Cell_ID"])
    s += ", "
    s += str(dict["P_STEM"])
    s += ", "
    s += str(dict["MAX_DIVISIONS"])
    s += ", "
    s += str(dict["CURRENT_DIVISIONS"])
    s += "\n"
   
    outputFile.write(s)
    outputFile.flush()
            
   
   

from PySteppablesExamples import MitosisSteppableBase
            

class ConstraintInitializerSteppable(SteppableBasePy):
    """
    The initializer class. 
    
    This should contain logic to set the initial values of all of the 
    initial cells. 
    """
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        field=self.getConcentrationField("Glu")
        field[:,:,:] = INITIAL_CONCENTRATION
        
        
        
        # create the output file, 
        # can optionally choose differnt name here.
        initOutputFile()
        
        for cell in self.cellList:
            cell.targetVolume=25.
            cell.lambdaVolume=40.0
            # access/modification of a dictionary attached to cell - 
            # make sure to decalare in main script that you will use such attribute
            cellDict=self.getDictionaryAttribute(cell)
            cellDict["P_STEM"]=0.10
            cellDict["MAX_DIVISIONS"]=INIT_MAX_DIVISIONS
            cellDict["CURRENT_DIVISIONS"]=0
            cellDict["Parent_ID"]=cell.id
            cellDict["Last_Stem_Cell_ID"]=cell.id
            
            writeCell(self, "INIT", 0, cell)
            AppendTree(self, cell, 0)


class GrowthSteppable(SteppableBasePy):
    """
    Handle the growth behavior. 
    
    Only increase the target volume if it is not zero. A zero target volume
    means that somethign else signaled this cell to die. 
    """
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        
        gluField=self.getConcentrationField("Glu")
        
        
        for cell in self.cellList:
            
            if cell.type != 0 and cell.targetVolume > 0:
                x=int(cell.xCOM)
                y=int(cell.yCOM)
                z=int(cell.zCOM)
            
                gluValue=gluField[x,y,z]
            
                gluAbs = GLU_ABSORBED_FRAC * gluValue
            
                gluField[x,y,z] = gluValue - gluAbs
                
            
                if (cell.targetVolume > 0) and (cell.targetVolume-cell.volume < 5):
                    cell.targetVolume+=1. * gluAbs / (GROWTH_THRESHOLD + gluAbs)        
  

class MitosisSteppable(MitosisSteppableBase):
    """
    Handle the motisis behavior.
    
    
    Stem cells may create either somatic of stem cell children, but
    somatic cells can only create other somatic cells. 
    
    Somatic cells have a finite number of cell divisions, but stem cells
    can divide an unlimited number of times. 
    """
    
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        self.mcs = 0
    
    def step(self,mcs):
        print "INSIDE MITOSIS STEPPABLE"
        self.mcs = mcs
        
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume > 50 and cell.type != 0:
                
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell
        
        parentCell.targetVolume=25
        childCell.targetVolume=25
        childCell.lambdaVolume=parentCell.lambdaVolume
        
        pdict=self.getDictionaryAttribute(parentCell)
        cdict=self.getDictionaryAttribute(childCell)
        
        print("pdict: " + str(pdict))
        print("cdict: " + str(cdict))
        
        # copy parent dict attributes to child
        
        
        
        if random.random() < 0.01:
#       Mutation
            x=random.gauss(0.0, .5)
            
            if x>0:
                pdict["MAX_DIVISIONS"]+=x
            else:
                pdict["MAX_DIVISIONS"]=pdict["MAX_DIVISIONS"]/(1 - x/INIT_MAX_DIVISIONS)
                
            writeCell(self, "MUTATION_DIVISIONS", self.mcs, parentCell)
            
            # add a new max divisions to the tree
            c = nodes[parentCell.id]
            c.attr["MAX_DIVISIONS"][self.mcs] = pdict["MAX_DIVISIONS"]


        if random.random() < 0.01:
#       Mutation
            x=random.gauss(0.0, .05)
            
            if pdict["P_STEM"]+x < 1 and pdict["P_STEM"]+x > 0:
                pdict["P_STEM"]+=x 
                
            writeCell(self, "MUTATION_P_STEM", self.mcs, parentCell)
            
            # add a new pstem mutation to the tree
            c = nodes[parentCell.id]
            c.attr["P_STEM"][self.mcs] = pdict["P_STEM"]
        
        
        for k,v in pdict.items():
            cdict[k] = v
            
            
        if parentCell.type==1:
            # parent cell was stem
            if random.random() < pdict["P_STEM"]:
                # set child type to stem
                childCell.type = 1
            else:
                # set child type to differentiated
                childCell.type=2
        else:
            # parent cell was differetiated, child type is always 
            # differentiated
            childCell.type=2
            
            pdivisions = pdict["CURRENT_DIVISIONS"]
            
            maxDivisions = pdict["MAX_DIVISIONS"]
            
            print("parent divisions: " + str(pdivisions))
            
            cdict["CURRENT_DIVISIONS"] = pdivisions + 1
            pdict["CURRENT_DIVISIONS"] = pdivisions + 1
    
        
       
            
            if pdivisions >= maxDivisions:
                # kill the cell
                parentCell.targetVolume=0
                childCell.targetVolume=0
                
                writeCell(self, "MITOSIS_DEATH", self.mcs, parentCell)
                writeCell(self, "MITOSIS_DEATH", self.mcs, childCell)
                
                # only kill the parent node, as this one exists in the
                # tree, and the child node has not been inserted yet
                KillNode(parentCell, self.mcs)
                
                # both parent and child are dead, so no need for
                # further processing
                return

                
        cdict["Parent_ID"]=parentCell.id
        #Keep track of child cell's immediate parent
        #If parent stem cell keep tack of last setm cell ancestor
        if parentCell.type == 1:
            cdict["Last_Stem_Cell_ID"]=parentCell.id
        else:
            cdict["Last_Stem_Cell_ID"]=pdict["Last_Stem_Cell_ID"]
        
        # write the newly created child cell out to the file
        writeCell(self, "MITOSIS", self.mcs, childCell)
        AppendTree(self, childCell, self.mcs)
            
            
class DeathSteppable(SteppableBasePy):
    """
    All cells have a probablity of random death.
    """
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def step(self,mcs):
        gluField=self.getConcentrationField("Glu")
        
        for cell in self.cellList:
            
            
            if cell.type != 0:
                x=int(cell.xCOM)
                y=int(cell.yCOM)
                z=int(cell.zCOM)
                gluValue=gluField[x,y,z]   
               
                if gluValue < DEATH_THRESHOLD:
                    writeCell(self, "DEATH", mcs, cell)
                    KillNode(cell, mcs)
                    
                    cell.targetVolume=0
                    cell.lambdaVolume=0
                    cell.type = 0
                    
            
        
        if mcs % KILL_PERIOD == 0:
        
            for cell in self.cellList:
                # pick rate of cell death, 1/1000
                if cell.xCOM > 100:
                    #if cell.type == 2 or random.random() < 0.5 :
                    
                    #Need to write to log file 
                    writeCell(self, "DEATH", mcs, cell)
                    KillNode(cell, mcs)
                    
                    cell.targetVolume=0
                    cell.lambdaVolume=0
                    cell.type = 0

 


class PlotSteppable(SteppableBasePy):
    """
    A stepable to create plots. 
    
    Plots must be created in the 'start' function, but data is added to the
    plots in the 'step' function. 
    """
    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        
    def start(self):
        
        # create a plot window, this code is added by the 
        # CC3DPython -> Scientific Plots -> Setup menu item
        self.pW=self.addNewPlotWindow(_title='Statistics',
            _xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Variables')
            
        self.pW.addPlot('Stem Count',_style='Lines',_color='red', _size=2)
        
        # create 3 histogram plot windows, this code is added by the
        # CC3DPython -> Scientific Plots -> Add Histogram Plot (start fcn) menu item
        
        #initialize setting for Histogram
        self.maxDivPlt=self.addNewPlotWindow(_title='Max Divisions',_xAxisTitle='max divisions',_yAxisTitle='cell count')
        self.maxDivPlt.addHistogramPlot(_plotName='Max Hist',_color='green',_alpha=100)# _alpha is transparency 0 is transparent, 255 is opaque        
    
        #initialize setting for Histogram
        self.currDivPlt=self.addNewPlotWindow(_title='Current Divisions',_xAxisTitle='current divisions',_yAxisTitle='cell count')
        self.currDivPlt.addHistogramPlot(_plotName='Curr Hist',_color='blue',_alpha=100)
        
        #initialize setting for Histogram
        self.pstemPlt=self.addNewPlotWindow(_title='Stem Child Prob',_xAxisTitle='Prob of Stem Child',_yAxisTitle='cell count')
        self.pstemPlt.addHistogramPlot(_plotName='Stem Prob',_color='red',_alpha=100)

  

    def step(self,mcs):
        n_type1 = 0
        n_type2 = 0
        divisions = 0.0
        
        if mcs < 500 :
            return
        
        for cell in self.cellList:
            if cell.type == 1:
                n_type1 += 1
            elif cell.type == 2:
                n_type2 += 1
                
            
                
        if n_type2 != 0 :
             self.pW.addDataPoint("Stem Count",mcs,float(n_type1)/n_type2)  
      

        # add data to the histogram plots, needs numpy to generate
        # histogram data
        import numpy
 
        # make a list of all of the max divisions attributes, but
        # we only want to get the max divisions for somatic cells
        max_divs = []
        for cell in self.cellList:
            if cell.type == 2:
                max_divs.append(self.getDictionaryAttribute(cell)["MAX_DIVISIONS"])

        # make a list of all of the current divisions, again, only
        # interested in the ones for somatic cells
        curr_divs = []
        for cell in self.cellList:
            if cell.type == 2:
                curr_divs.append(self.getDictionaryAttribute(cell)["CURRENT_DIVISIONS"])

            
        # make a list of the probability to make a stem child, here
        # we're only interested in the stem cells
        p_stem = []
        for cell in self.cellList:
            if cell.type == 1:
                p_stem.append(self.getDictionaryAttribute(cell)["P_STEM"])

        # use the numpy histogram function to bin the data
        (max_n, max_bins) = numpy.histogram(max_divs, bins=15)    
        (curr_n, curr_bins) = numpy.histogram(curr_divs, bins=15)  
        (pstem_n, pstem_bins) = numpy.histogram(p_stem, bins=15) 
  
                
        self.maxDivPlt.addHistPlotData('Max Hist',max_n, max_bins)
        self.currDivPlt.addHistPlotData('Curr Hist',curr_n, curr_bins)
        self.pstemPlt.addHistPlotData('Stem Prob',pstem_n, pstem_bins)
        
        # display the regular plot window
        self.pW.showAllPlots() 
        
        
        # need to go through all the histogram plot windows, and 
        # tell them to re-display
        self.currDivPlt.showAllHistPlots()
        self.maxDivPlt.showAllHistPlots()
        self.pstemPlt.showAllHistPlots()
   
   
    
