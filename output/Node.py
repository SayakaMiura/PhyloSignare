import os
import sys
from pathlib import Path
class Node:
        tree_age = 1
    
        def __init__(self, name, number_of_mutations):
            self.parent = None
            self.age = -1
            self.isNumbered = False
            self.children = []
            self.name = name
            self.number_of_mutations = number_of_mutations
            self.isChild = False
            self.numbered_name = "Normal"
    
        def __str__(self):
            return ', label = "' + 'B' + str(self.age) + ': ' + self.number_of_mutations + '"'
    
        # Looks for phylosigfinder for appropriate node and determines how many presences there are
        @staticmethod
        def read(node, path):
    
            # Opens phylosigfinder file
    
            if Path(path / Path(node.name + '_MutCount_PhyloSigFinder.txt')).exists():    
                phylosigfinder = open(path / Path(node.name + '_MutCount_PhyloSigFinder.txt'))
                phylosigfinder = phylosigfinder.readlines()
                presences = []
            
                        # Remove the first line
                phylosigfinder.pop(0)
            
                        # Split each line into an array
                for i in range(len(phylosigfinder)):
                    phylosigfinder[i] = phylosigfinder[i].split()
            
                        # Detects which signatures are present
                for line in phylosigfinder:
                    if line[1] == '1':
                        presences.append(line[0])
                return presences
            else:
                return []
    
        def add_child(self, node):
            node.parent = self
            self.children.append(node)
    
        def set_age(self):
            self.age = Node.tree_age
            Node.tree_age += 1
            self.isNumbered = True
            print(self.name + " age: " + str(self.age))
    
    
    # Modify tree and node

    

    
    