import pandas as pd
import numpy as np
import os

blocks = []

class Block:
    def __init__(self, pulser):
        self.pulser = pulser
        self.IDs = []
        self.dt = []

    def AddData(self, ID, dt):
        self.IDs.append(ID)
        self.dt.append(dt)

    def __str__(self):
        text = "PUlSER: " + str(self.pulser) + "\n"
        for ID, deltaT in zip(self.IDs, self.dt):
            text += "ID: " + str(ID) + " / " + str(deltaT) + "\n"

        return text


def returnBlock(ID):
    global blocks
    for block in blocks:
        if block.pulser == ID:
            return block

    new_block = Block(ID)
    blocks.append(new_block)
    return new_block

def lineToArray(line):
    cur_number = ""
    array = []
    for char in line:
        if char == " ":
            if cur_number != "":
                array.append(float(cur_number))
                cur_number = ""
        else:
            cur_number += char

    return array

def read_AST(file_name):
    cwd = os.getcwd()
    filePath = cwd + "\\DadosExperimentais\\Arquivos\\" + file_name + ".txt"
    file = open(filePath, "r")

    flag = True
    active = False
    while flag:
        line = file.readline()
        if line[:6] == "PULSER":
            ID = int(file.readline())
            file.readline()
            current_block = returnBlock(ID)
            active = True
        elif line == "":
            flag = False
        elif line == ' \n':
            pass
        else:
            if active:
                data = lineToArray(line)
                current_block.AddData(data[0], data[2])

    return blocks

"""
blocks = read_AST("AST")
for block in blocks:
    print(block)
"""