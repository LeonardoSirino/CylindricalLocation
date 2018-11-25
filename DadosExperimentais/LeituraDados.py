import os

blocks = []
cur_block = None
max_dt = 1500


class Block:
    def __init__(self, pulser):
        self.pulser = pulser - 1
        self.IDs = []
        self.dt = []
        self.X = 0
        self.Y = 0
        self.tf = 0
        self.t0 = 0

    def AddData(self, ID, dt):
        self.IDs.append(ID - 1)
        self.dt.append(dt)

    def Newhit(self, Ch, t):
        global max_dt

        if self.t0 == 0:
            self.t0 = t
            self.tf = t

        dt = t - self.t0

        try:
            self.IDs.index(Ch - 1)
            found = True
        except ValueError:
            found = False

        if not found and (t - self.tf) <= max_dt:
            self.IDs.append(Ch - 1)
            self.dt.append(dt)
            self.tf = t

        if (t - self.tf) > max_dt:
            self.CloseBlock()
            
            global cur_block
            cur_block.Newhit(Ch, t)


    def CloseBlock(self):
        global cur_block
        global blocks

        if len(cur_block.IDs) >= 3:
            blocks.append(cur_block)

        new_block = Block(-1)
        new_block.SetCoord(self.X, self.Y)
        
        cur_block = new_block

    def SetCoord(self, x, y):
        self.X = x
        self.Y = y

    def __str__(self):
        text = "PUlSER: " + str(self.pulser) + " --- X: " + str(round(self.X, 2)) + " Y: " + str(round(self.Y, 2)) + "\n"
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

    if cur_number != "":
        array.append(float(cur_number))

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
                if data[3] > 0:  # Filtro de amplitude
                    current_block.AddData(data[0], data[2])

    return blocks


def read_LineDisplayGroup(file_name):
    cwd = os.getcwd()
    filePath = cwd + "\\DadosExperimentais\\Arquivos\\" + file_name + ".txt"
    file = open(filePath, "r")

    flag = True
    k = 0
    while flag:
        line = file.readline()
        if line[:5] == "* Gp#":
            current_block = returnBlock(k)
            ki = line.find("[")
            kn = line.find("]")
            temp = line[(ki + 1):kn]
            channels = temp.split(",")
            channels = [int(x) for x in channels]

            ki = line.find("[", kn)
            kn = line.find("]", kn + 1)
            temp = line[(ki + 1):kn]
            dts = temp.split(",")
            dts = [int(x) for x in dts]
            dts = [0] + dts

            for (ch, dt) in zip(channels, dts):
                current_block.AddData(ch, dt)

            ki = line.find("x =")
            kn = line.find(",", ki)
            x = float(line[(ki + 3):kn])
            ki = line.find("y =")
            kn = line.find("dT", ki)
            y = float(line[(ki + 3):kn])

            current_block.SetCoord(x, y)
            k += 1

        elif line == "":
            flag = False

    return blocks


def read_LineDisplay(file_name, y):
    global cur_block
    cwd = os.getcwd()
    filePath = cwd + "\\DadosExperimentais\\Arquivos\\" + file_name + ".txt"
    file = open(filePath, "r")
    x = -250
    dx = 250

    flag = True
    k = 0
    while flag:
        line = file.readline()
        if line[:3] == "  1":
            data = lineToArray(line)
            t = data[1] * 1E6
            Ch = data[2]
            cur_block.Newhit(Ch, t)
        elif line[:3] == "128":
            x += dx
            if cur_block == None:
                cur_block = Block(-1)
                cur_block.SetCoord(x, y)
            else:
                cur_block.CloseBlock()
            
            cur_block.SetCoord(x, y)
        elif line == "":
            flag = False
        else:
            pass

    global blocks
    if len(cur_block.IDs) >= 3:
        blocks.append(cur_block)

    return blocks

"""
blocks = read_LineDisplay("Linha1_line_display", 0)
for block in blocks:
    print(block)
    pass

print(str(len(blocks)) + " blocos")
"""