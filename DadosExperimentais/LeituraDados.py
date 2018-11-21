import os

blocks = []


class Block:
    def __init__(self, pulser):
        self.pulser = pulser - 1
        self.IDs = []
        self.dt = []
        self.X = 0
        self.Y = 0

    def AddData(self, ID, dt):
        self.IDs.append(ID - 1)
        self.dt.append(dt)

    def SetCoord(self, x, y):
        self.X = x
        self.Y = y

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
                if data[3] > 0:  # Filtro de amplitude
                    current_block.AddData(data[0], data[2])

    return blocks


def read_LineDisplay(file_name):
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


"""
blocks = read_AST("AST_Samos")
for block in blocks:
    print(block)
    pass

"""