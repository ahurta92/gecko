from collections import OrderedDict

import json
import pandas as pd
import re


class daltonToJson:
    def __init__(self):
        self.outfile = None
        self.quad_response = None
        self.dipole_dict = None
        self.polar_dict = None
        self.calcSetup = {}
        self.calcRes = {}
        self.calcTask = {}
        self.simulationTime = {"cpuTime": 0.0, "wallTime": 0.0, "units": "second"}
        self.symmetry = {}
        self.calculations = []
        self.taskNumber = 1
        self.setupCount = 0
        self.subTask = False

        self.xx = 0
        self.xy = 0
        self.xz = 0
        self.yx = 0
        self.yy = 0
        self.yz = 0
        self.zx = 0
        self.zy = 0
        self.zz = 0

    def convert(self, streamIn):
        # Takes in the file and converts to JSON
        # This ordered dictionary allows me to read in one module at a time
        my_tasks = OrderedDict(
            [
                (
                    "a direct, restricted step, second order MCSCF program  ",
                    self.readScfDft,
                ),
                ("Linear Response calculation", self.readResponse),
                ("*                   SUMMARY OF COUPLED CLUSTER CALCULATION                    *",self.readExcitedCC),
                ("Singlet electronic excitation energies", self.readExcited),
                ("Nuclear contribution to dipole moments", self.readDipole),
                ("SYMGRP: Point group information", self.readSymGrp),
            ]
        )
        collectingInput = False
        line = streamIn.readline()
        while line:
            for myKey in my_tasks.keys():
                if line.find(myKey) >= 0:
                    # print(myKey)
                    myIndex = list(my_tasks).index(myKey)
                    if myIndex >= 0:
                        self.calcTask = {}
                        self.calcSetup = {}
                        self.calcRes = {}
                        self.taskNumber += 1
                        self.setupCount += 1
                    my_tasks[myKey](line, streamIn)
            if line.find("Total CPU  time used in DALTON:") >= 0:
                line2 = streamIn.readline()
                self.simulationTime = self.readTiming(line, line2)
                break
            line = streamIn.readline()

        return json.dumps(
            {
                "simulation": {
                    "calculations": self.calculations,
                    "simulationTime": self.simulationTime,
                    "symmetry": self.symmetry,
                }
            },
            indent=2,
            separators=(",", ": "),
            ensure_ascii=False,
        )

    def readScfDft(self, line, streamIn):
        # print('READING SCF')
        self.calcSetup["numberOfElectrons"] = 0
        self.calcSetup["molecularSpinMultiplicity"] = 1

        # for each one of these functions I read in a line and capture the value.
        # Value gets places in correct dictionary
        def closedShell(line):
            self.calcSetup["numberOfElectrons"] = int(line.split()[6])

        def doCharge(line):
            self.calcSetup["charge"] = int(float(line.split()[6]))

        def waveFuncSCF(line):
            if line.split()[5] == "HF":
                self.calcSetup["waveFunctionType"] = "HF"
                self.calcSetup["waveFunctionTheory"] = "Hartree-Fock"
            elif line.split()[5] == "KS-DFT":
                self.calcSetup["waveFunctionType"] = "KS-DFT"
                self.calcSetup["waveFunctionTheory"] = "DFT"

        def totalEn(line):
            self.calcRes["totalEnergy"] = {
                "value": float(line.split(":")[1]),
                "units": "Hartree",
            }

        def totalELectronEn(line):
            self.calcRes["electronEnergy"] = {
                "value": float(line.split(":")[1]),
                "units": "Hartree",
            }

        def nucEn(line):
            self.calcRes["nuclearRepulsionEnergy"] = {
                "value": float(line.split(":")[1]),
                "units": "Hartree",
            }

        self.calcTask["calculationType"] = "energyCalculation"
        # Here I can define what I use to match what I want to capture
        # if the line matches to one of these inputs we will take the line
        # and capture value and save in dictionary using above functions
        scfInp = {
            "@    Wave function type": waveFuncSCF,
            "@    Number of closed shell electrons": closedShell,
            "@    Total charge of the molecule": doCharge,
            "@    Final DFT energy": totalEn,
            "@    Final HF energy:": totalEn,
            "@    Nuclear repulsion:": nucEn,
            "@    Electronic energy:": totalELectronEn,
        }
        line = streamIn.readline()
        while line:
            for scfKey in scfInp.keys():
                if line.find(scfKey) >= 0:
                    # print(scfKey)
                    scfInp[scfKey](line)
            if line.find("CPU and wall time for SCF :") >= 0:
                # print(line)
                time = line.split(":")[1]
                time = time.split()
                self.calcTask["calculationTime"] = {
                    "cpuTime": float(time[0]),
                    "wallTime": float(time[1]),
                    "units": "second",
                }
            elif line.find("End of Wave Function Section") >= 0:
                # print(line)
                break
                # print(time)
            line = streamIn.readline()
        self.calcTask["calculationResults"] = self.calcRes
        self.calcTask["calculationSetup"] = self.calcSetup
        self.calculations.append(self.calcTask)

    def readQuadResponse(self, outfile):
        #print("READING QUAD RESPONSE from {}".format(outfile))
        self.outfile = outfile

        with open(self.outfile, "r") as file:
            rows = []
            for line in file:
                pattern = "  Results from quadratic response calculation"
                match = re.search(pattern, line)
                if match:
                    print(line)
                    line = file.readline()
                    line = file.readline()
                    for line in file:
                        if line.startswith("@ B-freq"):
                            # Create a regex pattern that captures three letters (X, Y, or Z) within parentheses,
                            # separated by semicolons, and the floating point numbers in the string
                            pattern = r"B-freq = ([\d\.]+)\s+C-freq = ([\d\.]+)\s+beta\((X|Y|Z);(X|Y|Z),(X|Y|Z)\) = (.*)"

                            # Use the search function from the re module to get the match object
                            match = re.search(pattern, line)

                            if match:
                                (
                                    b_freq,
                                    c_freq,
                                    letter1,
                                    letter2,
                                    letter3,
                                    beta_value,
                                ) = match.groups()
                                # Check if beta_value is another beta function call
                                # not a good idea
                                # if "beta" in beta_value:
                                # If it is, replace it with the original letters
                                # beta_value = f"beta({letter1};{letter2},{letter3})"

                                a_freq = -(float(b_freq) + float(c_freq))
                                # Append the data to the DataFrame
                                rows.append(
                                    pd.Series(
                                        {
                                            "A-freq": a_freq,
                                            "B-freq": float(b_freq),
                                            "C-freq": float(c_freq),
                                            "A": letter1,
                                            "B": letter2,
                                            "C": letter3,
                                            "Beta Value": (
                                                beta_value
                                                if "," in beta_value
                                                else float(beta_value)
                                            ),
                                        }
                                    )
                                )
                        else:
                            break
            try:
                beta_json = pd.concat(rows, axis=1).transpose()
            except ValueError as verror:
                print("ValueError: {}".format(verror))
                print("rows: {}".format(rows))
                beta_json = pd.DataFrame()
            # if beta_json is not empty process

            return beta_json

    def readAlphaFromQUADResponse(self, outfile):
        # print("READING QUAD RESPONSE from {}".format(outfile))
        self.outfile = outfile

        with open(self.outfile, "r") as file:
            rows = []
            for line in file:
                if line.startswith("@ QRLRVE:  <<"):
                    #   print(line)
                    line_split = line.split()
                    #   print(line_split)
                    i = line_split[3][0]
                    j = line_split[5][0]
                    omega = line_split[8].split(")")[0]
                    alpha = line_split[9]
                    #   print(i,j,omega,alpha)
                    rows.append(
                        pd.Series(
                            {"ij": i + j, "omega": float(omega), "alpha": float(alpha)}
                        )
                    )
                    # Check if beta_value is another beta function call
                    # not a good idea
                    # if "beta" in beta_value:
                    # If it is, replace it with the original letters
                    # beta_value = f"beta({letter1};{letter2},{letter3})"

                    # a_freq = -(float(b_freq) + float(c_freq))
                    # Append the data to the DataFrame
                    # rows.append(pd.Series({
                    #    "A-freq": a_freq,
                    #    "B-freq": float(b_freq),
                    #    "C-freq": float(c_freq),
                    #    "A": letter1,
                    #    "B": letter2,
                    #    "C": letter3,
                    #    "Beta Value": beta_value if "," in beta_value else float(
                    #        beta_value)
                    # }))
            try:
                alpha_json = pd.concat(rows, axis=1).transpose()
            except ValueError as verror:
                print("ValueError: {}".format(verror))
                print("rows: {}".format(rows))
                alpha_json = pd.DataFrame()
            # if beta_json is not empty process

            return alpha_json

    def readResponse(self, line, streamIn):
        self.calcTask["calculationType"] = "LinearResponse"

        def frequencies(line, streamIn):
            # print("READING FREQ")

            freq = line.split()[2:]
            num_freq = int(line.split()[0])
            if len(freq) != num_freq:
                nextline = streamIn.readline()
                next_freq = nextline.split()
                freq = freq + next_freq
            self.calcSetup["frequencies"] = []
            self.calcSetup["numFrequencies"] = int(num_freq)
            # empty list of size 3
            self.calcRes["secondOrderProp"] = {}
            self.calcRes["secondOrderProp"]["values"] = []

            self.polar_dict = {
                "xx": [0] * int(num_freq),
                "xy": [0] * int(num_freq),
                "xz": [0] * int(num_freq),
                "yx": [0] * int(num_freq),
                "yy": [0] * int(num_freq),
                "yz": [0] * int(num_freq),
                "zx": [0] * int(num_freq),
                "zy": [0] * int(num_freq),
                "zz": [0] * int(num_freq),
            }

            # need to get the values at frequency
            for f in freq:
                self.calcSetup["frequencies"].append(float(f.replace("D", "E")))

        def secondOrderProp(line, streamIn):
            lineSP = line.split()
            opA = lineSP[2]
            opB = lineSP[4]

            opKey = " ; ".join([opA, opB])

            def grab_val(line):
                return float(line[7])

            def grab_operator(line):
                return " ".join(line)[1:6]

            if opKey == "XDIPLEN ; XDIPLEN":
                self.polar_dict["xx"][self.xx] = grab_val(lineSP)
                self.xx += 1
                # print(opKey,self.xx,lineSP)
            elif opKey == "XDIPLEN ; YDIPLEN":
                self.polar_dict["xy"][self.xy] = grab_val(lineSP)
                self.polar_dict["yx"][self.yx] = self.polar_dict["xy"][self.xy]
                self.xy += 1
                self.yx += 1
            elif opKey == "XDIPLEN ; ZDIPLEN":
                self.polar_dict["xz"][self.xz] = grab_val(lineSP)
                self.polar_dict["zx"][self.zx] = self.polar_dict["xz"][self.xz]
                self.xz += 1
                self.zx += 1
            elif opKey == "YDIPLEN ; YDIPLEN":
                self.polar_dict["yy"][self.yy] = grab_val(lineSP)
                self.yy += 1
                # print(opKey,self.yy)
            elif opKey == "YDIPLEN ; ZDIPLEN":
                self.polar_dict["yz"][self.yz] = grab_val(lineSP)
                self.polar_dict["zy"][self.zy] = self.polar_dict["yz"][self.yz]
                self.yz += 1
                self.zy += 1
            elif opKey == "ZDIPLEN ; ZDIPLEN":
                self.polar_dict["zz"][self.zz] = grab_val(lineSP)
                self.zz += 1
                # print(opKey,self.zz)

        scfInp = OrderedDict(
            [
                ("B-frequencies", frequencies),
                (
                    ">> =",
                    secondOrderProp,
                ),
            ]
        )
        # if the line matches to one of these inputs we will take the line
        line = streamIn.readline()
        while line:
            # SCF converged
            if line.find("Total CPU  time used in RESPONSE:") >= 0:
                line2 = streamIn.readline()
                self.calcTask["calculationTime"] = self.readTaskTimes(line, line2)
                break
            else:
                for scfKey in scfInp.keys():
                    if line.find(scfKey) >= 0:
                        scfInp[scfKey](line, streamIn)
            line = streamIn.readline()
        self.calcTask["calculationResults"] = self.polar_dict
        self.calcTask["calculationSetup"] = self.calcSetup
        self.calculations.append(self.calcTask)
    def readExcitedCC(self, line, streamIn):
        self.calcTask["calculationType"] = "SingletExcitationEnergy"
        # if the line matches to one of these inputs we will take the line

        # skip 8 lines
        for i in range(8):
            line = streamIn.readline()

        # remove any '' from the list
        line = list(filter(None, line.split(" ")))
        self.scf_energy= float(line[-2])
        self.mp2_energy=list(filter(None, streamIn.readline().split(" ")))[-2]
        self.cc2_energy=list(filter(None, streamIn.readline().split(" ")))[-2]
        #print("SCF Energy:", self.scf_energy)
        #print("MP2 Energy:", self.mp2_energy)
        #print("CC2 Energy:", self.cc2_energy)


        for i in range (6):
            line = streamIn.readline()
        line = streamIn.readline()
        #print(line)

        self.calcRes["frequencies"] = []
        self.calcRes["eV"] = []
        self.calcRes["Sym"] = []
        self.calcRes["Mode"] = []

        self.calcRes["cc2_energy"] = self.cc2_energy
        self.calcRes["mp2_energy"] = self.mp2_energy
        self.calcRes["scf_energy"] = self.scf_energy

        while line:
            line=streamIn.readline()

            if line.find("                                                                                ") >= 0:
                print("END OF CC Calculation")
                break
            elif line.find("+-----------------------------------------------------------------------------+") >= 0:
                pass
                #print("We skip this line")
            else:
                line=list(filter(None,line.split(" ")))
                line=list(filter(lambda x: x!="|",line))
                # if at this point len list is > 1 then we have a result that we would like to capture
                if len(line)>1:
                    #print(line)
                    self.calcRes["Sym"].append(str(line[0][1:]))
                    self.calcRes["Mode"].append(int(line[1]))
                    self.calcRes["frequencies"].append(float(line[2]))
                    self.calcRes["eV"].append(float(line[3]))
            

        while line:
            # SCF converged
            line = streamIn.readline()
            if line.find("Total CPU  time used in ABACUS:") >= 0:
                #print("We break out")
                line2 = streamIn.readline()
                self.calcTask["calculationTime"] = self.readTaskTimes(line, line2)
                break
        self.calcTask["calculationResults"] = self.calcRes
        self.calculations.append(self.calcTask)


    def readExcited(self, line, streamIn):
        # print("READING RESPONSE")
        self.calcTask["calculationType"] = "SingletExcitationEnergy"
        # if the line matches to one of these inputs we will take the line
        line = streamIn.readline()
        while line:
            # SCF converged
            if line.find("Total CPU  time used in ABACUS:") >= 0:
                line2 = streamIn.readline()
                self.calcTask["calculationTime"] = self.readTaskTimes(line, line2)
                break
            # find the first line to start reading
            if line.find("==============") >= 0:
                self.calcRes["frequencies"] = []
                self.calcRes["Sym"] = []
                self.calcRes["Mode"] = []
                line = streamIn.readline()
                while line:
                    if line.find("==============") >= 0:
                        break
                    # find the first ---- line
                    if line.find("------------") >= 0:
                        self.calcRes["frequencies"].append([])
                        self.calcRes["Sym"].append([])
                        self.calcRes["Mode"].append([])
                        line = streamIn.readline()
                        while line:
                            if line.find("-----------") >= 0:
                                self.calcRes["frequencies"].append([])
                                self.calcRes["Sym"].append([])
                                self.calcRes["Mode"].append([])
                            elif line.find("========") == -1:
                                self.calcRes["frequencies"][-1].append(
                                    float(line.split()[2])
                                )
                                self.calcRes["Sym"][-1].append(int(line.split()[0]))
                                self.calcRes["Mode"][-1].append(int(line.split()[1]))
                            else:
                                break
                            line = streamIn.readline()
                    if line.find("==============") >= 0:
                        break
                    line = streamIn.readline()
                    # read if not ------
            line = streamIn.readline()
        self.calcTask["calculationResults"] = self.calcRes
        self.calculations.append(self.calcTask)

    def readSymGrp(self, line, streamIn):
        # print("READING SYMGRP")
        # Two options for symmetry... 
        # Either is no symmetry aka C1 
        #  - full_point_group = 'C1'
        #  - Representation point group: C1
        #  - Irrep name for each symmetry: symmetry:  ['A']
        # else there is symmetry and we need to read the point group and irreps

        line = streamIn.readline()
        line = streamIn.readline()
        group_line = streamIn.readline()
        full_point_group = group_line.split(" ")
        # remove all empty strings
        full_point_group = list(filter(None, full_point_group))[-2]
        print("full point group:", full_point_group)
        if(full_point_group == "C1"):
            self.symmetry["point_group"] = full_point_group
            self.symmetry["representation"] = full_point_group
            self.symmetry["irreps"] = ["A"]
            print("Symmetry:", self.symmetry)
        else:
            symrep_line = streamIn.readline()  # Representation
            symrep_line = list(filter(None, symrep_line.split(" ")))[-2]
            print("Representation point group:", full_point_group)
            line = streamIn.readline()  # Representation
            irrep_line = streamIn.readline()  # Irrep name for each symmetry:
            irrep_line = list(filter(None, irrep_line.split(" ")))
            # search for the index of the string symmetry:
            sym_index = irrep_line.index("symmetry:")
            # collect from the index to the end of the list -1
            sym_list = irrep_line[sym_index + 1 : -1]
            # grab every other element starting at 1
            sym_list = sym_list[1::2]
            print("sym_list:", sym_list)
            self.symmetry["point_group"] = full_point_group
            self.symmetry["representation"] = symrep_line
            self.symmetry["irreps"] = sym_list
            print("Symmetry:", self.symmetry)




    def readDipole(self, line, streamIn):
        self.calcTask["calculationType"] = "Dipole"
        self.dipole_dict = {"x": 0, "y": 0, "z": 0}
        # if the line matches to one of these inputs we will take the line
        line = streamIn.readline()
        line = streamIn.readline()
        line = streamIn.readline()
        line = streamIn.readline()
        while line:
            if line.find("Total CPU  time used in HERMIT:") >= 0:
                break
            else:
                line_split = line.split()
                if len(line_split) > 0 and (
                    (line_split[0] == "x")
                    or (line_split[0] == "y")
                    or (line_split[0] == "z")
                ):
                    opKey = line_split[0]
                    dipole_value = line_split[1]
                    if opKey:
                        self.dipole_dict[opKey] = dipole_value
            line = streamIn.readline()

        self.calcTask["calculationResults"] = self.dipole_dict
        self.calculations.append(self.calcTask)

    def readExcited(self, line, streamIn):
        # print("READING RESPONSE")
        self.calcTask["calculationType"] = "SingletExcitationEnergy"
        # if the line matches to one of these inputs we will take the line
        line = streamIn.readline()
        while line:
            # SCF converged
            if line.find("Total CPU  time used in ABACUS:") >= 0:
                line2 = streamIn.readline()
                self.calcTask["calculationTime"] = self.readTaskTimes(line, line2)
                break
            # find the first line to start reading
            if line.find("==============") >= 0:
                self.calcRes["frequencies"] = []
                self.calcRes["Sym"] = []
                self.calcRes["Mode"] = []
                line = streamIn.readline()
                while line:
                    if line.find("==============") >= 0:
                        break
                    # find the first ---- line
                    if line.find("------------") >= 0:
                        self.calcRes["frequencies"].append([])
                        self.calcRes["Sym"].append([])
                        self.calcRes["Mode"].append([])
                        line = streamIn.readline()
                        while line:
                            if line.find("-----------") >= 0:
                                self.calcRes["frequencies"].append([])
                                self.calcRes["Sym"].append([])
                                self.calcRes["Mode"].append([])
                            elif line.find("========") == -1:
                                self.calcRes["frequencies"][-1].append(
                                    float(line.split()[2])
                                )
                                self.calcRes["Sym"][-1].append(int(line.split()[0]))
                                self.calcRes["Mode"][-1].append(int(line.split()[1]))
                            else:
                                break
                            line = streamIn.readline()
                    if line.find("==============") >= 0:
                        break
                    line = streamIn.readline()
                    # read if not ------
            line = streamIn.readline()
        self.calcTask["calculationResults"] = self.calcRes
        self.calculations.append(self.calcTask)

    def readTiming(self, line1, line2):
        vars = line1.split(":")[1]
        vars = vars.split()[0]
        vars2 = line2.split(":")
        vars2 = vars.split()[0]

        return {"cpuTime": float(vars), "wallTime": float(vars2), "units": "second"}

    def readTaskTimes(self, line1, line2):
        vars = line1.split(":")[1]
        vars = vars.split()[0]
        vars2 = line2.split(":")
        vars2 = vars.split()[0]

        return {"cpuTime": float(vars), "wallTime": float(vars2), "units": "second"}
