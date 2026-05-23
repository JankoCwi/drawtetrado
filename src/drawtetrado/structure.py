import json
import sys
import math
import subprocess

from drawtetrado.svg_painter import Point, ConnType, ConnFlow


class Nucleotide:
    def FindConnections(self, used):
        connected_to = ""
        curr_min = sys.maxsize
        for name, nucl in used.items():
            if nucl["index"] > self.index and \
               curr_min > nucl["index"] and \
               self.chain == nucl["chain"]:
                curr_min = nucl["index"]
                connected_to = name
        return connected_to

    def Block(self, width, height, angle):
        self.coords = []
        self.coords.append(Point(0.0, 0.0))

        if math.tan(math.radians(angle)) != 0:
            shift = height / math.tan(math.radians(angle))
        else:
            shift = height / math.tan(math.radians(0.01))

        self.coords.append(Point(shift, -height))
        self.coords.append(Point(shift + width, -height))
        self.coords.append(Point(width, 0.0))
        self.center = Point((shift + width) / 2.0, -height / 2.0)

    def ShiftBlock(self, shift):
        for i in range(4):
            self.coords[i] += shift
        self.center += shift

    def CalculateCoordinates(self, conf):
        shift = Point(0, 0)

        sin_val = math.sin(math.radians(conf.angle))
        cos_val = math.cos(math.radians(conf.angle))
        tan_val_rev = math.tan(math.radians(90.0 - conf.angle))

        shift.y -= self.tetrade_no * (
            sin_val * (conf.longer + conf.shorter + conf.spacing) +
            conf.tetrade_spacing
        )

        if self.GetOnzPlusMinus() == "-":
            width_0 = conf.longer
            height_0 = sin_val * conf.shorter
            width_1 = conf.shorter
            height_1 = sin_val * conf.longer
        else:
            width_0 = conf.shorter
            height_0 = sin_val * conf.longer
            width_1 = conf.longer
            height_1 = sin_val * conf.shorter

        if self.position == 0:
            self.Block(width_0, height_0, conf.angle)

        elif self.position == 1:
            self.Block(width_1, height_1, conf.angle)
            shift.y -= (width_1 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) * cos_val

        elif self.position == 2:
            self.Block(width_0, height_0, conf.angle)
            shift.y -= (width_0 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) + \
                        ((width_0 + conf.spacing) * sin_val) * tan_val_rev

        elif self.position == 3:
            self.Block(width_1, height_1, conf.angle)
            shift.x += width_0 + conf.spacing

        self.ShiftBlock(shift)

    def GetOnzPlusMinus(self):
        return "-" if self.onz_full[-1] == "-" else "+"

    def GetOnz(self):
        res = self.onz_full[0:1].lower() + "_"
        res += "minus" if self.GetOnzPlusMinus() == "-" else "plus"

        if res not in (
            "o_plus", "o_minus",
            "n_plus", "n_minus",
            "z_plus", "z_minus"
        ):
            return "onz_default"
        return res

    def __init__(self, data, used_nucl, tetr_no, tetr_onz, pos):
        self.number = data["number"]
        self.short_name = data["shortName"]
        self.full_name = data["fullName"]
        self.chain = data["chain"]
        self.index = data["index"]
        self.onz_full = tetr_onz
        self.onz = self.GetOnz()
        self.bond = data["glycosidicBond"]
        self.tetrade_no = tetr_no
        self.position = pos
        self.connected_to = self.FindConnections(used_nucl)
        self.connected_from = ""

        self.coords = []
        self.center = Point(0, 0)

        self.connection_type = ConnType.UNKNOWN
        self.flow_out = ConnFlow.UNKNOWN
        self.flow_in = ConnFlow.UNKNOWN

        self.priority_conn = -1
        self.priority_edge = -1
        self.priority_nucl = -1


class Quadruplex:
    def UsedNucleotides(self, tetrad, nucl):
        used = {}
        for _, data in tetrad.items():
            used[data["nt1"]] = nucl[data["nt1"]]
            used[data["nt2"]] = nucl[data["nt2"]]
            used[data["nt3"]] = nucl[data["nt3"]]
            used[data["nt4"]] = nucl[data["nt4"]]
        return used

    def PrepareNucleotides(self, structure, quadruplex_id, tetrad_id):
        nucl = structure.nucleotides

        if tetrad_id >= 0:
            tetrads = structure.single_tetrads[quadruplex_id][tetrad_id]
        else:
            tetrads = structure.tetrads[quadruplex_id]

        tetrads_order = structure.tetrads_order[quadruplex_id]
        used_nucl = self.UsedNucleotides(tetrads, nucl)

        self.tetrads = []
        self.nucl_quad = {}

        tetr_no = 0

        for tetrad_name in tetrads_order:
            if tetrad_name not in tetrads:
                continue

            tetrad = tetrads[tetrad_name]

            nt1 = tetrad["nt1"]
            nt2 = tetrad["nt2"]
            nt3 = tetrad["nt3"]
            nt4 = tetrad["nt4"]
            onz = tetrad["onz"]

            self.nucl_quad[nt1] = Nucleotide(nucl[nt1], used_nucl, tetr_no, onz, 0)
            self.nucl_quad[nt2] = Nucleotide(nucl[nt2], used_nucl, tetr_no, onz, 1)
            self.nucl_quad[nt3] = Nucleotide(nucl[nt3], used_nucl, tetr_no, onz, 2)
            self.nucl_quad[nt4] = Nucleotide(nucl[nt4], used_nucl, tetr_no, onz, 3)

            self.tetrads.append([nt1, nt2, nt3, nt4])
            tetr_no += 1

    def GetChainFirstLast(self):
        chains = {}

        for _, nucl in self.nucl_quad.items():
            if nucl.chain not in chains:
                chains[nucl.chain] = {"first": "", "last": "", "val": sys.maxsize}

        for name, nucl in self.nucl_quad.items():
            if chains[nucl.chain]["val"] > nucl.index:
                chains[nucl.chain]["val"] = nucl.index
                chains[nucl.chain]["first"] = name

        for chain, data in chains.items():
            curr = data["first"]
            nxt = self.nucl_quad[curr].connected_to

            while nxt != "":
                self.nucl_quad[nxt].connected_from = curr
                curr = nxt
                nxt = self.nucl_quad[nxt].connected_to

            data["last"] = curr

        return chains

    def __init__(self, structure, quadruplex_id, tetrad_id=-1):
        self.nucl_quad = {}
        self.tetrads = []

        self.PrepareNucleotides(structure, quadruplex_id, tetrad_id)
        self.chains = self.GetChainFirstLast()

        if tetrad_id >= 0:
            self.tracts = [structure.tracts[quadruplex_id][tetrad_id]]
        else:
            self.tracts = structure.tracts[quadruplex_id]

    def GetNucleotidesPositions(self):
        lst = [-1] * len(self.nucl_quad)

        for name, nucl in self.nucl_quad.items():
            if nucl.connected_to != "":
                conn = self.nucl_quad[nucl.connected_to]
                lst[nucl.tetrade_no * 4 + nucl.position] = (
                    conn.tetrade_no * 4 + conn.position
                )
            else:
                lst[nucl.tetrade_no * 4 + nucl.position] = -1

        return lst

    def GetSameRotations(self):
        lst = [-1] * len(self.tetrads)
        return lst

    def GetAlignments(self):
        lst = [-1] * len(self.nucl_quad)
        return lst

    def Optimize(self):
        import svg_optimizer as optimizer_mod

        optimized = optimizer_mod.solve(
            self.GetNucleotidesPositions(),
            self.GetSameRotations(),
            self.GetAlignments()
        )

        optimized = list(map(int, optimized))

        for _, nucl in self.nucl_quad.items():
            base = nucl.tetrade_no * 4

            for x in range(4):
                if optimized[base + x] == nucl.position:
                    nucl.position = x
                    break

        for i, tetrad in enumerate(self.tetrads):
            base = i * 4
            new_tetrad = [None] * 4

            for j in range(4):
                new_tetrad[j] = tetrad[optimized[base + j]]

            self.tetrads[i] = new_tetrad


class Structure:
    def __init__(self):
        self.nucleotides = {}
        self.tetrads = []
        self.tetrads_order = []
        self.single_tetrads = []
        self.tracts = []

    def addNucleotide(self, name, data):
        self.nucleotides[name] = data

    def fromFile(self, path):
        with open(path) as file:
            return self.fromJsonDict(json.load(file))

    def fromString(self, json_string):
        return self.fromJsonDict(json.loads(json_string))

    def fromJsonDict(self, json_dict):
        for data in json_dict["nucleotides"]:
            self.addNucleotide(data["fullName"], data)

        for helice in json_dict["helices"]:
            single_tetrads_local = []
            tetrad_unordered = {}
            tracts_all = []

            for quad in helice["quadruplexes"]:
                tetrad_unordered_local = {}

                for data in quad["tetrads"]:
                    tetrad_unordered[data["id"]] = data
                    tetrad_unordered_local[data["id"]] = data

                single_tetrads_local.append(tetrad_unordered_local)

                tracts_all.append(quad.get("tracts", []))

            tetrad_ordered = []
            tetrad_pairs = helice.get("tetradPairs", [])[:]

            if tetrad_pairs:
                pair = tetrad_pairs.pop(0)
                tetrad_ordered.extend([pair["tetrad1"], pair["tetrad2"]])

            while tetrad_pairs:
                found = False
                for i, pair in enumerate(tetrad_pairs):
                    for j, t in enumerate(tetrad_ordered):
                        if t == pair["tetrad1"]:
                            tetrad_ordered.insert(j + 1, pair["tetrad2"])
                            found = True
                            break
                        if t == pair["tetrad2"]:
                            tetrad_ordered.insert(j, pair["tetrad1"])
                            found = True
                            break
                    if found:
                        tetrad_pairs.pop(i)
                        break

            tetrad_ordered.reverse()

            if len(tetrad_ordered) > 1:
                self.tetrads_order.append(tetrad_ordered)
                self.tetrads.append(tetrad_unordered)

            self.single_tetrads.append(single_tetrads_local)
            self.tracts.append(tracts_all)

        return self
