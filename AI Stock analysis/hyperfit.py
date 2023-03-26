import numpy as np
import sys
import random

# Constants
nstock = 4
maxyear = 2500
mpar = 16

# Custom data structures
class PopData:
    def __init__(self, year, pop):
        self.year = year
        self.pop = pop

class Range:
    def __init__(self, low, high):
        self.low = low
        self.high = high

class ParData:
    def __init__(self):
        self.var = [Range(0, 0) for _ in range(mpar)]

class Stock:
    def __init__(self, ecosphere, humansphere, ignorance, knowledge):
        self.ecosphere = ecosphere
        self.humansphere = humansphere
        self.ignorance = ignorance
        self.knowledge = knowledge

# Remaining code to read population and parameter data, and process the simulation
def initialize_stocks(stocks, pardata):
    stocks["ecosphere"] = pardata["var"][7]["low"]
    stocks["humansphere"] = pardata["var"][6]["low"]
    stocks["ignorance"] = pardata["var"][5]["low"]
    stocks["knowledge"] = pardata["var"][8]["low"]

def simulate(pardata, stocks, trajectory):
    for i in range(maxyear):
        pE = equation_pE(ecosphere=stocks["ecosphere"], humansphere=stocks["humansphere"])
        domestication = (
            equation_g(
                I_0=pardata["var"][5]["low"],
                t=pardata["var"][9]["low"],
                pE=pE,
                u=pardata["var"][10]["low"],
                w=pardata["var"][12]["low"],
                p=pardata["var"][13]["low"],
                py=pardata["var"][14]["low"],
                sy=pardata["var"][15]["low"],
                year=i
            ) * stocks["humansphere"]
        )
        rewilding = stocks["humansphere"] * stocks["ignorance"]
        learning = stocks["knowledge"] * pardata["var"][4]["low"]
        obsolescence = stocks["knowledge"] * equation_r(v=pardata["var"][11]["low"], pE=pE)
        
        stocks["humansphere"] = stocks["humansphere"] + domestication - rewilding
        stocks["ecosphere"] = stocks["ecosphere"] - domestication + rewilding
        stocks["ignorance"] = stocks["ignorance"] - learning + obsolescence
        stocks["knowledge"] = stocks["knowledge"] + learning - obsolescence
        
        trajectory[i] = equation_population(stocks, pardata)

def rmsd(data_pop, trajectory, startyear, endyear):
    x = 0.0
    for y in range(startyear, endyear+1):
        x += (data_pop[y] - trajectory[y])**2
    x = np.sqrt(x / (endyear - startyear + 1))
    return x

def read_parameters(lines):
    pardata = {"var": [{"low": 0.0, "high": 0.0} for _ in range(16)]}
    
    for line in lines:
        if line.startswith("VARIABLE"):
            var, x, y = line[8:].split()
            x, y = float(x), float(y)
            if y == x:
                y = 0.0
            if y != 0 and x > y:
                x, y = y, x
                
            var_map = {"a": 0, "b": 1, "c": 2, "d": 3, "k": 4, "I_0": 5, "H_0": 6, "E_0": 7,
                       "K_0": 8, "t": 9, "u": 10, "v": 11, "w": 12, "p": 13, "py": 14, "sy": 15}
            idx = var_map.get(var, -1)
            if idx >= 0:
                pardata["var"][idx]["low"] = x
                pardata["var"][idx]["high"] = y
                
    return pardata

def read_pop_data(lines):
    popdata = []
    data_pop = np.zeros(maxyear+1)
    
    for line in lines:
        year, pop = map(float, line.split(','))
        year = int(year)
        if year > maxyear:
            break
        popdata.append({"year": year, "pop": pop})
    
    ndat = len(popdata)
    myear = popdata[-1]["year"]
    year = popdata[0]["year"]
    oldpop = popdata[0]["pop"]
    
    for j in range(ndat):
        for i in range(year+1, popdata[j]["year"]+1):
            data_pop[i] = (oldpop * (popdata[j]["year"] - i) + popdata[j]["pop"] * (i - year)) / (popdata[j]["year"] - year)
        year = popdata[j]["year"]
        oldpop = popdata[j]["pop"]
    
    return data_pop

def equation_g(I_0, t, pE, u, w, p, py, sy, year):
    pp = 1 - p
    g = I_0 + math.log(2) / t
    g *= (1 - np.exp(u * pE))

    if pE < w:
        q = year - sy
        if q <= 0.0:
            x = 1.0
        elif q >= py:
            x = np.exp(-10.0 * (w - pE))
            x = pp * (1 - x) + x
        else:
            x = np.exp(-10.0 * (w - pE))
            q /= py
            x = (q * pp + (1 - q)) * (1 - x) + x
        g *= x

    return g

def equation_growth(pop, data_pop):
    x = data_pop - pop
    x /= pop
    return x

def equation_logData(data_pop):
    return math.log(data_pop + 1) / math.log(10)

def equation_logHumans(pop):
    return math.log(pop + 1) / math.log(10)

def equation_pE(ecosphere, humansphere):
    x = ecosphere / (ecosphere + humansphere)
    return x

def equation_population(stocks, pardata):
    pE = equation_pE(ecosphere=stocks['ecosphere'], humansphere=stocks['humansphere'])
    cce = equation_cce(pE=pE, a=pardata['var'][0]['low'])
    cch = equation_cch(pE=pE, d=pardata['var'][3]['low'], Knowledge=stocks['knowledge'], Io=pardata['var'][5]['low'])
    cc = equation_cc(cch, cce, b=pardata['var'][1]['low'], c=pardata['var'][2]['low'])
    x = cc * stocks['humansphere']
    return x

def equation_r(v, pE):
    x = 0.5 * np.exp(v * pE)
    return x

def equation_cc(cch, cce, b, c):
    return (b + c * cch) * cce

def equation_cch(pE, d, Knowledge, Io):
    pH = 1.0 - pE
    cch = pH * (1 - np.exp(d * Knowledge))
    return cch


def equation_cce(pE, a):
    if 2 * a >= (1 + pE):
        cce = 0.0
    else:
        cce = pE ** (0.5 / (1 + pE - 2 * a))
    return cce

def output_traj(trajectory, data_pop):
    for i in range(1, 2501):
        print(f"{i}, {data_pop[i - 1]}, {trajectory[i - 1]}")

def output_rmsd(rmsd):
    print(f"{rmsd}")
    
# start simulation
if __name__ == "__main__":
    # Initialization
    popdata = [PopData(0, 0.0) for _ in range(200)]
    pardata = ParData()
    newpardata = ParData()
    stocks = Stock(0, 0, 0, 0)
    newstocks = Stock(0, 0, 0, 0)
    data_pop = np.zeros(maxyear)
    trajectory = np.zeros(maxyear)
    varnames = "a   b   c   d   k   I_0 H_0 E_0 K_0 t   u   v   w   p   py  sy  "

    # Read command line arguments
    argc = len(sys.argv)

    if argc < 3:
        print("Usage: python hyperfit.py <popdata> <parameters> [maxgrid] [binaryfile] [startyear] [endyear] [cutoff]")
        sys.exit(1)

    popfile = sys.argv[1]
    parfile = sys.argv[2]
    maxgrid = int(sys.argv[3]) if argc >= 4 else 5
    binaryfile = sys.argv[4] if argc >= 5 else None
    startyear = int(sys.argv[5]) if argc >= 6 else 1970
    endyear = int(sys.argv[6]) if argc >= 7 else 2010
    rmscut = float(sys.argv[7]) if argc >= 8 else None
    sample = True if rmscut is not None else False

    # Read population data
    with open(popfile, "r") as f:
        lines = f.readlines()
        data_pop = read_pop_data(lines)

    # Read parameter data
    with open(parfile, "r") as f:
        lines = f.readlines()
        pardata = read_par_data(lines)
    
    # Initialize stocks
    stocks.humansphere = pardata.var[4].low
    stocks.knowledge = pardata.var[7].low
    stocks.ecosphere = pardata.var[6].low
    stocks.population = pardata.var[8].low

    # Initialize trajectory
    trajectory[0] = equation_logHumans(stocks.population)
    
    # start simulatio
    for i in range(1, maxyear):
        pE = equation_pE(ecosphere=stocks.ecosphere, humansphere=stocks.humansphere)
        r = equation_r(v=pardata.var[11].low, pE=pE)
        g = equation_g(I_0=pardata.var[9].low, t=pardata.var[10].low, pE=pE, u=pardata.var[12].low, w=pardata.var[13].low, p=pardata.var[14].low, py=pardata.var[15].low, sy=pardata.var[16].low, year=i)
        trajectory[i] = g * trajectory[i - 1] + (1 - g) * equation_logHumans(stocks.population)
        stocks.population = equation_population(stocks=stocks, pardata=pardata)
        stocks.ecosphere = stocks.ecosphere * (1 - r)
        stocks.humansphere = stocks.humansphere * (1 - r) + stocks.population
        stocks.knowledge = stocks.knowledge * (1 - r) + r

    # Output trajectory
    output_traj(trajectory, data_pop)