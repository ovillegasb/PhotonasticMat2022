"""This is a file that loads routine configurations used in Notebooks."""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.constants import e, h, c, k
from scipy.stats import norm

""" Matplotlib tools. """

# Matplotlib parameters to change.
MplotParam = {
    "figure.dpi": 100,
    "figure.figsize": [9, 6]
}

# Box on text in figures.
boxText = {
        "facecolor": "0.85",
        "edgecolor": "k",
        "boxstyle": "round"
    }

# aggregation for pandas
aggregation = {
        "Rg": ["mean", "std"],
        "k2": ["mean", "std"],
        "dmax": ["mean", "std"],
        "bin": "count"
    }

""" Some Definitions. """

# Vector with wavelength values (in nm)
w_nm = np.linspace(100, 800, num=1401)

# Hatree to kJ/mol
Ha_to_kJmol = 2625.5

# Hatree to J
Ha_to_J = 4.359748199e-18

# PROP figures pol
POLPROP = {
    "Rg": {
        "color": "#2a90a6",
        "binwidth": .1,
        "xlim": [0.5, 3.0],
        "label": "Radius of gyration (nm)",
        "units": "nm"
    },
    "dmax": {
        "color": "#5b557b",
        "binwidth": .3,
        "xlim": [1.0, 9.0],
        "label": "Max. distance (nm)",
        "units": "nm"
    },
    "k2": {
        "color": "#fea6ad",
        "binwidth": .05,
        "xlim": [0, 1],
        "label": "Shape anisotropy",
        "units": ""
    },
    "rdf": {
        "color": "#3293f0",
        "xlim": [0.25, 2.5],
        "label": "g(r) PC - POL",
        "units": ""

    },
    "distance": {
        "color": "#66c0c0",
        "binwidth": .1,
        "xlim": [0, 4.0],
        "label": "distance",
        "units": "nm"
    }
}


""" Some Functions. """


def rot(x):
    """Rotate the negative angles by 360 degrees."""
    if x < 0:
        return x + 360
    else:
        return x


def toTime(frame, freq=1, t0=0.0):
    """Convert regular frames to time (in ps)."""
    return frame * freq + t0 


def binned_distance(x, rmin, rmax, binwidth):
    """Binned distances."""
    bins = np.arange(rmin, rmax, binwidth)
    index = np.digitize(x, bins)
    try:
        return bins[index]
    except IndexError:
        return np.NaN


def save_fig(name, outfig):
    """Save the figure."""
    plt.savefig(
        f"{outfig}/{name}",
        dpi=300,
        bbox_inches="tight",
        format="png",
        facecolor='white',
        transparent=False
    )


def format_time(seconds):
    """Convert seconds to formated time."""
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{seconds:02d}"


class dihedral_trans:
    """Object representing trans dihedrals."""
    
    def __init__(self, data):
        """Load a vector of angles in degrees."""
        self.data = data
        self.continuity_study()

    def continuity_study(self):
        """Study the data curve."""
        dat = self.data
        mean_val = 0.0
        mean_evol = []
        std_evol = []

        for i, val in enumerate(dat):
            if i == 0:
                mean_val += val
                std_evol.append(0)
            else:
                mean_val = (mean_val + val) / 2
                std_evol.append(np.std(dat[0:i]))

            mean_evol.append(mean_val)

        self.mean_evol = mean_evol
        self.std_evol = std_evol

        # Calculates the absolute value of the gradient
        self.grad_dat = abs(np.gradient(std_evol))

        # Look for peaks that are above a gradient of 1
        peaks, _ = find_peaks(self.grad_dat, height=1.0)
        self.peaks = peaks

        # The first peak found after 15% of the data is taken
        self.start = int(15 * len(self.grad_dat) / 100)
        self.taked_peak = peaks[peaks > self.start][0]
        self.chage_position = np.where(
            self.grad_dat == self.grad_dat[self.taked_peak].max()
            )[0][0]

    def show_evol(self):
        """Plot the evolution of the mean value and the standard deviation."""
        fig, ax = plt.subplots(nrows=2)

        ax[0].set_title("Evolution of the average value")
        ax[0].plot(self.mean_evol)

        ax[1].set_title("Evolution of standard deviation")
        ax[1].plot(self.std_evol)

        plt.tight_layout()
        plt.show()

    def show_grad(self):
        """Plot gradiend."""
        fig, ax = plt.subplots()

        ax.set_title("Gradient")
        ax.plot(self.grad_dat, "green")
        ax.plot(self.peaks, self.grad_dat[self.peaks], "X", color="orange")

        plt.axhline(y=1.0, ls="--", color="gray")
        plt.axvline(x=self.start, ls="--", color="gray")

        plt.show()


def read_UVVis(file):
    """Read the UV-Vis spectrum data from a Gaussian output."""
    ExcitedE = re.compile(r"""
        ^\sExcited\sState\s+(?P<Nstate>\d+):\s+\w+-\w\s+(?P<E>\d+.\d+)\seV\s+(?P<wl>\d+.\d+)\snm\s+f=(?P<f>\d+.\d+)
    """, re.X)

    Transitions = re.compile(r"""
        ^\s+(?P<states>\d+\s->\s\d+)\s+(?P<p>[+-]?\d+\.\d+)
    """, re.X)
    
    state = 0
    UVVis = {}
    with open(file, "r") as f:
        UVVis["states"] = {}
        UVVis["transition"] = {}
        for line in f:
            if ExcitedE.match(line):
                m = ExcitedE.match(line)
                state = int(m.groupdict()["Nstate"])
                # print(m.groupdict())
                UVVis["states"][state] = m.groupdict()
                UVVis["transition"][state] = []
        
            elif Transitions.match(line):
                m = Transitions.match(line)
                # print(m.groupdict())
                UVVis["transition"][state].append(m.groupdict())
                
    transition = UVVis["transition"]
    for n in transition:
        tab = pd.DataFrame(transition[n])
        tab = tab.astype({"p": np.float64})
        tab["test"] = tab["p"]**2
    
        ou = tab[tab["test"] == tab["test"].max()].values[0][0:-1]
        UVVis["states"][n]["t"] = ou[0]
        UVVis["states"][n]["p"] = ou[1]
        
    table = []
    for i in UVVis["states"]:
        table.append(UVVis["states"][i])
        
    table = pd.DataFrame(table)
    table = table.astype({
        "E": np.float64,
        "wl": np.float64,
        "f": np.float64,
        "Nstate": np.int64
    })

    return table


def epsilon_i(wl, wl_i, f, sigma=0.4):
    """
    Return the epsilon value for a frequency.

    Parameters:
    -----------
    wl : array

    sigma : float
        eV
    """    
    # Convert to cm-1
    sigma_cm = (sigma * 1e-2 * e / h / c)
    sigma_nm = (sigma * 1e-9 * e / h / c)
    
    return 1.3062974e8 * (f / sigma_cm) * np.exp(-((((1/wl) - (1/wl_i))/(sigma_nm))**2))


def epsilon_tot(wl, table, sigma=0.4, normalize=False):
    """Return the epsilon value for all frequency."""
    try:
        eps = np.zeros(len(wl))
    except TypeError:
        eps = 0.0

    for i in table.index:
        eps += epsilon_i(wl, table.loc[i, "wl"], table.loc[i, "f"])
        
    if normalize:
        eps = eps / eps.max()
        
    return eps


def dnorm(values):
    """."""
    return norm.pdf(values, values.mean(), values.std())


def get_spectre_info(logs_files, dtype="cont", sigma=0.4):
    """."""
    data = {}
    info = []
    
    for i, file in enumerate(logs_files):
        data[i] = read_UVVis(file)
        metadata = file.split("/")[4:]
        # print(i, metadata[0], metadata[1], metadata[-1].split(".")[-2].split("_")[-1])
        if dtype=="GMX":
            info.append({
                "index": i,
                "isomer": metadata[0].split(".")[0].split("_")[1],
                "replica": int(metadata[1].replace("md", "")),
                "frame": int(metadata[-1].split(".")[-2].split("_")[-1])
            })
        else:
            try:
                info.append({
                    "index": i,
                    "isomer": metadata[0],
                    "replica": int(metadata[1].split("_")[-1]),
                    "frame": int(metadata[-1].split(".")[-2].split("_")[-1])
                })
            except ValueError:
                metadata = [
                    file.split("/")[2].split(".")[0].split("_")[1],
                    0,
                    file.split("/")[-1].split(".")[0].split("_")[-1]
                ]
                info.append({
                    "index": i,
                    "isomer": metadata[0],
                    "replica": metadata[1],
                    "frame": int(metadata[2])
                })
        # print(file.split("/")[-1].split(".")[-2].split("_")[-1])
        
    print("Number of files analyzed:", len(data))
    
    spectra = []
    w_s1 = []
    w_s2 = []
    w_s3 = []
    w_s4 = []
    w_s5 = []
    w_s6 = []

    f_s1 = []
    f_s2 = []
    f_s3 = []
    f_s4 = []
    f_s5 = []
    f_s6 = []

    for i in data:
        spectra.append(epsilon_tot(w_nm, data[i], sigma))
        w_s1.append(data[i].loc[0, "wl"])
        w_s2.append(data[i].loc[1, "wl"])
        w_s3.append(data[i].loc[2, "wl"])
        w_s4.append(data[i].loc[3, "wl"])
        w_s5.append(data[i].loc[4, "wl"])
        w_s6.append(data[i].loc[5, "wl"])

        f_s1.append(data[i].loc[0, "f"])
        f_s2.append(data[i].loc[1, "f"])
        f_s3.append(data[i].loc[2, "f"])
        f_s4.append(data[i].loc[3, "f"])
        f_s5.append(data[i].loc[4, "f"])
        f_s6.append(data[i].loc[5, "f"])
        
    spectra = np.array(spectra).mean(axis=0)
    w_s1 = np.array(w_s1)
    w_s2 = np.array(w_s2)
    w_s3 = np.array(w_s3)
    w_s4 = np.array(w_s4)
    w_s5 = np.array(w_s5)
    w_s6 = np.array(w_s6)

    f_s1 = np.array(f_s1)
    f_s2 = np.array(f_s2)
    f_s3 = np.array(f_s3)
    f_s4 = np.array(f_s4)
    f_s5 = np.array(f_s5)
    f_s6 = np.array(f_s6)

    dfinfo = pd.DataFrame(info)
    dfinfo.set_index("index", inplace=True)
    dfinfo["s1"] = w_s1
    dfinfo["s2"] = w_s2
    dfinfo["s3"] = w_s3
    dfinfo["s4"] = w_s4
    dfinfo["s5"] = w_s5
    dfinfo["s6"] = w_s6

    dfinfo["f1"] = f_s1
    dfinfo["f2"] = f_s2
    dfinfo["f3"] = f_s3
    dfinfo["f4"] = f_s4
    dfinfo["f5"] = f_s5
    dfinfo["f6"] = f_s6
    
    return spectra, dfinfo


def showInfo():
    """Print information about definitions in this submodule."""
    print("Functions or objects in nbconfig")
    print("--------------------------------")
    print("""
    MplotParam : (dict)
        Matplotlib parameters to change.

    boxText : (dict)
        Box on text in figures.

    rot : (function)
        Rotate the negative angles by 360 degrees.

    toTime : (function)
        Convert regular frames to time (in ps).

    binned_distance : (function)
        Binned distances (in nm).

    """)

