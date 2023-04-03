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


def get_spectre_info(logs_files):
    """."""
    data = {}
    
    for i, file in enumerate(logs_files):
        data[i] = read_UVVis(file)
        
    print("Number of files analyzed:", len(data))
    
    spectra = []
    w_s1 = []
    w_s2 = []
    for i in data:
        spectra.append(epsilon_tot(w_nm, data[i]))
        w_s1.append(data[i].loc[0, "wl"])
        w_s2.append(data[i].loc[1, "wl"])
        
    spectra = np.array(spectra).mean(axis=0)
    w_s1 = np.array(w_s1)
    w_s2 = np.array(w_s2)
    
    return spectra, w_s1, w_s2


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
