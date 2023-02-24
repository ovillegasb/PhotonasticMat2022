"""This is a file that loads routine configurations used in Notebooks."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

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
        "dmax": ["mean", "std"]
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
