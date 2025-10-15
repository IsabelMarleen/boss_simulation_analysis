import numpy as np
from os import walk
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re

# path = "/hps/nobackup/goldman/ipoetzsch/sim_2mil_out_filt_2/out_boss_runs_sim_2mil/masks/boss.npz"
# mask = np.load(path)

# mask["chr21"].shape #(450906, 2) instead of expected (45090682, 2) 
# for each position in genome for forward and reverse because of default window size 100

def concat_covs(path):
    def natural_key(string_):
        """See https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
    filenames = next(walk(path), (None, None, []))[2]
    filenames = sorted(filenames, key=natural_key)[1:]  # sort in correct order and remove boss.npz to avoid duplication
    covs = [np.load(path+f) for f in filenames]
    cov = [np.vstack(cov)]

    np.savetxt(cov)

    return cov

def plot_masks(masks, out_path):
    Masks_fw = masks[0]#.astype(bool)
    Masks_rev = masks[1]#.astype(bool)

    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    # ax.eventplot(masks[:15, :3], orientation="vertical", linewidth=0.75)
    # ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
    #    ylim=(0, 8), yticks=np.arange(1, 8))
    im = ax.imshow(Masks_fw, aspect='auto', cmap=plt.cm.gray, origin='lower')  # (0:black, 1:white)

    # get the colors of the values, according to the 
    # colormap used by imshow
    colors = [ im.cmap(im.norm(value)) for value in Masks_fw]
    # create a patch (proxy artist) for every color
    values = np.unique(Masks_fw.ravel())
    patches = [ mpatches.Patch(color=str(np.unique(colors)[0]), label="Reject"),
                mpatches.Patch(facecolor=str(np.unique(colors)[1]), label="Accept", edgecolor='black')]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    plt.xlabel('genome positions in 100bp windows')
    plt.ylabel('batches of 1000 fragments')
    fig.savefig(out_path, bbox_inches="tight")

    return



if __name__ == "__main__":
    # input to this script is the logfile produced by the run
    cov = concat_covs(path=sys.argv[1])