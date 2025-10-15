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

def concat_masks(path, log_path, seq_speed=400, pores=512):
    log_data = np.genfromtxt(log_path, delimiter=',', dtype=None, encoding='utf8')
    log_data = log_data[log_data[:,0] == "boss"]
    dump_filter = (np.int64(np.int64(log_data[:, 4]) / 1000 -1))
    def natural_key(string_):
        """See https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
    filenames = next(walk(path), (None, None, []))[2]
    filenames = sorted(filenames, key=natural_key)[1:]  # sort in correct order and remove boss.npz to avoid duplication
    filenames_filt = np.array(filenames)[dump_filter].tolist()
    fw = [np.load(path+f)["chr21"][:,0] for f in [filenames[0]] + filenames_filt]
    rev = [np.load(path+f)["chr21"][:,1] for f in [filenames[0]] + filenames_filt]
    masks = [np.vstack(fw), np.vstack(rev)]

    time = np.hstack((np.array(0.0), np.float64(log_data[:, 1])/seq_speed/pores/60))

    return masks, time

def plot_masks(masks_str, time, out_path, cutout_lower=None,cutout_higher=None):
    fig, ax = plt.subplots(1,2, width_ratios=[5, 1], layout="constrained", sharey = True)
    fig.set_figwidth(16)
    if cutout_higher and cutout_lower is not None:
        X, Y = np.meshgrid(np.arange(cutout_lower, cutout_higher), time)
        # im = ax[0].imshow(Masks_fw[:, cutout_lower:cutout_higher], aspect='auto', cmap=plt.cm.gray, origin='lower', interpolation='none', extent=[cutout_lower, cutout_higher])  # (0:black, 1:white)
        im = ax[0].pcolormesh(X, Y, masks_str[:, cutout_lower:cutout_higher], cmap=plt.cm.gray)  # (0:black, 1:white)
        ax[0].hlines(time, xmin=cutout_lower, xmax=cutout_higher)
    else:
        X, Y = np.meshgrid(np.arange(masks_str.shape[1]), time)
        # im = ax[0].imshow(Masks_fw, aspect='auto', cmap=plt.cm.gray, origin='lower', interpolation='none')  # (0:black, 1:white)
        im = ax[0].pcolormesh(X, Y, masks_str, cmap=plt.cm.gray)  # (0:black, 1:white)
        ax[0].hlines(time, xmin=0, xmax=masks_str.shape[1])
    # get the colors of the values, according to the 
    # colormap used by imshow
    colors = [ im.cmap(im.norm(value)) for value in masks_str]
    # create a patch (proxy artist) for every color
    values = np.unique(masks_str.ravel())
    patches = [ mpatches.Patch(color=str(np.unique(colors)[0]), label="Reject"),
                mpatches.Patch(facecolor=str(np.unique(colors)[1]), label="Accept", edgecolor='black')]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    ax[0].set_xlabel('genome positions in 100bp windows')
    ax[0].set_ylabel('pseudominutes')
    ax[0].set_ylim(bottom=0)
    # ax[0].grid
    # ax[0].set(yticks=np.round(time), yticklabels=np.round(time));

    if cutout_higher and cutout_lower is not None:
        ax[1].plot(np.average(masks_str[:, cutout_lower:cutout_higher], axis=1), time)
        print(f"{np.average(masks_str[:, cutout_lower:cutout_higher])}")
    else:
        ax[1].plot(np.average(masks_str, axis=1), time)
        print(f"{np.average(masks_str)}") # 0.33399976935325765 for accept unmapped third correction: 0.43869772746426083 for reject unmapped
    ax[1].set_xlabel('prop. of accepted reads')
    ax[1].set_xlim(left=0)
    ax[1].set_ylabel('')

    fig.savefig(out_path, bbox_inches="tight")

    return


if __name__ == "__main__":
    # input to this script is the logfile produced by the run
    masks, time = concat_masks(path=sys.argv[1], log_path=sys.argv[2])
    Masks_fw = masks[0]
    Masks_rev = masks[1]
    plot_masks(Masks_fw, time, sys.argv[3])
    plot_masks(Masks_rev, time, sys.argv[4])