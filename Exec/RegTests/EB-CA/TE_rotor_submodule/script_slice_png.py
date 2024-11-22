import numpy as np
from matplotlib import pyplot as plt
import yt
import glob
import os
import re
from tqdm import tqdm
from multiprocessing import Pool

# VD: Print only warnings and errors
yt.utilities.logger.set_log_level("warning")

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)

def save_png(pltfiles, idx, field, f_min, f_max, angle):
    ds = yt.load(pltfiles[idx])

    # Create a slice of a field along the x axis
    plt_num = idx
    p = yt.SlicePlot(ds, angle, [("boxlib", field)])
    p.set_zlim(field, zmin=f_min, zmax=f_max)
    p.save(f"./figures/fig_{field}_{angle}_{plt_num:04d}.png")

def f_min_max(pltfiles, i, field):
    ds = yt.load(pltfiles[i])
    ad = ds.all_data()
    ray = ds.ortho_ray(0, (0, 0))
    ray_sort = np.argsort(ray[("boxlib", "x")])
    srt = np.array(ray[("boxlib", field)][ray_sort])

    return srt


fdir = "./figures/"
field = "magvel"
angle = "z"

pltfiles = natural_sort(glob.glob(os.path.join(fdir, "plt*")))
n_plots = len(pltfiles)

min_each = np.zeros((n_plots, len(field)))
max_each = np.zeros((n_plots, len(field)))

jobs1 = []
with Pool(processes=100) as pool:
    for i in tqdm(range(n_plots)):
        job1 = pool.apply_async(f_min_max, args=(pltfiles, i, field))
        jobs1.append(job1)

    for job1 in tqdm(jobs1):
        srt = job1.get()

f_max = np.max(srt)
f_min = np.min(srt)

jobs2 = []
with Pool(processes=100) as pool:
    for idx in tqdm(range(n_plots)):
        job2 = pool.apply_async(save_png, args=(pltfiles, idx, field, f_min, f_max, angle))
        jobs2.append(job2)

    for job2 in tqdm(jobs2):
        job2.get()
