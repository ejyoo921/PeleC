import numpy as np
from matplotlib import pyplot as plt
import yt
import argparse
import glob
import os
import re
#import cv2

# VD: Print only warnings and errors
yt.utilities.logger.set_log_level("warning")

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)

fdir = "./"
field = "magvel"

pltfiles = natural_sort(glob.glob(os.path.join(fdir, "plt*")))
n_plots = len(pltfiles)

min_each = np.zeros((n_plots, len(field)))
max_each = np.zeros((n_plots, len(field)))

for i in range(n_plots):
    ds = yt.load(pltfiles[i])
    ad = ds.all_data()
    ray = ds.ortho_ray(0, (0, 0))
    ray_sort = np.argsort(ray[("boxlib", "x")])
    srt = np.array(ray[("boxlib",   field)][ray_sort])

f_max = np.max(srt)
f_min = np.min(srt)

for idx in range(n_plots):
    ds = yt.load(pltfiles[idx])

    # Create a slice of a field along the x axis
    plt_num = idx
    p = yt.SlicePlot(ds, "x", [("boxlib", field)])

    p.set_zlim(field, zmin=f_min, zmax=f_max)
    

    #p.save("./figures/fig_"+field+str(plt_num)+".png")
    p.save(f"./figures/fig_{field}_{plt_num:04d}.png")

#image_folder = './figures/'
#video_name = './movies/mv_'+field+'.avi'

#images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
#frame = cv2.imread(os.path.join(image_folder, images[0]))
#height, width, layers = frame.shape

#video = cv2.VideoWriter(video_name, 0, 1, (width,height))

#for image in images:
#    video.write(cv2.imread(os.path.join(image_folder, image)))

#cv2.destroyAllWindows()
#video.release()
