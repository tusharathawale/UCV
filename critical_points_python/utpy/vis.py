""" This module will hold a collection of "show" functions that can be either
be shown onscreen or saved to a file. Some do not give this option as they are
used mostly as a precursor for another function call. 
"""

# Use this backend if you are saving directly to file, otherwise you may need
# to enable an interactive backend, see here:
# https://matplotlib.org/faq/usage_faq.html#what-is-a-backend

# import matplotlib
# matplotlib.use('Agg')

from matplotlib import colors, cm
import matplotlib.pyplot as plt
from skimage import filters, morphology
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from itertools import cycle

import nglpy as ngl
import topopy

from utpy.utils import *

# I played with several different color schemes, uncomment the others to try
# them out instead:
color_list = [
    [228, 26, 28],
    [55, 126, 184],
    [77, 175, 74],
    [152, 78, 163],
    [255, 127, 0],
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
    [153, 153, 153],
]

# color_list = [
#     [141, 211, 199],
#     [255, 255, 179],
#     [190, 186, 218],
#     [251, 128, 114],
#     [128, 177, 211],
#     [253, 180, 98],
#     [179, 222, 105],
#     [252, 205, 229],
#     [217, 217, 217],
# ]
# color_list = [
#     [105, 239, 123],
#     [149, 56, 144],
#     [192, 222, 164],
#     [14, 80, 62],
#     [153, 222, 249],
#     [24, 81, 155],
#     [218, 185, 255],
#     [66, 30, 200],
#     [183, 211, 33],
# ]
# color_list = [
#     [251, 180, 174],
#     [179, 205, 227],
#     [204, 235, 197],
#     [222, 203, 228],
#     [254, 217, 166],
#     [255, 255, 204],
#     [229, 216, 189],
#     [253, 218, 236],
#     [242, 242, 242],
# ]

color_list = np.array(np.array(plt.cm.tab20.colors) * 255, dtype=int)
ccycle = cycle(color_list)


def overlay_alpha_image_lazy(background_rgb, overlay_rgba, alpha):
    """ An alpha blending technique scraped from the internet that performs a 
        "lazy" blending approach. This is from here:
        https://gist.github.com/pthom/5155d319a7957a38aeb2ac9e54cc0999

        cf https://en.wikipedia.org/wiki/Alpha_compositing#Alpha_blending
        If the destination background is opaque, then
        out_rgb = overlay_rgb * overlay_alpha + background_rgb * (1 - overlay_alpha)

    Args:
        background_rgb (numpy.ndarray): The first image as a 3D array where the
        first two dimensions are spatial and the third is a three channel RGB.
        overlay_rgba (numpy.ndarray): The second image blended on top of the
        first as a 3D array where the first two dimensions are spatial and the
        third is a three channel RGB.
        alpha (float): The amount of blending.

    Returns:
        out_rgb (numpy.ndarray): The blended image as a 3D array where the first
        two dimensions are spatial and the third is a three channel RGB.

    """
    overlay_alpha = overlay_rgba[:, :, 3].astype(np.float) / 255.0 * alpha
    overlay_alpha_3 = np.dstack((overlay_alpha, overlay_alpha, overlay_alpha))
    overlay_rgb = overlay_rgba[:, :, :3].astype(np.float)
    background_rgb_f = background_rgb.astype(np.float)
    out_rgb = overlay_rgb * overlay_alpha_3 + background_rgb_f * (1.0 - overlay_alpha_3)
    out_rgb = out_rgb.astype(np.uint8)
    return out_rgb


def overlay_alpha_image_precise(background_rgb, overlay_rgba, alpha, gamma_factor=2.2):
    """ An alpha blending technique scraped from the internet that performs a 
        "precise" blending approach. This is from here:
        https://gist.github.com/pthom/5155d319a7957a38aeb2ac9e54cc0999

        cf https://en.wikipedia.org/wiki/Alpha_compositing#Alpha_blending
        If the destination background is opaque, then
        out_rgb = overlay_rgb * overlay_alpha + background_rgb * (1 - overlay_alpha)

        cf minute physics brilliant clip "Computer color is broken" :
        https://www.youtube.com/watch?v=LKnqECcg6Gw
        the RGB values are gamma-corrected by the sensor (in order to keep
        accuracy for lower luminancy), we need to undo this before averaging.

    Args:
        background_rgb (numpy.ndarray): The first image as a 3D array where the
        first two dimensions are spatial and the third is a three channel RGB.
        overlay_rgba (numpy.ndarray): The second image blended on top of the
        first as a 3D array where the first two dimensions are spatial and the
        third is a three channel RGB.
        alpha (float): The amount of blending.
        gamma_factor (float): The gamma correction factor parameter.

    Returns:
        out_rgb (numpy.ndarray): The blended image as a 3D array where the first
        two dimensions are spatial and the third is a three channel RGB.

    """
    overlay_alpha = overlay_rgba[:, :, 3].astype(np.float) / 255.0 * alpha
    overlay_alpha_3 = np.dstack((overlay_alpha, overlay_alpha, overlay_alpha))

    overlay_rgb_squared = np.float_power(
        overlay_rgba[:, :, :3].astype(np.float), gamma_factor
    )
    background_rgb_squared = np.float_power(
        background_rgb.astype(np.float), gamma_factor
    )
    out_rgb_squared = overlay_rgb_squared * overlay_alpha_3 + background_rgb_squared * (
        1.0 - overlay_alpha_3
    )
    out_rgb = np.float_power(out_rgb_squared, 1.0 / gamma_factor)
    out_rgb = out_rgb.astype(np.uint8)
    return out_rgb


def show_image(image):
    """ Draw a two dimensional scalar field to a pseudocolor image canvas.

    Args:
        image (numpy.ndarray): A 2D scalar field grid of function values.

    Returns:
        None

    """
    plt.figure()
    img = plt.imshow(
        image,
        cmap=plt.cm.viridis,
        norm=colors.LogNorm(vmin=image.min(), vmax=image.max()),
    )
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.colorbar(img, orientation="horizontal")


def show_colormapped_image(image, my_dir, screen=False, filename="height.png"):
    """ Visualize heightmap of a 2D image

    Args:
        image (numpy.ndarray): A 2D scalar field grid of function values.
        my_dir (str): The directory where the file will be saved.
        screen (bool): Flag specifying whether the image should be shown
            onscreen.
        filename (str): Default filename
    Returns:
        None

    """
    plt.figure()
    img = plt.imshow(image, cmap="cividis", vmin=image.min(), vmax=image.max())
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.colorbar(img, orientation="horizontal")
    plt.gca().set_ylim(0, image.shape[0])
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_msc(grid,
             my_dir,
             persistence=None,
             n_clusters=None,
             screen=False,
             filename="msc.png"):
    """
        Visualize Morse complex of a 2D scalar field
    
    Args:
        grid (numpy.ndarray): A 2D scalar field grid of function values.
        my_dir (str): The directory where the file will be saved.
        persistence (float): the amount of simplfication to perform before
            visualizing Morse complex.
        n_clusters (int): the number of partitions to simplify to if the
            persistence is omitted.
        screen (bool): Flag specifying whether the image should be shown
            onscreen.
        filename (str): Flag specifying whether the image should be shown
            onscreen.
    Returns:
        None

    """
    # grid is a 2D array of hXw dimension representing discrete scalar field
    # X is a 2D array of pairs [(x1,y1), (x2,y2), ..]
    # Y is 1D array of function values for indices in array in X
    X, Y = massage_data(grid)
    h, w = grid.shape

    # Build morse complex for discrete field
    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    tmc.build(X, Y)

    # Find persistence value (if no input persistence value is provided), such that with simplification using that value causes number of
    # partitions equal to input number of partitions (n_clusters)
    if persistence is None:
        for p in tmc.persistences:
            if len(tmc.get_partitions(p).keys()) <= n_clusters:
                persistence = p
                break

    # Morse complex partitions of a discrete field
    partitions = tmc.get_partitions(persistence)
    # Keys contains keys used to represent Morse complex partitions
    keys = partitions.keys()

    # remap key values to serial numbers [1,2,..,n_clusters]
    keyMap = {}
    for i, k in enumerate(keys):
        keyMap[k] = i

    #colorList = [
    #                [31, 119, 180], 1f77b4
    #                [148, 103, 189], 9467bd
    #                [214, 39, 40],d62728
    #                [255, 187, 120],ffbb78
    #                [174, 199, 232],aec7e8
    #                [152, 223, 138],98df8a
    #                [255, 152, 150],ff9896
    #                [44, 160, 44],2ca02c
    #                [255, 127, 14]ff7f0e
    #                 ]

    # color list corresponding to Ackley function in the paper
    #colorList = [
    #    "#d62728",
    #    "#aec7e8",
    #    "#9467bd",
    #    "#ff7f0e",
    #    "#98df8a",
    #    "#2ca02c",
    #    "#ff9896",
    #    "#ffbb78",
    #    "#1f77b4"
    #]

    colorList = [
        "#1f78b4",
        "#33a02c",
        "#e31a1c",
        "#ff7f00",
        "#6a3d9a",
        "#b15928",
        "#a6cee3",
        "#b2df8a",
        "#fb9a99",
        "#fdbf6f",
        "#cab2d6",
        "#ffff99",
        "#cccccc",
    ]

# colorList = [
#        "#1f78b4",
#        "#ff7f00",
#        "#e31a1c",
#        "#33a02c",
#        "#6a3d9a",
#        "#b15928",
#        "#a6cee3",
#        "#b2df8a",
#        "#fb9a99",
#        "#fdbf6f",
#        "#cab2d6",
#        "#ffff99",
#        "#cccccc",
#    ]
#    colorList = np.array(np.array(plt.cm.tab20.colors))

    # Generate cycle of colors
    ccycle = cycle(colorList)

    # usedColors contains list of colors to be used to color Morse complex cells
    uniqueCount = len(keys)
    usedColors = []
    for i, c in zip(range(uniqueCount), ccycle):
        usedColors.append(c)
    # Generate colormap
    cmap = colors.ListedColormap(usedColors)

    bounds = np.array([keyMap[k] for k in keys]) - 0.5
    bounds = bounds.tolist()
    bounds.append(bounds[-1] + 1)
    plt.figure()

    # color mesh to store colors of Morse complex info for visualization
    color_mesh = np.zeros((h, w))
    # color mesh contain remapped keys for Morse complex cells unlike 2D mesh partitions
    for key, indices in partitions.items():
        for idx in indices:
            color_mesh[idx // w, idx % w] = keyMap[key]

    # Do color mapping based remapped keys [1,2,..,num_clusters] and color map
    img = plt.imshow(color_mesh, cmap=cmap, interpolation="nearest", origin="lower")
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.colorbar(
        img,
        cmap=cmap,
        ticks=[range(uniqueCount)],
        boundaries=bounds,
        orientation="horizontal",
    )
    # plt.contour(grid, cmap=cm.viridis)
    # Save figure
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    # If screen is set to true, display result upon running script
    if screen:
        plt.show()
    plt.close()


def show_msc_boundaries(grid,
                        persistence=None,
                        n_clusters=None,
                        color="#000000"):
    """ TODO
    
    Args:
        grid (numpy.ndarray): A 2D scalar field grid of function values.
        persistence (float): the amount of simplfication to perform before
            plotting.
        n_clusters (int): the number of partitions to simplify to if the
            persistence is omitted.
        color (str): hex representation of a color.
    Returns:
        None

    """
    X, Y = utpy.utils.massage_data(grid)
    h, w = grid.shape

    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    tmc.build(X, Y)

    if persistence is None:
        for p in tmc.persistences:
            if len(tmc.get_partitions(p).keys()) <= n_clusters:
                persistence = p
                break

    partitions = tmc.get_partitions(persistence)
    keys = partitions.keys()

    keyMap = {}
    levels = []
    for i, k in enumerate(keys):
        keyMap[k] = i
        levels.append(i + 0.5)

    color_mesh = np.zeros((h, w))
    for key, indices in partitions.items():
        for idx in indices:
            color_mesh[idx // w, idx % w] = keyMap[key]

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    for i in keyMap.values():
        plt.contour(color_mesh == i, colors=color, levels=levels, linewidths=1)
    plt.gca().set_aspect("equal")


def show_method_comparison(ground_truth,
                           ensemble,
                           n_clusters,
                           assignments):
    """
    """
    plt.figure()
    mean_realization = np.mean(ensemble, axis=2)
    # median_realization = np.median(ensemble, axis=2)
    show_combined_overlay(ensemble, assignments, 0.2, False)
    show_msc_boundaries(ground_truth, n_clusters=n_clusters, color="#4daf4a")
    show_msc_boundaries(mean_realization, n_clusters=n_clusters, color="#e41a1c")
    plt.plot([-1, -0.5], [0, 1], color="#4daf4a", linewidth=1, label="ground truth")
    plt.plot([-1, -0.5], [0, 1], color="#e41a1c", linewidth=1, label="mean")
    plt.gca().set_xlim(0, 39)
    plt.gca().set_ylim(0, 39)
    _ = plt.legend(
        bbox_to_anchor=(0, 1.02, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=3,
    )


def show_persistence_charts(ensemble,
                            my_dir,
                            screen=False,
                            filename="composite_persistence_charts.png"):
    """
        Analyze effect of persistence simplification on count of Morse complex cells. Specifically plot
        number of Morse complex cells against level of persistence simplfication for each instance of Ensemble,
        and visualize them all together.
        
        Args:
            ensemble: An ensemble of 2D scalar fields
            my_dir: destination to store visualization of persistence chart
            screen: flag to indicate if the result should be displayed on screen during
            execution of function or not
            filename: name of output file
            
        Returns:
            Plot of number of Morse complex cells against level of persistence simplification overlaid on each other for all instances.z
        
    """
    print('here')
    plt.figure()

    all_ps = []
    all_counts = []
    for i in range(ensemble.shape[2]):
        graph = ngl.EmptyRegionGraph(**graph_params)
        tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    
        X, Y = massage_data(ensemble[:, :, i])
        tmc.build(X, Y)
        ps = [0]

        count = len(np.unique(list(tmc.get_partitions(0).keys())))
        counts = [count]
        #eps = 1e-6
        eps = 0
        for i, p in enumerate(tmc.persistences):
            ps.append(p)
            counts.append(count)
            count = len(np.unique(list(tmc.get_partitions(p + eps).keys())))
            ps.append(p)
            counts.append(count)
        all_ps.append(ps)
        all_counts.append(counts)
        plt.plot(ps, counts, alpha=0.2, c="#1f78b4")

    ax = plt.gca()
    ax.set_ylabel('# 2-Cells', fontsize=20)
    ax.set_xlabel('Persistence Simplification Level', fontsize=20)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    # ax.set_ylim(0, 25)
    # plt.axhline(2, 0.10, 0.20, linestyle='dashed', color='#000000')
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight", dpi =200)
    if screen:
        plt.show()
    plt.close()
    return all_ps, all_counts


def show_persistence_chart(realization,
                           my_dir,
                           screen=False,
                           filename="persistence_chart.png"):
    """
    """
    plt.figure()

    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)

    X, Y = massage_data(realization)
    tmc.build(X, Y)
    ps = [0]

    count = len(np.unique(list(tmc.get_partitions(0).keys())))
    counts = [count]
    eps = 1e-6
    for i, p in enumerate(tmc.persistences):
        ps.append(p)
        counts.append(count)
        count = len(np.unique(list(tmc.get_partitions(p + eps).keys())))
        ps.append(p)
        counts.append(count)

    plt.plot(ps, counts, c="#1f78b4")

    ax = plt.gca()
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return ps, counts


def show_survival_count(ensemble,
                        my_dir,
                        screen=False,
                        filename="survival_count.png"):
    """
    """
    all_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        counts, _, _, _ = count_persistence(ensemble[:, :, i])
        all_counts += counts

    plt.figure()

    img = plt.imshow(
        (all_counts - np.min(all_counts)) / (np.max(all_counts) - np.min(all_counts)),
        cmap=cm.viridis,
    )
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    k = (np.max(all_counts) + np.min(all_counts)) / 2.0
    # print("Contour visualized for survival count: {}".format(k))
    # plt.contour(all_counts, levels=[k], colors='#FFFF00')
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_counts


my_cmap = cm.viridis
#my_cmap = colors.ListedColormap(cm.tab20.colors[:16])

def tushar_show_survival_count(ensemble,
                               my_dir,
                               screen=False,
                               filename="tushar_survival_count.png"):
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        # count_persistence function computes survival for single instance and is
        # defined in utils.py
        weighted_counts = np.squeeze(tushar_count_persistence(ensemble[:, :, i]))
        all_weighted_counts += weighted_counts
    
    plt.figure()
    # img = plt.imshow((all_weighted_counts - np.min(all_weighted_counts))/
    #                  (np.max(all_weighted_counts) - np.min(all_weighted_counts)), cmap=my_cmap)
    img = plt.imshow(all_weighted_counts / np.max(all_weighted_counts), cmap=my_cmap)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts

def show_weighted_survival_count(ensemble,
                                 my_dir,
                                 screen=False,
                                 filename="weighted_survival_count.png"):
    """
    """
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        _, weighted_counts, _, _ = count_persistence(ensemble[:, :, i])
        all_weighted_counts += weighted_counts

    plt.figure()
    # img = plt.imshow((all_weighted_counts - np.min(all_weighted_counts))/
    #                  (np.max(all_weighted_counts) - np.min(all_weighted_counts)), cmap=my_cmap)
    img = plt.imshow(
        (all_weighted_counts - np.min(all_weighted_counts))
        / (np.max(all_weighted_counts) - np.min(all_weighted_counts)),
        cmap=my_cmap,
    )
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts


def show_weighted_instability_count(ensemble,
                                    my_dir,
                                    screen=False,
                                    filename="weighted_instability_count.png"):
    """
    """
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        _, _, weighted_counts, _ = count_persistence(ensemble[:, :, i])
        all_weighted_counts += weighted_counts

    plt.figure()
    img = plt.imshow(
        (all_weighted_counts - np.min(all_weighted_counts))
        / (np.max(all_weighted_counts) - np.min(all_weighted_counts)),
        cmap=my_cmap,
    )
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts


def show_weighted_consumption_count(ensemble,
                                    my_dir,
                                    screen=False,
                                    filename="weighted_consumption_count.png"):
    """
        Visualize Survival Map using weighted consumption count technique.
    """
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        # count_persistence function computes survival for single instance and is
        # defined in utils.py
        _, _, _, weighted_counts = count_persistence(ensemble[:, :, i])
        all_weighted_counts += weighted_counts

    plt.figure()
    # img = plt.imshow((all_weighted_counts - np.min(all_weighted_counts))/
    #                  (np.max(all_weighted_counts) - np.min(all_weighted_counts)), cmap=my_cmap)
    img = plt.imshow(all_weighted_counts / np.max(all_weighted_counts), cmap=my_cmap)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts

# Newly added by Tushar
def visualize_variance_weighted_consumption_count(ensemble,
                                    my_dir,
                                    screen=False,
                                    filename="variance_weighted_consumption_count.png"):
    """
        Visualize Survival Map using weighted consumption count technique.
        """
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    weighted_counts_ensemble = np.zeros(ensemble[:,:,:].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        # count_persistence function computes survival for single instance and is
        # defined in utils.py
        #_, _, _, weighted_counts = count_persistence(ensemble[:, :, i])
        weighted_counts = np.squeeze(tushar_count_persistence(ensemble[:, :, i]))
        weighted_counts_ensemble[:,:,i] = weighted_counts;
        #all_weighted_counts += weighted_counts

    # Compute variance
    all_weighted_counts = np.var(weighted_counts_ensemble, axis=2)
    plt.figure()
    # img = plt.imshow((all_weighted_counts - np.min(all_weighted_counts))/
    #                  (np.max(all_weighted_counts) - np.min(all_weighted_counts)), cmap=my_cmap)
    img = plt.imshow(all_weighted_counts / np.max(all_weighted_counts), cmap=my_cmap)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    if ensemble.shape[1] > 2 * ensemble.shape[0]:
        plt.colorbar(img, orientation="horizontal")
    else:
        plt.colorbar(img, orientation="vertical")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts

def show_max_consumption(ensemble,
                         my_dir,
                         screen=False,
                         filename="max_consumption.png"):
    """
    """
    all_weighted_counts = np.zeros(ensemble[:, :, 0].shape)
    for i in range(ensemble.shape[2]):
        # for i in range(1):
        weighted_counts = max_consumption(ensemble[:, :, i])
        all_weighted_counts += weighted_counts

    plt.figure()
    img = plt.imshow(all_weighted_counts / np.max(all_weighted_counts), cmap=cm.viridis)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_ylim(0, all_weighted_counts.shape[0])
    plt.colorbar(img, orientation="horizontal")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()
    return all_weighted_counts


def show_median_counts(ensemble,
                       my_dir,
                       screen=False,
                       filename="_median.png"):
    """
    """
    survival_counts = np.zeros(ensemble.shape)
    weighted_survival_counts = np.zeros(ensemble.shape)
    weighted_instability_counts = np.zeros(ensemble.shape)
    weighted_consumption_counts = np.zeros(ensemble.shape)
    for i in range(ensemble.shape[2]):
        survival_counts[:, :, i], weighted_survival_counts[
            :, :, i
        ], weighted_instability_counts[:, :, i], weighted_consumption_counts[
            :, :, i
        ] = count_persistence(
            ensemble[:, :, i]
        )

    for name, counts in zip(
        [
            "survival_count",
            "weighted_survival_count",
            "weighted_instability_count",
            "weighted_consumption_count",
        ],
        [
            survival_counts,
            weighted_survival_counts,
            weighted_instability_counts,
            weighted_consumption_counts,
        ],
    ):
        plt.figure()
        img = plt.imshow(
            (np.median(counts, axis=2) - np.max(counts))
            / (np.max(counts) - np.max(counts)),
            cmap=my_cmap,
        )
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().set_ylim(0, counts.shape[0])
        if ensemble.shape[0] > 2 * ensemble.shape[1]:
            plt.colorbar(img, orientation="horizontal")
        else:
            plt.colorbar(img, orientation="vertical")
        plt.savefig("{}/{}".format(my_dir, name + filename), bbox_inches="tight")
        if screen:
            plt.show()
        plt.close()


def show_variance(counts,
                  my_dir,
                  screen=False,
                  filename="weighted_surival_count_variance.png"):
    """
    """
    mean_images = []
    max_radius = 5

    image = (counts - np.min(counts)) / (np.max(counts) - np.min(counts))
    eps = 1e-16
    for i in range(1, max_radius):
        mean_images.append(filters.rank.mean(image, selem=morphology.disk(i)))
        if screen:
            show_image(mean_images[-1] + eps)
            plt.gca().set_title("Mean r={}".format(i))
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)
            plt.gca().set_ylim(0, counts.shape[0])
            plt.show()
            plt.close()

    image = 255 * image
    variance_images = []
    eps = 1e-16
    for i, mean_image in enumerate(mean_images):
        variance_images.append(np.power(image - mean_image, 2) + eps)
        show_image(variance_images[-1])
        plt.gca().set_title("Variance r={}".format(i))
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().set_ylim(0, counts.shape[0])
        if screen:
            plt.show()
        plt.savefig("{}/{}_{}".format(my_dir, i, filename), bbox_inches="tight")
        plt.close()


def show_probabilities_colormap(ensemble,
                                assignments,
                                my_dir,
                                screen=False,
                                filename="partition_probabilities.png"):
    """
        Visualize cluster membership probabilities for domain in separate subplot, where each subplot represents one cluster ID.
        
        ensemble: ensemble dataset
        assignments: function definition representing funtion assign_labels (in utils.py) with fixed parameter values for maxima_map, n_clusters, and persistence. See example inside analyze function in pipeline.py
        my_dir: Directory to store result
        screen: Whether show result on screen or not
        filename: name of rusult file
        
    """
    ps = []
    fields = []
    
    # Fields has a size same as ensemble after running the loop below, in which each plane of fields represents partitions of domain represented by cluster IDs
    for i in range(ensemble.shape[2]):
        field, p = assignments(ensemble[:, :, i])
        ps.append(p)
        fields.append(field)

    ps = np.array(ps)
    fields = np.array(fields)

    num_partitions = len(np.unique(fields[0]))
    shape = (num_partitions,) + fields[0].shape
    label_images = np.zeros(shape)

    dim1 = int(np.ceil(np.sqrt(num_partitions)))
    dim2 = int(np.floor(num_partitions / dim1 + 0.5))

    # 2D array of subplots to plot partition probabilities for single cluster
    fig, axes = plt.subplots(dim1, dim2, tight_layout=True)

    axes = axes.flatten()
    # Visualize in separate subplot each Morse complex 2-cell by looping through cell ID
    for i in range(num_partitions):
        # test image is a 3D array where fields has cluster ID i
        test_image = fields == i
        # Sum of 3D array along x axim to get a 2D array representing
        # counts/probabilities of position (x,y) belonging to cluster ID i
        # if hXw is size of a 2D array, label images has size= num_partitionsXhXW
        label_images[i] = np.sum(test_image, axis=0)
        
        #label_images[i] is a 2D image
        img = axes[i].imshow(label_images[i])
        axes[i].set_ylim(0, ensemble.shape[0])
        axes[i].get_xaxis().set_visible(False)
        axes[i].get_yaxis().set_visible(False)
    for i in range(num_partitions, dim1 * dim2):
        axes[i].set_visible(False)

    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_blended_overlay(ensemble,
                         assignments,
                         my_dir,
                         gamma=2.2,
                         screen=False,
                         filename="uncertain_assignment_blended_overlay.png"):
    """
        Visualization of Probabilistic Map
        
        Args:
            ensemble: Ensemble dataset
            Assignments:function definition representing funtion assign_labels (in utils.py) with fixed parameter values for maxima_map, n_clusters, and persistence. See example inside analyze function in pipeline.py
            my_dir: Directory to store result
            screen: Whether show result on screen or not
            filename: name of rusult file
            
    """
    ps = []
    fields = []
    
    # count stores number of ensemble realizations. Variable will be used to compute cluster probabilities later
    count = ensemble.shape[2]
    
    # Fields has a size same as ensemble after running the loop below, in which each plane of fields represents partitions of domain represented by cluster IDs
    for i in range(count):
        field, p = assignments(ensemble[:, :, i])
        ps.append(p)
        fields.append(field)

    ps = np.array(ps)
    fields = np.array(fields)

    num_partitions = len(np.unique(fields[0]))
    shape = (num_partitions,) + fields[0].shape
    label_images = np.zeros(shape)

    for i in range(num_partitions):
        # Collect positions (x,y,z) which flow to cluster ID i across all ensemble realizations
        test_image = fields == i
        # Sum all planes of test_image to understand how many times position (x,y) flows to cluster ID relative to count variable
        label_images[i] = np.sum(test_image, axis=0)

    # Store color image for each plane representing probability of flowing to cluster ID i
    colored_images = []
    for i, c in zip(range(num_partitions), ccycle):
        colored_image = np.zeros(label_images[0].shape + (4,))
        # For each cluster ID, color is fixed below. It is stored in first three channels RGB.
        colored_image[:, :, 0] = c[0] / 255.0
        colored_image[:, :, 1] = c[1] / 255.0
        colored_image[:, :, 2] = c[2] / 255.0
        colored_images.append(colored_image)

    for i, label_image in enumerate(label_images):
        # This is critical component. It denotes probability of flowing to cluster ID i.
        # It is stored in fourth channel Alpha.
        colored_images[i][:, :, 3] = label_image / count

    # Plot composit of colored image planes, where each plane denotes colormapped (through alpha channel) probability of belonging to cluster ID i
    plt.figure()
    composite_image = 255 * np.ones(colored_images[0].shape)[:, :, :-1]
    for colored_image in colored_images:
        composite_image = overlay_alpha_image_precise(
            composite_image, 255 * colored_image, 1.0, gamma
        )

    # Visualize composit image (probabilistic map)
    plt.imshow(composite_image)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    # plt.gca().set_ylim(ensemble.shape[0], 0)
    plt.gca().set_ylim(0, ensemble.shape[0])
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_contour_overlay(ensemble,
                         assignments,
                         my_dir,
                         colored=False,
                         screen=False,
                         filename="uncertain_assignment_contour_overlay.png"):
    """
        Visualizze certain regions of probabilistic map along with 0.5 probability contours
        
        Args:
        ensemble: Ensemble dataset
        Assignments:function definition representing funtion assign_labels (in utils.py) with fixed parameter values for maxima_map, n_clusters, and persistence. See example inside analyze function in pipeline.py
        my_dir: Directory to store result
        colored: flag to denote if you want to only display uncertain regions
        screen: Whether show result on screen or not
        filename: name of rusult file
    """
    ps = []
    fields = []
    
    # count stores number of ensemble realizations. Variable will be used to compute cluster probabilities later
    count = ensemble.shape[2]
    
    # Fields has a size same as ensemble after running the loop below, in which each plane of fields represents partitions of domain represented by cluster IDs
    for i in range(count):
        
        field, p = assignments(ensemble[:, :, i])
        ps.append(p)
        fields.append(field)

    ps = np.array(ps)
    fields = np.array(fields)

    num_partitions = len(np.unique(fields[0]))
    shape = (num_partitions,) + fields[0].shape
    label_images = np.zeros(shape)

    for i in range(num_partitions):
         # Collect positions (x,y,z) which flow to cluster ID i across all ensemble realizations
        test_image = fields == i
         # Sum all planes of test_image to understand how many times position (x,y) flows to cluster ID relative to count variable
        label_images[i] = np.sum(test_image, axis=0)

    # Store color image for each plane representing probability of flowing to cluster ID i
    colored_images = []
    for i, c in zip(range(num_partitions), ccycle):
        colored_image = np.zeros(label_images[0].shape + (4,))
        # For each cluster ID, color is fixed below. It is stored in first three channels RGB.
        colored_image[:, :, 0] = c[0] / 255.0
        colored_image[:, :, 1] = c[1] / 255.0
        colored_image[:, :, 2] = c[2] / 255.0
        colored_images.append(colored_image)

    for i, label_image in enumerate(label_images):
        # This is critical component. It denotes probability of flowing to cluster ID i.
        # It is stored in fourth channel Alpha.
        colored_images[i][:, :, 3] = label_image / count

    plt.figure()
    for i, color in zip(range(num_partitions), ccycle):
        my_color = "#{:>02}{:>02}{:>02}".format(*[hex(c).split("x")[-1] for c in color])
        if colored:
            plt.contourf(
                colored_images[i][:, :, 3], levels=[1e-6, 1], colors=my_color, alpha=0.5
            )
        else:
            # Visualize certain regions
            plt.contourf(
                colored_images[i][:, :, 3],
                levels=[0.99999, 1],
                colors=my_color,
                alpha=0.5,
            )
        # Plot contours for each cluster
        plt.contour(
            colored_images[i][:, :, 3],
            levels=[0.0, 0.5, 1],
            colors=my_color,
            linewidths=[1, 0.5, 1.0],
            linestyles=["solid", "dashed", "solid"],
        )
    plt.gca().set_xlim(0, ensemble.shape[1])
    # plt.gca().set_ylim(ensemble.shape[0], 0)
    plt.gca().set_ylim(0, ensemble.shape[0])
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_aspect("equal")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_certain_regions(ensemble,
                         assignments,
                         my_dir,
                         contours=False,
                         screen=False,
                         filename="certain_assignment.png"):
    """
        Visualize certain regions of Probabilistic Map
        
        Args:
        ensemble: Ensemble dataset
        Assignments:function definition representing funtion assign_labels (in utils.py) with fixed parameter values for maxima_map, n_clusters, and persistence. See example inside analyze function in pipeline.py
        my_dir: Directory to store result
        screen: Whether show result on screen or not
        filename: name of rusult file
        
    """
    ps = []
    fields = []
    
     # count stores number of ensemble realizations. Variable will be used to compute cluster probabilities later
    count = ensemble.shape[2]
    
    # Fields has a size same as ensemble after running the loop below, in which each plane of fields represents partitions of domain represented by cluster IDs
    for i in range(count):
        field, p = assignments(ensemble[:, :, i])
        ps.append(p)
        fields.append(field)

    ps = np.array(ps)
    fields = np.array(fields)

    num_partitions = len(np.unique(fields[0]))
    shape = (num_partitions,) + fields[0].shape
    label_images = np.zeros(shape)

    for i in range(num_partitions):
        # Collect positions (x,y,z) which flow to cluster ID i across all ensemble realizations
        test_image = fields == i
        
        # Sum all planes of test_image to understand how many times position (x,y) flows to cluster ID i relative to count variable
        label_images[i] = np.sum(test_image, axis=0)

    # Store color image for each plane representing probability of flowing to cluster ID i
    colored_images = []
    for i, c in zip(range(num_partitions), ccycle):
        colored_image = np.zeros(label_images[0].shape + (4,))
        # For each cluster ID, color is fixed below. It is stored in first three channels RGB.
        colored_image[:, :, 0] = c[0] / 255.0
        colored_image[:, :, 1] = c[1] / 255.0
        colored_image[:, :, 2] = c[2] / 255.0
        colored_images.append(colored_image)

    for i, label_image in enumerate(label_images):
        # This is critical component. It denotes probability of flowing to cluster ID i.
        # It is stored in fourth channel Alpha.
        colored_images[i][:, :, 3] = label_image / count

    # Visualize certain regions of colored image planes, where each plane denotes colormapped (through alpha channel) probability of belonging to cluster ID i
    plt.figure()
    for i, color in zip(range(num_partitions), ccycle):
        my_color = "#{:>02}{:>02}{:>02}".format(*[hex(c).split("x")[-1] for c in color])
        # Draw filled contours for isovlue of 1 for alpha channel, i.e., certain regions
        plt.contourf(
            colored_images[i][:, :, 3], levels=[0.99999, 1], colors=my_color, alpha=0.5
        )
        # if flag for contours is set, draw contours of 3 probability values (not filled)
        # 0.5 probability contour essentially represent expected Morse complex 1-cells
        if contours:
            plt.contour(
                colored_images[i][:, :, 3],
                levels=[0.0, 0.5, 1],
                colors=my_color,
                linewidths=[1, 0.5, 1.0],
                linestyles=["solid", "dashed", "solid"],
            )
        # else:
        #     plt.contour(colored_images[i][:, :, 3], levels=[0.0, 1], colors=my_color, linewidths=[
        #                 1, 1.0], linestyles=['solid', 'dashed', 'solid'])
    plt.gca().set_xlim(0, ensemble.shape[1])
    # plt.gca().set_ylim(ensemble.shape[0], 0)
    plt.gca().set_ylim(0, ensemble.shape[0])
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_aspect("equal")
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_combined_overlay(ensemble,
                          assignments,
                          my_dir,
                          gamma=2.2,
                          contours=False,
                          screen=False,
                          filename="uncertain_region_assignments.png"):
    """
        Visualization of uncertain regions of Probabilistic Map
        
        Args:
        ensemble: Ensemble dataset
        Assignments:function definition representing funtion assign_labels (in utils.py) with fixed parameter values for maxima_map, n_clusters, and persistence. See example inside analyze function in pipeline.py
        my_dir: Directory to store result
        contours: Visualize contours or not
        screen: Whether show result on screen or not
        filename: name of rusult file
    """
    ps = []
    fields = []
    
     # count stores number of ensemble realizations. Variable will be used to compute cluster probabilities later
    count = ensemble.shape[2]
    
    # Fields has a size same as ensemble after running the loop below, in which each plane of fields represents partitions of domain represented by cluster IDs
    for i in range(count):
        field, p = assignments(ensemble[:, :, i])
        ps.append(p)
        fields.append(field)

    ps = np.array(ps)
    fields = np.array(fields)

    num_partitions = len(np.unique(fields[0]))
    shape = (num_partitions,) + fields[0].shape
    label_images = np.zeros(shape)

    for i in range(num_partitions):
         # Collect positions (x,y,z) which flow to cluster ID i across all ensemble realizations
        test_image = fields == i
         # Sum all planes of test_image to understand how many times position (x,y) flows to cluster ID i relative to count variable
        label_images[i] = np.sum(test_image, axis=0)

     # Store color image for each plane representing probability of flowing to cluster ID i
    colored_images = []
    for i, c in zip(range(num_partitions), ccycle):
        colored_image = np.zeros(label_images[0].shape + (4,))
         # For each cluster ID, color is fixed below. It is stored in first three channels RGB.
        colored_image[:, :, 0] = c[0] / 255.0
        colored_image[:, :, 1] = c[1] / 255.0
        colored_image[:, :, 2] = c[2] / 255.0
        colored_images.append(colored_image)

    for i, label_image in enumerate(label_images):
        # This is critical component. It denotes probability of flowing to cluster ID i.
        # It is stored in fourth channel Alpha.
        colored_images[i][:, :, 3] = label_image / count

     # Visualize uncertain regions of colored image planes, where each plane denotes colormapped (through alpha channel) probability of belonging to cluster ID i
    plt.figure()

    # Generate a composite image (probabilistic map)
    composite_image = 255 * np.ones(colored_images[0].shape)[:, :, :-1]
    for colored_image in colored_images:
        composite_image = overlay_alpha_image_precise(
            composite_image, 255 * colored_image, 1.0, gamma
        )
    for i, color in zip(range(num_partitions), ccycle):
        my_color = "#{:>02}{:>02}{:>02}".format(*[hex(c).split("x")[-1] for c in color])
        # Fill certain portions of image, i.e., alpha = 1, using white color. Thus we are left
        # with visualization of uncertain regions.
        plt.contourf(
            colored_images[i][:, :, 3], levels=[0.99999, 1], colors="#FFFFFF", alpha=1
        )
        
        # Visualize 0.5 probability contours, i.e., expected Morse complex boundaries.
        if contours:
            plt.contour(
                colored_images[i][:, :, 3],
                levels=[0.5],
                colors="#000000",
                linewidths=[1.0],
                linestyles=["solid"],
            )

    # visualize composite image with certain regions filled with white color
    plt.imshow(composite_image)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    # plt.gca().set_ylim(ensemble.shape[0], 0)
    plt.gca().set_ylim(0, ensemble.shape[0])
    plt.savefig("{}/{}".format(my_dir, filename), bbox_inches="tight")
    if screen:
        plt.show()
    plt.close()


def show_label_boundaries(labels, ax, color="#000000"):
    """
    """    
    h, w = labels.shape

    keyMap = {}
    levels = []
    for i, k in enumerate(np.unique(labels)):
        keyMap[k] = i
        levels.append(i + 0.5)

    color_mesh = np.zeros((h, w))
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            color_mesh[i, j] = keyMap[labels[i, j]]

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    for i in keyMap.values():
        ax.contour(color_mesh == i, colors=color, levels=levels, linewidths=1)
    ax.set_aspect("equal")
    # draw(filename)
