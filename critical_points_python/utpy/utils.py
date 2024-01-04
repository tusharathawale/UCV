from glob import glob
import scipy.io
import sklearn.cluster
from sklearn.preprocessing import MinMaxScaler
import numpy as np

import nglpy as ngl
import topopy
import flatpy
import math

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

graph_params = {
    "index": None,
    "max_neighbors": 10,
    "relaxed": False,
    "beta": 1,
    "p": 2.0,
}


def load_data(foo="ackley", noise_level=0.3):
    """ TODO
    
    Args:
    Returns:
        None

    """
    assignment_map = {
        "4peaks": assignments_4peaks,
        "ackley": assignments_ackley,
        "salomon": assignments_salomon,
    }
    base_name = "data/" + foo

    ground_truth = scipy.io.loadmat(base_name + "/groundtruth.mat")["gt"]
    uncertain_realizations = scipy.io.loadmat(
        "{}/{}_uncertain.mat".format(base_name, str(noise_level))
    )["noisyEnsemble"]

    assignments = assignment_map[foo]
    return ground_truth, uncertain_realizations, assignments


def load_ensemble(name="matVelocity"):
    """
        Load an ensemble of realizations, where each realization is stored as a MAT file.
        Note: Make sure that the name of the MAT file is same as the name of stored matrix.
        
    Args:
        name: Directory name in which realizations (MAT files) are stored.
    
    Returns:
        A tuple of realizations

    """
    base_name = "data/" + name
    files = glob("{}/*.mat".format(base_name))
    uncertain_realizations = []
    for filename in files:
        token = filename.rsplit("/", 1)[1].split(".")[0]
        # Temperature data is a bit special in that it has multiple
        # entries, let's just take the first for now.
        if name == "matTemperature":
            uncertain_realizations.append(scipy.io.loadmat(filename)[token][:, :, 0].T)
        else:
            uncertain_realizations.append(scipy.io.loadmat(filename)[token].T)
    return np.array(uncertain_realizations).T

def massage_data(grid):
    """ Takes 2D scalar field as input, and generates a 2D array X and 1D array Y. Array X holds index positions for grid vertices, and Array Y holds function values for corresponding index positions in X.
    
    Args:
        grid: 2D array representing scalar field
        
    Returns:
        X: a 2D array of (x,y) pairs: [(x1,y1), (x2,y2),..]
        Y: data values for each index pair in array X

    """
    X = []
    Y = []
    for row, vals in enumerate(grid):
        for col, val in enumerate(vals):
            X.append([col, row])
            Y.append(val)
    return np.array(X, dtype=np.float64), np.array(Y, dtype=np.float64)

def count_persistence(grid):
    """
    The function for computing Survival Map in the submitted paper using different techniques.
    Args:
        grid: A 2D scalar field representing single instance of ensemble
        
    Returns:
        consumption_counts: Returns survival counts in a consumption_counts array denoting survival of local gradient flows upon persistence simplification. Other survival metrics such as plain counts etc were also returned and could be potentially useful.
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

    # Sorted heirarchy of max-saddle pairs based on persistence values
    # p: persistence value
    # k: dying index
    # x: surviving index
    partitions = tmc.get_partitions()
    sorted_hierarchy = sorted(
        [(p, k, x) for k, (p, x, s) in tmc.merge_sequence.items()]
    )

    field = np.zeros(Y.shape, dtype=int)
    for k, v in partitions.items():
        field[v] = k

    counts = np.zeros(Y.shape)
    weighted_counts = np.zeros(Y.shape)
    unstable_counts = np.zeros(Y.shape)
    consumption_counts = np.zeros(Y.shape)
    last_persistence = 0
    for persistence, dying_index, surviving_index in sorted_hierarchy:
        next_field = np.array(field)
        next_field[np.where(next_field == dying_index)] = surviving_index

        counts[np.where(field == next_field)] += 1
        weighted_counts[np.where(field == next_field)] += persistence - last_persistence
        unstable_counts[np.where(field != next_field)] += persistence - last_persistence
        consumption_counts[np.where(field == surviving_index)] += (
            persistence - last_persistence
        )
        field = next_field
        last_persistence = persistence

    return (
        counts.reshape(grid.shape),
        weighted_counts.reshape(grid.shape),
        unstable_counts.reshape(grid.shape),
        consumption_counts.reshape(grid.shape),
    )

def tushar_count_persistence(grid):
    """
        The function for computing Survival Map in the submitted paper using different techniques.
        Args:
        grid: A 2D scalar field representing single instance of ensemble
        
        Returns:
        consumption_counts: Returns survival counts in a consumption_counts array denoting survival of local gradient flows upon persistence simplification. Other survival metrics such as plain counts etc were also returned and could be potentially useful.
        """
    
    # grid is a 2D array of hXw dimension representing discrete scalar field
    # X is a 2D array of pairs [(x1,y1), (x2,y2), ..]
    # Y is 1D array of function values for indices in array in X
    X, Y = massage_data(grid)
    h, w = grid.shape
    
    #print('hello')
    # Build morse complex for discrete field
    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    tmc.build(X, Y)
    
    # Sorted heirarchy of max-saddle pairs based on persistence values
    # p: persistence value
    # k: dying index
    # x: surviving index
    partitions = tmc.get_partitions()
    sorted_hierarchy = sorted(
                              [(p, k, x) for k, (p, x, s) in tmc.merge_sequence.items()]
                              )
        
    field = np.zeros(Y.shape, dtype=int)
    for k, v in partitions.items():
        field[v] = k

    counts = np.zeros(Y.shape)
    last_persistence = 0
    
    # Loop through max-saddle pairs such that their persistence values are in increasing order
    for persistence, dying_index, surviving_index in sorted_hierarchy:
        
        print('persistence',persistence)
        print('dying index', dying_index)
        print('surviving index', surviving_index)
        
        # Initialize next_field to be same as input field (Morse complex partitions)
        next_field = np.array(field)
        
        # Apply one level of topological simplification to next_field
        next_field[np.where(next_field == dying_index)] = surviving_index

        # New implementation based on reviewer's comment: The gradient flow reversal count is equal to the depth in the tree that represents merger of Morse complex cells.
        #if(surviving_index!=dying_index):
        #    counts[np.where((field == dying_index) & (counts == 0))] = persistence
        #else:
        #    counts[np.where((field == dying_index) & (counts == 0))] = 1.1*last_persistence
        
        # Increment counts at positions where gradient flows did not change their direction upon simplification
        #counts[np.where(field == next_field)] += 1
        
        #counts[np.where(field != dying_index)] += 1
        #counts[np.where(field == dying_index)] += 1
        
        #counts[np.where(field != dying_index)] += persistence - last_persistence
        
        #if np.count_nonzero(counts[np.where(field == surviving_index)]) > 0:
        #    counts[np.where(field == surviving_index)] += persistence -last_persistence
        #else:
        #    counts[np.where(field == surviving_index)] += persistence
        
        #weighted_counts[np.where(field == next_field)] += persistence - last_persistence
        #unstable_counts[np.where(field != next_field)] += persistence - last_persistence

        #counts[np.where(field == surviving_index)] += (persistence-last_persistence)
        counts[np.where(field == surviving_index)] += persistence

        field = next_field
        last_persistence = persistence

    return (
        counts.reshape(grid.shape),
        #weighted_counts.reshape(grid.shape),
        #unstable_counts.reshape(grid.shape),
        #consumption_counts.reshape(grid.shape),
    )



def max_consumption(grid):
    """
        The function for computing Survival Map using a maximum measure (could be potentially useful).
        This one was not used in the submitted paper.
    Args:
        grid: A 2D scalar field representing single instance of ensemble
    
    Returns:
        consumption_maxs: Returns maximum survival counts in a consumption_maxs array denoting maximum survival of local gradient flows upon persistence simplification.

    """
    X, Y = massage_data(grid)
    h, w = grid.shape

    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    tmc.build(X, Y)

    partitions = tmc.get_partitions()
    sorted_hierarchy = sorted(
        [(p, k, x) for k, (p, x, s) in tmc.merge_sequence.items()]
    )

    field = np.zeros(Y.shape, dtype=int)
    for k, v in partitions.items():
        field[v] = k

    consumption_maxs = np.zeros(Y.shape)
    last_persistence = 0
    for persistence, dying_index, surviving_index in sorted_hierarchy:
        next_field = np.array(field)
        next_field[np.where(next_field == dying_index)] = surviving_index
        idxs = np.where(field == surviving_index)
        consumption_maxs[idxs] = np.maximum(
            consumption_maxs[idxs], (persistence - last_persistence)
        )
        field = next_field
        last_persistence = persistence

    return consumption_maxs.reshape(grid.shape)


def get_persistence_from_count(ensemble, n_clusters):
    """
    Compute average persistence value for an ensemble given number of partitions for each instance
    
    Args:
        ensemble: An ensemble of 2D scalar fields
        n_clusters: number of partitions for Morse complex of each instance
    
    Returns:
        Average of persistence for all instances given fixed number of partitions for Morse complex

    """
    persistences = []
    for i in range(ensemble.shape[2]):
        graph = ngl.EmptyRegionGraph(**graph_params)
        tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)

        X, Y = massage_data(ensemble[:, :, i])
        tmc.build(X, Y)
        for p in tmc.persistences:
            if len(tmc.get_partitions(p).keys()) <= n_clusters:
                persistences.append(p)
                break
    return np.average(persistences)


def get_count_from_persistence(ensemble, persistence):
    """
        Get average number of Morse complex cells across ensmble members
    
    Args:
        ensemble: ensemble of 2D scalar fields
        
    Returns:
        Average number of Morse complex cells across ensemble members

    """
    max_counts = []
    for i in range(ensemble.shape[2]):
        graph = ngl.EmptyRegionGraph(**graph_params)
        tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
        X, Y = massage_data(ensemble[:, :, i])
        tmc.build(X, Y)
        max_counts.append(len(tmc.get_partitions(persistence).keys()))
    return int(np.average(max_counts))


def create_assignment_map(ensemble, n_clusters=None, persistence=None):
    """
        Gather local maxima of ensemble members (simplified such that each member has number of local maxima = n_clusters) and cluster local maxima into given number of clusters
    
    Args:
        ensemble: Ensemble of 2D scalar fields
        n_clusters: simplification level specified by number of clusters to divide local maxima into
        persistence: simplification level specified by persistence value
        
    Returns:
    maxima_map : 2D array representing cluster number for local maxima at position (x,y)
    Cluster numbers vary from [0,.., n_clusters-1]
    """
    
    # Inputting number of clusters or persistence simplification level is must
    # In Probabilistic map technique, we provide number of clusters same as number of mandatory maxima.
    if n_clusters is None and persistence is None:
        raise ValueError("Must specify either n_clusters or persistence")
    max_points = list()
    max_member = list()

    max_counts = []

    # For each ensemble realization store local maxima in array max_points
    for i in range(ensemble.shape[2]):
        
        # compute Morse complex cells
        graph = ngl.EmptyRegionGraph(**graph_params)
        tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
        X, Y = massage_data(ensemble[:, :, i])
        tmc.build(X, Y)

        # When number of Morse complex cells (n_clusters) is provided as input
        if n_clusters is not None:
            for p in tmc.persistences:
                # simplify member until number of Morse complex cells = n_clusters
                if len(tmc.get_partitions(p).keys()) <= n_clusters:
                    # keys are such that they represent indices of local maxima
                    # hence they are unique too, and represent number of Morse complex
                    # cells
                    for key in tmc.get_partitions(p).keys():
                        # Store ensmemble member ID (in max_member) and local maxima
                        # indices (max_points) for that member ID
                        max_points.append((int(X[key, 0]), int(X[key, 1])))
                        max_member.append(i)
                    break
        else:
            max_counts.append(len(tmc.get_partitions(persistence).keys()))
            for key in tmc.get_partitions(persistence).keys():
                max_points.append((int(X[key, 0]), int(X[key, 1])))
                max_member.append(i)


    if n_clusters is None and persistence is not None:
        n_clusters = int(np.average(max_counts))

    # Array of local maxima
    maxima = np.array(max_points)

    # Cluster local maxima into number of clusters specified by user using sklearn.
    # In paper, number of clusters = no of mandatory maxima
    maxima = MinMaxScaler().fit_transform(maxima)
    clustering = sklearn.cluster.MeanShift().fit(maxima)
    clustering = sklearn.cluster.MiniBatchKMeans(n_clusters=n_clusters).fit(maxima)
    # clustering = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters).fit(maxima)
    # clustering = sklearn.cluster.KMeans(n_clusters=n_clusters).fit(maxima)
    # clustering = sklearn.cluster.DBSCAN(eps=0.3, min_samples=3).fit(maxima)
    # clustering = sklearn.cluster.SpectralClustering(n_clusters=n_clusters).fit(maxima)
    unique_labels = np.unique(clustering.labels_)

    # Maxima map is a 2D array. Locations of local maxima (max_points) are assigned id of a cluster to which local maxima belong. Cluster id were computed using sclearn.
    maxima_map = {}
    for i in range(len(maxima)):
        maxima_map[max_points[i]] = clustering.labels_[i]

    return maxima_map


def assign_labels(grid, maxima_map, n_clusters=None, persistence=None):
    """
        Assign cluster label to domain positions (x,y) depending upon which cluster (x,y) flows to.
        
    Args:
        grid: 2D discrete scalar function
        maxima_map: stores mapping of local maximum to cluster ID (cluster IDs are derived through global analysis of local maxima of all ensemble members) for input n_clusters and persistence vlaue.
        
    Returns:
        field: 2D array. Each array position (x,y) has a cluster label assigned [0,..,n_clusters] to which position (x,y) flows to. The 2D array corresponds to grid (ensemble realization given as input). Note maxima_map is constant across all ensemble realizations, as it was computed through global analysis of all ensemble realizations.
    """
    
    if n_clusters is None and persistence is None:
        raise ValueError("Must specify either n_clusters or persistence")

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
    #partitions equal to input number of partitions (n_clusters)
    if n_clusters is None and persistence is not None:
        partitions = tmc.get_partitions(persistence)
    else:
        for p in tmc.persistences:
            if len(tmc.get_partitions(p).keys()) <= n_clusters:
                persistence = p
                partitions = tmc.get_partitions(p)
                break

    # Field stores cluster index to which each location (x,y) of a 2D field flows to
    # cluster indices are are stored in maxima map
    field = np.zeros(Y.shape, dtype=int)
    # k is index of local maximum of a discrete field (linearized) and v is set of all indices which flow to local maximum k. Since maxima_map(X[k]) stores a cluster ID for local maximum, cluster ID is assigned to positions field[v]
    for k, v in partitions.items():
        field[v] = maxima_map[(int(X[k, 0]), int(X[k, 1]))]

    return field.reshape(grid.shape), persistence


def generate_ensemble(foo, noise_level, count=50, noise_model="uniform"):
    """
        Generate an ensemble of uncertain 2D scalar fields for synthetic dataset.
    
    Args:
        foo: Synthetic function [it is a funtion definition]
        noise_level: Noise will be injected in range [-noise_level, noise_level]
        count: Number of ensemble realizations to be generated
        noise_model: uniform/nonparametric
        
    Returns:
        ground_truth: Synthetic function sampled on a discrete grid
        uncertain_realizations: Ensemble realizations

    """
    
    # Discrete grid
    xi = np.arange(0, 1, 0.02)
    #xi = np.arange(0, 1, 0.2)
    xv, yv = np.meshgrid(xi, xi)
    print("xv yv",xv,yv)
    # grid vertex indices in linear sequence
    X = np.vstack((xv.flatten(), yv.flatten())).T
    
    # Compute ground truth function on indices in linear sequence
    # foo defined in flatpy library
    Z = foo(X)
    
    # Build Morse complex on ground truth (on indices in linear sequence)
    # Morse complex functions defined in topopy library
    graph = ngl.EmptyRegionGraph(**graph_params)
    tmc = topopy.MorseComplex(graph=graph, gradient="steepest", normalization=None)
    tmc.build(X, Z)
    
    partitions = tmc.get_partitions()
    # Build a sorted hierarchy of maximum-saddle pairs based on persistence value for pairs
    # Each entry in sorted hierarchy is [persistence value (p), number of domain partitions for persistence value p]
    sorted_hierarchy = sorted(
        [(p, len(tmc.get_partitions(p))) for k, (p, x, s) in tmc.merge_sequence.items()]
    )

    # This is handled outside this function at the analyze_synthetic level
    # z_range = max(Z) - min(Z)
    # noise = 0.5*z_range*noise_level
    noise = noise_level
    
    # ground truth function computed on 2D discrete grid
    zv = Z.reshape(xv.shape)
    ground_truth = zv
    
    # Initialize ensemble of 2D scalar fields with zeros
    uncertain_realizations = np.zeros(shape=ground_truth.shape + (count,))

    # Set random seed
    np.random.seed(0)
    
    # Add noise for each realization depending upon noise model
    for i in range(count):
        if noise_model == "uniform":
            uncertain_realizations[:, :, i] = flatpy.utils.add_uniform_noise(
                ground_truth, noise
            )
    
        elif noise_model == "nonparametric":
            uncertain_realizations[
                :, :, i
            ] = flatpy.utils.add_nonparametric_uniform_noise(
                ground_truth, noise, 0.4, 3 * noise
            )
        elif noise_model == "variable":
            uncertain_realizations[:, :, i] = flatpy.utils.add_nonuniform_noise(
                ground_truth, noise
            )

    # Return ground truth and ensemble members
    return ground_truth, uncertain_realizations

# Minimum agreement level of number of 2-cells across members.
# Parameter can range in 0 and 1
def autotune_from_persistence(all_ps, all_counts, agreement=0.7):
    """ TODO
    
    Args:
    Returns:
        None

    """
    unique_persistences = np.array(sorted(set([p for ps in all_ps for p in ps])))
    unique_counts = np.zeros(shape=(len(unique_persistences), len(all_ps)))
    for row, p in enumerate(unique_persistences):
        for col in range(len(all_ps)):
            index = 0
            while index < len(all_ps[col]) - 1 and all_ps[col][index] < p:
                index += 1
            unique_counts[row, col] = all_counts[col][index]

    first_saved = None
    for p, c in zip(unique_persistences, unique_counts):
        mu = np.mean(c)
        sigma = np.std(c)
        counts = np.bincount(np.array(c, dtype=int))
        max_count = np.argmax(counts)
        if first_saved is None and counts[max_count] >= agreement * len(all_ps):
            first_saved = (p, max_count)
    plt.figure()

    x = unique_persistences
    y = np.mean(unique_counts, axis=1)
    sigma = np.std(unique_counts, axis=1)

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(0, sigma.max())
    lc = LineCollection(segments, cmap="viridis", norm=norm)
    # Set the values used for colormapping
    lc.set_array(sigma)
    lc.set_linewidth(4)
    line = plt.gca().add_collection(lc)
    plt.colorbar(line)
    plt.gca().set_xlim(0, np.max(unique_persistences))
    plt.gca().set_ylim(0, np.max(unique_counts))
    plt.show()

    # plt.figure()
    # df = pd.DataFrame({"Member": memberships, "Persistence": unique_persistences, "Maximum_Count": counts})
    # sns.lineplot(x="Persistence", y="Maximum_Count", data=df)
    # plt.show()

    return first_saved


def autotune_from_survival_count(counts):
    """ TODO
    
    Args:
    Returns:
        None

    """
    # Perform image-based intensity segmentation here
    clustering = sklearn.cluster.DBSCAN(eps=0.3, min_samples=3).fit(counts)
    return
