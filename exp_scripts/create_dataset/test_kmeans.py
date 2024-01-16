from sklearn.cluster import KMeans
import numpy as np

data_values = [
204,
269,
236,
336,
222,
492,
261,
219,
320,
245,
393,
207,
333,
401,
215,
183,
421,
225,
326,
228]

data_reshaped = np.array(data_values).reshape(-1,1)


# Create a KMeans object
kmeans = KMeans(n_clusters=3)

# Fit the model to the data
kmeans.fit(data_reshaped)

# Predict the cluster labels for each data point
labels = kmeans.predict(data_reshaped)

# Print the cluster labels
print(labels)