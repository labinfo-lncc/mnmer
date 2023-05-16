# %%
# Common imports
import numpy as np
import os
import pandas as pd
import time
import sys

# %%
FULL_FILE = "/home/aandrade/mnmer/mockup/mockup_13mer.csv"
print(FULL_FILE)

# %%
features = pd.read_csv(FULL_FILE)
features.head()

# %%
(n_obs, n_feat) = features.shape
print("num of contigs, num of features")
(n_obs, n_feat)

# %%
data   = features.drop("seqid", axis=1) # drop labels for training set
N_data = len(data)
label  = features["seqid"].copy()

# %%
for i in range(len(label)):
    target = label[i]
    position = target.find("_")
    if position > 0:
        new_target = target[:position]
        label[i] = new_target       

# %%
n_trueclasses = np.unique(label).size
print("num of real clusters") 
n_trueclasses

# %%
print("num of contigs per clust")
label.value_counts()

# %%
SIZE = 20000
SIZE = np.min([SIZE,N_data])
X = data[:SIZE]
y = label[:SIZE]

# %% [markdown]
# AGGLOMERATIVE CLUSTERING
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics

start = time.time()
model = AgglomerativeClustering(n_clusters=n_trueclasses, linkage="ward")
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="agglomerative.clust"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

mainMetrics = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette

from sklearn.cluster import KMeans

start = time.time()
model = KMeans(n_clusters=n_trueclasses, init='k-means++', n_init='auto')
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="KMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics1 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics1, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics1

from sklearn.cluster import MiniBatchKMeans

start = time.time()
model = MiniBatchKMeans(n_clusters=n_trueclasses, n_init='auto')
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="MiniBatchKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics2 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics2, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics2

from sklearn.cluster import BisectingKMeans

start = time.time()
model = BisectingKMeans(n_clusters=n_trueclasses, random_state=0)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="BisectingKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics3 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics3, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics3

from sklearn.cluster import SpectralClustering

start = time.time()
model = SpectralClustering(n_clusters=n_trueclasses)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="Spectral"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics4 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics4, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics4

from sklearn.mixture import GaussianMixture

start = time.time()
model = GaussianMixture(n_components=10)
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="GaussianMixture"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics5 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics5, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics5

from sklearn.cluster import AffinityPropagation

start = time.time()
model = AffinityPropagation(damping=0.8)
yhat  = model.fit_predict(X)
cluster_centers_indices = model.cluster_centers_indices_
n_clusters_ = len(cluster_centers_indices)
duration = time.time()-start

db = "mockup"
mn = "{1,3}-mer"
algorithm="AffinityPropagation"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics6 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics6, mainMetrics])
mainMetrics.to_csv("metrics_mockup_13mer.csv", sep=',')

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics6,FULL_FILE,data,label,X,y

# %%
FULL_FILE = "/home/aandrade/mnmer/mockup/mockup_22mer.csv"
print(FULL_FILE)

# %%
features = pd.read_csv(FULL_FILE)
features.head()

# %%
(n_obs, n_feat) = features.shape
print("num of contigs, num of features")
(n_obs, n_feat)

# %%
data   = features.drop("seqid", axis=1) # drop labels for training set
N_data = len(data)
label  = features["seqid"].copy()

# %%
for i in range(len(label)):
    target = label[i]
    position = target.find("_")
    if position > 0:
        new_target = target[:position]
        label[i] = new_target       

# %%
n_trueclasses = np.unique(label).size
print("num of real clusters") 
n_trueclasses

# %%
print("num of contigs per clust")
label.value_counts()

# %%
SIZE = 20000
SIZE = np.min([SIZE,N_data])
X = data[:SIZE]
y = label[:SIZE]

# %% [markdown]
# AGGLOMERATIVE CLUSTERING
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics

start = time.time()
model = AgglomerativeClustering(n_clusters=n_trueclasses, linkage="ward")
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="agglomerative.clust"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

mainMetrics = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette

from sklearn.cluster import KMeans

start = time.time()
model = KMeans(n_clusters=n_trueclasses, init='k-means++', n_init='auto')
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="KMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics1 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics1, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics1

from sklearn.cluster import MiniBatchKMeans

start = time.time()
model = MiniBatchKMeans(n_clusters=n_trueclasses, n_init='auto')
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="MiniBatchKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics2 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics2, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics2

from sklearn.cluster import BisectingKMeans

start = time.time()
model = BisectingKMeans(n_clusters=n_trueclasses, random_state=0)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="BisectingKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics3 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics3, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics3

from sklearn.cluster import SpectralClustering

start = time.time()
model = SpectralClustering(n_clusters=n_trueclasses)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="Spectral"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics4 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics4, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics4

from sklearn.mixture import GaussianMixture

start = time.time()
model = GaussianMixture(n_components=10)
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="GaussianMixture"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics5 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics5, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics5

from sklearn.cluster import AffinityPropagation

start = time.time()
model = AffinityPropagation(damping=0.8)
yhat  = model.fit_predict(X)
cluster_centers_indices = model.cluster_centers_indices_
n_clusters_ = len(cluster_centers_indices)
duration = time.time()-start

db = "mockup"
mn = "{2,2}-mer"
algorithm="AffinityPropagation"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics6 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics6, mainMetrics])
mainMetrics.to_csv("metrics_mockup_22mer.csv", sep=',')

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics6,FULL_FILE,data,label,X,y

# %%
FULL_FILE = "/home/aandrade/mnmer/mockup/mockup_31mer.csv"
print(FULL_FILE)

# %%
features = pd.read_csv(FULL_FILE)
features.head()

# %%
(n_obs, n_feat) = features.shape
print("num of contigs, num of features")
(n_obs, n_feat)

# %%
data   = features.drop("seqid", axis=1) # drop labels for training set
N_data = len(data)
label  = features["seqid"].copy()

# %%
for i in range(len(label)):
    target = label[i]
    position = target.find("_")
    if position > 0:
        new_target = target[:position]
        label[i] = new_target       

# %%
n_trueclasses = np.unique(label).size
print("num of real clusters") 
n_trueclasses

# %%
print("num of contigs per clust")
label.value_counts()

# %%
SIZE = 20000
SIZE = np.min([SIZE,N_data])
X = data[:SIZE]
y = label[:SIZE]

# %% [markdown]
# AGGLOMERATIVE CLUSTERING
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics

start = time.time()
model = AgglomerativeClustering(n_clusters=n_trueclasses, linkage="ward")
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="agglomerative.clust"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

mainMetrics = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette

from sklearn.cluster import KMeans

start = time.time()
model = KMeans(n_clusters=n_trueclasses, init='k-means++', n_init='auto')
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="KMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics1 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics1, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics1

from sklearn.cluster import MiniBatchKMeans

start = time.time()
model = MiniBatchKMeans(n_clusters=n_trueclasses, n_init='auto')
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="MiniBatchKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics2 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics2, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics2

from sklearn.cluster import BisectingKMeans

start = time.time()
model = BisectingKMeans(n_clusters=n_trueclasses, random_state=0)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="BisectingKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics3 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics3, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics3

from sklearn.cluster import SpectralClustering

start = time.time()
model = SpectralClustering(n_clusters=n_trueclasses)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="Spectral"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics4 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics4, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics4

from sklearn.mixture import GaussianMixture

start = time.time()
model = GaussianMixture(n_components=10)
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="GaussianMixture"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics5 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics5, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics5

from sklearn.cluster import AffinityPropagation

start = time.time()
model = AffinityPropagation(damping=0.8)
yhat  = model.fit_predict(X)
cluster_centers_indices = model.cluster_centers_indices_
n_clusters_ = len(cluster_centers_indices)
duration = time.time()-start

db = "mockup"
mn = "{3,1}-mer"
algorithm="AffinityPropagation"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics6 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics6, mainMetrics])
mainMetrics.to_csv("metrics_mockup_31mer.csv", sep=',')

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics6,FULL_FILE,data,label,X,y

# %%
FULL_FILE = "/home/aandrade/mnmer/mockup/mockup_40mer.csv"
print(FULL_FILE)

# %%
features = pd.read_csv(FULL_FILE)
features.head()

# %%
(n_obs, n_feat) = features.shape
print("num of contigs, num of features")
(n_obs, n_feat)

# %%
data   = features.drop("seqid", axis=1) # drop labels for training set
N_data = len(data)
label  = features["seqid"].copy()

# %%
for i in range(len(label)):
    target = label[i]
    position = target.find("_")
    if position > 0:
        new_target = target[:position]
        label[i] = new_target       

# %%
n_trueclasses = np.unique(label).size
print("num of real clusters") 
n_trueclasses

# %%
print("num of contigs per clust")
label.value_counts()

# %%
SIZE = 20000
SIZE = np.min([SIZE,N_data])
X = data[:SIZE]
y = label[:SIZE]

# %% [markdown]
# AGGLOMERATIVE CLUSTERING
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics

start = time.time()
model = AgglomerativeClustering(n_clusters=n_trueclasses, linkage="ward")
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="agglomerative.clust"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

mainMetrics = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette

from sklearn.cluster import KMeans

start = time.time()
model = KMeans(n_clusters=n_trueclasses, init='k-means++', n_init='auto')
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="KMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics1 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics1, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics1

from sklearn.cluster import MiniBatchKMeans

start = time.time()
model = MiniBatchKMeans(n_clusters=n_trueclasses, n_init='auto')
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="MiniBatchKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics2 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics2, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics2

from sklearn.cluster import BisectingKMeans

start = time.time()
model = BisectingKMeans(n_clusters=n_trueclasses, random_state=0)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="BisectingKMeans"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics3 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics3, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics3

from sklearn.cluster import SpectralClustering

start = time.time()
model = SpectralClustering(n_clusters=n_trueclasses)
yhat = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="Spectral"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics4 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics4, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics4

from sklearn.mixture import GaussianMixture

start = time.time()
model = GaussianMixture(n_components=10)
yhat  = model.fit_predict(X)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="GaussianMixture"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics5 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics5, mainMetrics])

del db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette,metrics5

from sklearn.cluster import AffinityPropagation

start = time.time()
model = AffinityPropagation(damping=0.8)
yhat  = model.fit_predict(X)
cluster_centers_indices = model.cluster_centers_indices_
n_clusters_ = len(cluster_centers_indices)
duration = time.time()-start

db = "mockup"
mn = "{0,4}-mer"
algorithm="AffinityPropagation"
homogeneity = metrics.homogeneity_score(y,yhat)
completeness = metrics.completeness_score(y,yhat)
vmeasure = metrics.v_measure_score(y,yhat)
adjrand = metrics.adjusted_rand_score(y,yhat)
adjmutual = metrics.adjusted_mutual_info_score(y,yhat)
silhouette= metrics.silhouette_score(X, yhat, metric="sqeuclidean")

metrics6 = pd.DataFrame([[db,mn,algorithm,homogeneity,completeness,vmeasure,adjrand,adjmutual,silhouette]],columns=['db','mn','algorithm','homogeneity','commpleteness','vmeasure','adjrand','adjmutual','silhouette'])

mainMetrics = pd.concat([metrics6, mainMetrics])
mainMetrics.to_csv("metrics_mockup_40mer.csv", sep=',')


