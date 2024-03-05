import sys
import os
import numpy as np
import matplotlib.pyplot as pl
from sklearn.cluster import KMeans



def clustersPlot(czas_kmeans, kroki_kmeans, clusters):
    pl.rcParams["font.family"] = "Arial"
    pl.rcParams['figure.dpi'] = 300
    pl.rcParams["legend.fancybox"] = False
    pl.rcParams["legend.framealpha"] = 1
    pl.clf()
    pl.scatter(czas_kmeans,kroki_kmeans,alpha=0.01)
    pl.plot([np.mean(c) for c in clusters if len(c)>90],[1 for c in clusters if len(c)>90],linewidth=0,c='r',marker='.')
    pl.show()


def stepsClusterization(nazwa_pliku):
    nazwa_pliku = os.path.split(nazwa_pliku)[1]
    czas_kmeans = np.loadtxt(nazwa_pliku)
    kroki_kmeans = np.array([1 for i in range(len(czas_kmeans))])
    data = np.array([czas_kmeans,kroki_kmeans]).transpose()
    clusters = [[czas_kmeans[0]]]

    for i in range(1,len(czas_kmeans)):
        if (czas_kmeans[i]-clusters[-1][-1]) < 5e9:
            clusters[-1].append(czas_kmeans[i])
        else:
            clusters.append([czas_kmeans[i]])

    filtered_clusters = [c for c in clusters if len(c)>90]
    start_stop_times = np.array([[np.min(c),np.max(c)] for c in filtered_clusters]).transpose()
    times = np.hstack((start_stop_times[0],start_stop_times[1]))

    np.savetxt(f"czas_do_analizy_{nazwa_pliku}.txt", np.sort(times))
    # times = np.insert(times,0)
    # np.savetxt("czas_do_analizy_bone_2.txt",np.sort(times)/1E12)

    return czas_kmeans, kroki_kmeans, clusters



if __name__ == '__main__':
    nazwa_pliku = sys.argv[1]
    czas_kmeans, kroki_kmeans, clusters = stepsClusterization(nazwa_pliku)
    clustersPlot(czas_kmeans, kroki_kmeans, clusters)