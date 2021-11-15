#inspiring from from: https://towardsdatascience.com/using-weighted-k-means-clustering-to-determine-distribution-centres-locations-2567646fc31d
#apresiasi dan t kasih disampaikan kpd Bp. Parawendy Indarjo, Data scientist at Bukalapak
# source data ada di github repo, https://github.com/nasriaw/clustering-rice/blob/main/lokasi_desa.csv

#1. baca source data
df_desa = read.csv('/home/nasri/Dokumen/Belajar_komputer/bahasa_R/desa_kabmalang/running_R/lokasi_desa.csv')
head(df_desa)
dim(df_desa)
df_desa[3:6,1:5]

#2.Function haversine distance 
haversine_dist = function(point1, point2) { #each argument is a numeric vector with two elements (lon, lat)
  lon1 = point1[1] 
  lat1 = point1[2]
  lon2 = point2[1]
  lat2 = point2[2]
  
  R = 6371000 #earth radius in meters
  phi1 = lat1 * pi / 180 #convert to radian
  phi2 = lat2 * pi / 180 #convert to radian
  delta_phi = (lat2 - lat1) * pi / 180
  delta_lambda = (lon2 - lon1) * pi / 180
  
  a = (sin(delta_phi/2))^2 + cos(phi1) * cos(phi2) * ((sin(delta_lambda/2))^2)
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  
  distance = R * c #haversine distance between point1 and point 2 in meters
  return(round(distance, 2))
}

#3.function to compute total sum of squares error, SSE
SSE = function(desa_df, clustering) { 
#desa_df is a dataframe with columns: desa name, longitude, latitude, luas_sawah_ha
#clustering is vector of cluster assignment of desa
  
#preparing argument
desa_longlat_mat = as.matrix.data.frame(desa_df[,2:3])
  colnames(desa_longlat_mat) = NULL
  
#compute centroid
  centroid_long = vector()
  centroid_lat = vector()
  for (k in c(1:max(clustering))) {
    cluster_k = which(clustering == k) #desa index of cluster k
    centroid_long[k] = weighted.mean(desa_df$longitude[cluster_k], desa_df$luas_sawah_ha[cluster_k])
    centroid_lat[k] = weighted.mean(desa_df$latitude[cluster_k], desa_df$luas_sawah_ha[cluster_k])
  }
  mat_centroid = cbind(centroid_long, centroid_lat)
  colnames(mat_centroid) = NULL
  
#compute SSE
  sum_squares_err = 0
  for (k in c(1:max(clustering))) {
    cluster_k = which(clustering == k) #desa index of cluster k
    for (i in c(1:length(cluster_k))) {
      for (j in c(1:2)) { #long, lat
        square_err = (desa_longlat_mat[cluster_k[i],j] - mat_centroid[k,j])^2
        sum_squares_err = sum_squares_err + square_err
      }
    }
  }
  return(sum_squares_err)
}

#4. weighted_kmeans function
weighted_kmeans = function(df, K) {
  #df is a dataframe with columns: desa name, longitude, latitude, luas_sawah_ha (as weight)
  #K is number of clusters desired
  #initialize containers
  list_cluster = list()
  vect_sse = vector()
  #iterate algorithm 3 times (due to local minima)
  for (iter in c(1:3)) {
    #initialization
    init_centroids_index = sample(nrow(df),K)
    distance_matrix = matrix(data = NA, nrow = nrow(df), ncol = K)
    cluster = vector()
    centroid_long = vector()
    centroid_lat = vector()
    #compute distance between desa and initial centroids (haversine distance)
    for (k in c(1:K)) {
      for (i in c(1:nrow(df))) {
        desa_i = as.numeric(df[i,2:3])
        centroid_k = as.numeric(df[init_centroids_index[k],2:3])
        distance_matrix[i,k] = haversine_dist(desa_i,centroid_k)
      }
    }
    #initial cluster assignment for each desa
    for (i in c(1:nrow(df))) {
      cluster[i] = which.min(distance_matrix[i,])
    }
    #iteration baseline
    old_cluster = vector(length = length(cluster))
    new_cluster = cluster
    #iterations
    while (!all(old_cluster == new_cluster)) {
      #update old cluster assignment
      old_cluster = new_cluster
      #calculate centroids (weighted average)
      for (k in c(1:K)) {
        cluster_k = which(old_cluster == k) #desa index of cluster k
        centroid_long[k] = weighted.mean(df$longitude[cluster_k], df$luas_sawah_ha[cluster_k])
        centroid_lat[k] = weighted.mean(df$latitude[cluster_k], df$luas_sawah_ha[cluster_k])
      }
      df_centroid = as.data.frame(cbind(centroid_long, centroid_lat))
      #compute distance between desa and centroids
      for (k in c(1:K)) {
        for (i in c(1:nrow(df))) {
          desa_i = as.numeric(df[i,2:3])
          centroid_k = as.numeric(df_centroid[k,])
          distance_matrix[i,k] = haversine_dist(desa_i,centroid_k)
        }
      }
      #update cluster assignment for each desa
      for (i in c(1:nrow(df))) {
        cluster[i] = which.min(distance_matrix[i,])
      }
      #update new_cluster
      new_cluster = cluster
    }
    #algor iterations result
    list_cluster[[iter]] = cluster #cluster assignment
    vect_sse[iter] = SSE(df, cluster) #sum of squares error
  }
  #function result
  best_iter = which.min(vect_sse)
  return(list_cluster[[best_iter]]) #best cluster assignment (smallest SSE)
}

#5. find the optimal K using elbow method
#import library
library(ggplot2)
#initiate container
vect_see = vector()
#fill vect_see via looping
K=5
for (K in c(1:8)) {
  cluster_assignment = weighted_kmeans(df_desa, K)
  vect_see[K] = SSE(df_desa, cluster_assignment)
}

#create data frame to use ggplot
df_elbow = as.data.frame(cbind(number_cluster = c(2:8), SSE = vect_see[-1]))

#plotting SSE decline over K (number of clusters)
ggplot(data = df_elbow, aes(x = number_cluster, y = SSE)) + geom_line(color = 'skyblue') + geom_point(color = 'red', size = 3)

#6. plotting clustering result and centroids
K = 5 #optimal K from elbow method
cluster_assignment = weighted_kmeans(df_desa, K)
df_desa_cluster = as.data.frame(cbind(longitude = df_desa$longitude, latitude = df_desa$latitude, cluster = cluster_assignment))

#calculate centroids (weighted average)
centroid_long = vector()
centroid_lat = vector()

for (k in c(1:K)) {
  cluster_k = which(cluster_assignment == k) #desa index of cluster k
  centroid_long[k] = weighted.mean(df_desa$longitude[cluster_k], df_desa$luas_sawah_ha[cluster_k])
  centroid_lat[k] = weighted.mean(df_desa$latitude[cluster_k], df_desa$luas_sawah_ha[cluster_k])
}

#create data frame for centroid with dummy cluster number 
df_centroid = as.data.frame(cbind(longitude = centroid_long, latitude = centroid_lat, cluster = rep(K+1, length(centroid_lat))))

#append df_desa_cluster and df_centroid for ggplot
df_kmeans_result = rbind.data.frame(df_desa_cluster, df_centroid)

#the plot clusters
plot(df_kmeans_result$longitude[which(df_kmeans_result$cluster==1)],
     df_kmeans_result$latitude[which(df_kmeans_result$cluster==1)],
     col="violetred3",pch=19,
     xlim = c(min(df_kmeans_result$longitude)-0.1,max(df_kmeans_result$longitude)+0.1), 
     ylim=c(-8.45,-7.6),xlab = "longitude", ylab = "latitude",main = "Lokasi Pusat Klaster RMU hasil Perhitungan Weighted K-means")
points(df_kmeans_result$longitude[which(df_kmeans_result$cluster==2)],df_kmeans_result$latitude[which(df_kmeans_result$cluster==2)],col="cadetblue",pch=19)
points(df_kmeans_result$longitude[which(df_kmeans_result$cluster==3)],df_kmeans_result$latitude[which(df_kmeans_result$cluster==3)],col="antiquewhite",pch=19)
points(df_kmeans_result$longitude[which(df_kmeans_result$cluster==4)],df_kmeans_result$latitude[which(df_kmeans_result$cluster==4)],col="cadetblue1",pch=19)
points(df_kmeans_result$longitude[which(df_kmeans_result$cluster==5)],df_kmeans_result$latitude[which(df_kmeans_result$cluster==5)],col="antiquewhite3",pch=19)

#centroids
points(df_kmeans_result$longitude[which(df_kmeans_result$cluster==6)],df_kmeans_result$latitude[which(df_kmeans_result$cluster==6)],col="black",pch=7)

