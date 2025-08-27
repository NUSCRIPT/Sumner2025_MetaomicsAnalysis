
source("src/common/load_metadata.r")
library(proxy)
library(infotheo)

data_discrete <- permmeta %>%
    mutate(across(everything(), as.character))


# Get column names
features <- colnames(data_discrete)

# Initialize distance matrix
n <- length(features)
mi_matrix <- matrix(0, n, n)
colnames(mi_matrix) <- features
rownames(mi_matrix) <- features

# Compute mutual information between features
for (i in 1:n) {
    for (j in 1:n) {
        mi_matrix[i, j] <- mutinformation(data_discrete[[i]], data_discrete[[j]])
    }
}

# Convert MI matrix to a distance matrix (distance = max(MI) - MI)
mi_max <- max(mi_matrix)
dist_matrix <- mi_max - mi_matrix

# Convert to dist object
dist_object <- as.dist(dist_matrix)


hc <- hclust(dist_object, method = "average")  # You can use "complete", "single", etc.
plot(hc, main = "Hierarchical Clustering with Mutual Information Distance")
