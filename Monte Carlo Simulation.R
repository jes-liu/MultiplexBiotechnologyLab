rm(list=ls())

# Points----
num_points <- 40  # number of points for each type
num_type_points <- 25  # number of types of points

x_dim <- 50
y_dim <- 20

xpoints <- integer()
ypoints <- integer()
for (i in 1:x_dim) {
  xpoints <- append(xpoints, rep.int(i, y_dim), after=length(xpoints))
  ypoints <- append(ypoints, 1:y_dim, after=length(ypoints))
}

pal <- colorRampPalette(c("red", "yellow"))
pal <- pal(num_type_points)
rand_colors = sample(pal, x_dim*y_dim, replace=TRUE)

df_points <- data.frame(X = xpoints, Y = ypoints, Color = rand_colors)

plot(df_points$X, df_points$Y, xlim=c(1,x_dim), ylim=c(1,y_dim), col=df_points$Color,
     xlab="X", ylab="Y")

# Nearest Distance then plot distribution----
find_nearest_neighbors <- function(df) {
  colors <- unique(df$Color)
  
  dist <- integer()
  
  for (i in 1:length(colors)) {  # sort by color
    df_color <- df[df$Color == colors[i],]
    
    for (j in 1:nrow(df_color)) {  # going down each point one by one
      x0 <- df_color[j, 1]
      y0 <- df_color[j, 2]
      distances <- integer()
      df_color_rem <- df_color[-j,]

      for (k in 1:nrow(df_color_rem)) {  # finding distance for each point from initial point
        x <- df_color_rem[k, 1]
        y <- df_color_rem[k, 2]
        z = sqrt((x-x0)^2 + (y-y0)^2)
        distances <- append(distances, z)
      }
      min_dist = min(distances)
      dist <- append(dist, min_dist, after=as.numeric(rownames(df_color)[j])-1)
    }
  }
  return(dist)
}

distance <- find_nearest_neighbors(df_points)
hist(distance)
summary(distance)
