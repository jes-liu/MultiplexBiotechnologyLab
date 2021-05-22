rm(list=ls())

# Points----
num_points <- 400  # number of points for each type
num_type_points <- 25  # number of types of points

radius <- 1
space <- 2  # the gap between two circles

ppx <- 100  # points per x axis
ppy <- num_points*num_type_points/ppx  # points per y axis

x_dim <- ppx*2*radius + space*ppx
y_dim <- ppy*2*radius + space*ppy

xpoints <- integer()
ypoints <- integer()
for (i in 1:ppx) {
  coords <- seq(from=radius, to=y_dim, by=2*radius+space)
  xpoints <- append(xpoints, rep.int(coords[i], ppy), after=length(xpoints))
  ypoints <- append(ypoints, coords, after=length(ypoints))
}

pal <- colorRampPalette(c("red", "yellow"))
pal <- pal(num_type_points)
pal <- rep(pal, num_points)
rand_colors = sample(pal, ppx*ppy, replace=FALSE)

df_points <- data.frame(X = xpoints, Y = ypoints, Color = rand_colors)

plot(df_points$X, df_points$Y, xlim=c(0-1,x_dim), ylim=c(0-1,y_dim), 
     col=df_points$Color, xlab="X", ylab="Y", pch=19)

for (i in 1:(num_points*num_type_points)) {
  symbols(df_points$X[i], df_points$Y[i], circles=radius,
          add=TRUE, inches=FALSE, bg=df_points$Color[i])
}

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
  df$Distance = dist
  return(df)
}

df_points <- find_nearest_neighbors(df_points)
hist(df_points$Distance)
summary(df_points$Distance)
sd(df_points$Distance)

mean_by_color <- aggregate(df_points[, 4], list(df_points$Color), mean)
barplot(mean_by_color$x, col=mean_by_color$Group.1, xlab='Colors (by code)', ylab='mean')
sd(mean_by_color$x)

sorted_df <- df_points[order(-df_points$Distance),]
sorted_df[101,]$Distance
sorted_df[501,]$Distance
