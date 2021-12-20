rm(list=ls())

library(openxlsx)

dir = 'C:/Users/jesse/OneDrive/Documents/MATLAB/data.xlsx'

raw = data.frame(read.xlsx(dir, sheet=1))
raw = transform(raw, intensity = as.numeric(intensity),
                x.position = as.numeric(x.position), 
                y.position = as.numeric(y.position))
raw = raw[1:(nrow(raw)-1),]

hist(raw$intensity, main='Intensity frequency distribution', xlab='Intensities')


bright_threshold = 1300
n_beads = nrow(raw)

# take control's intensity
control = mean(raw[raw$intensity>1000 & raw$intensity<1100,]$intensity, na.rm=T)


# take surrounding's intensity
all_distance = integer()
all_surround = integer()

# not inclusive*
distance_1300_1400 = integer()
surround_1300_1400 = integer()

distance_1400_1500 = integer()
surround_1400_1500 = integer()

distance_1500_1600 = integer()
surround_1500_1600 = integer()

distance_1600_1700 = integer()
surround_1600_1700 = integer()

distance_1700_1800 = integer()
surround_1700_1800 = integer()

distance_1800_1900 = integer()
surround_1800_1900 = integer()

distance_1900_2000 = integer()
surround_1900_2000 = integer()

distance_2000_2500 = integer()
surround_2000_2500 = integer()

distance_2500_3000 = integer()
surround_2500_3000 = integer()

distance_3000_4000 = integer()
surround_3000_4000 = integer()

for (x in 1:n_beads) {
  print(x)
  if (raw$intensity[x] > bright_threshold) {
    all_distance = append(all_distance, 0)
    all_surround = append(all_surround, raw$intensity[x])
    
    
    if (raw$intensity[x]<1400) {
      distance_1300_1400 = append(distance_1300_1400, 0)
      surround_1300_1400 = append(surround_1300_1400, raw$intensity[x])
    } 
    else if (raw$intensity[x]<1500) {
      distance_1400_1500 = append(distance_1400_1500, 0)
      surround_1400_1500 = append(surround_1400_1500, raw$intensity[x])
    } 
    else if (raw$intensity[x]<1600) {
      distance_1500_1600 = append(distance_1500_1600, 0)
      surround_1500_1600 = append(surround_1500_1600, raw$intensity[x])
    } 
    else if (raw$intensity[x]<1700) {
      distance_1600_1700 = append(distance_1600_1700, 0)
      surround_1600_1700 = append(surround_1600_1700, raw$intensity[x])
    } 
    else if (raw$intensity[x]<1800) {
      distance_1700_1800 = append(distance_1700_1800, 0)
      surround_1700_1800 = append(surround_1700_1800, raw$intensity[x])
    } 
    else if (raw$intensity[x]<1900) {
      distance_1800_1900 = append(distance_1800_1900, 0)
      surround_1800_1900 = append(surround_1800_1900, raw$intensity[x])
    } 
    else if (raw$intensity[x]<2000) {
      distance_1900_2000 = append(distance_1900_2000, 0)
      surround_1900_2000 = append(surround_1900_2000, raw$intensity[x])
    } 
    else if (raw$intensity[x]<2500) {
      distance_2000_2500 = append(distance_2000_2500, 0)
      surround_2000_2500 = append(surround_2000_2500, raw$intensity[x])
    } 
    else if (raw$intensity[x]<3000) {
      distance_2500_3000 = append(distance_2500_3000, 0)
      surround_2500_3000 = append(surround_2500_3000, raw$intensity[x])
    } 
    else if (raw$intensity[x]<4000) {
      distance_3000_4000 = append(distance_3000_4000, 0)
      surround_3000_4000 = append(surround_3000_4000, raw$intensity[x])
    } 
    
    
    x0 = raw$x.position[x]
    y0 = raw$y.position[x]
    
    index = which(raw$x.position<=x0+8 & raw$x.position>=x0-8 & raw$y.position<=y0+8 & raw$y.position>=y0-8)
    # which returns index of TRUE rows
    
    index = index[index!=x]
    
    # find the distance and see if the bead qualifies to be within the range of the spillover
    for (y in 1:length(index)) {
      x1 = raw$x.position[index[y]]
      y1 = raw$y.position[index[y]]
      dist = sqrt((x1-x0)^2 + (y1-y0)^2)
      
      # to put the "distance" in the correct format to make parabola
      if (x1-x0 < 0) {
        dist = -dist
      }
      
      if (abs(dist) <= 8 & abs(dist) >= 2.5) {

        if (raw$intensity[x]<1400) {
          distance_1300_1400 = append(distance_1300_1400, dist)
          surround_1300_1400 = append(surround_1300_1400, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<1500) {
          distance_1400_1500 = append(distance_1400_1500, dist)
          surround_1400_1500 = append(surround_1400_1500, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<1600) {
          distance_1500_1600 = append(distance_1500_1600, dist)
          surround_1500_1600 = append(surround_1500_1600, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<1700) {
          distance_1600_1700 = append(distance_1600_1700, dist)
          surround_1600_1700 = append(surround_1600_1700, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<1800) {
          distance_1700_1800 = append(distance_1700_1800, dist)
          surround_1700_1800 = append(surround_1700_1800, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<1900) {
          distance_1800_1900 = append(distance_1800_1900, dist)
          surround_1800_1900 = append(surround_1800_1900, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<2000) {
          distance_1900_2000 = append(distance_1900_2000, dist)
          surround_1900_2000 = append(surround_1900_2000, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<2500) {
          distance_2000_2500 = append(distance_2000_2500, dist)
          surround_2000_2500 = append(surround_2000_2500, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<3000) {
          distance_2500_3000 = append(distance_2500_3000, dist)
          surround_2500_3000 = append(surround_2500_3000, raw$intensity[index[y]])
        } 
        else if (raw$intensity[x]<4000) {
          distance_3000_4000 = append(distance_3000_4000, dist)
          surround_3000_4000 = append(surround_3000_4000, raw$intensity[index[y]])
        } 
        
        all_distance = append(all_distance, dist)
        all_surround = append(all_surround, raw$intensity[index[y]])
      }
    }
  }
}

plot(distance_1300_1400, surround_1300_1400)

df_1300_1400 = data.frame(distance=distance_1300_1400, surround=surround_1300_1400)
df_1400_1500 = data.frame(distance=distance_1400_1500, surround=surround_1400_1500)
df_1500_1600 = data.frame(distance=distance_1500_1600, surround=surround_1500_1600)
df_1600_1700 = data.frame(distance=distance_1600_1700, surround=surround_1600_1700)
df_1700_1800 = data.frame(distance=distance_1700_1800, surround=surround_1700_1800)
df_1800_1900 = data.frame(distance=distance_1800_1900, surround=surround_1800_1900)
df_1900_2000 = data.frame(distance=distance_1900_2000, surround=surround_1900_2000)
df_2000_2500 = data.frame(distance=distance_2000_2500, surround=surround_2000_2500)
df_2500_3000 = data.frame(distance=distance_2500_3000, surround=surround_2500_3000)
df_3000_4000 = data.frame(distance=distance_3000_4000, surround=surround_3000_4000)


#######################
dir_1300_1400 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1300_1400.xlsx'
df_1300_1400_r = data.frame(x=df_1300_1400[df_1300_1400$surround<bright_threshold,]$distance, 
                     y=df_1300_1400[df_1300_1400$surround<bright_threshold,]$surround)
write.xlsx(df_1300_1400_r, dir_1300_1400, overwrite=TRUE)

dir_1400_1500 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1400_1500.xlsx'
df_1400_1500_r = data.frame(x=df_1400_1500[df_1400_1500$surround<bright_threshold,]$distance, 
                     y=df_1400_1500[df_1400_1500$surround<bright_threshold,]$surround)
write.xlsx(df_1400_1500_r, dir_1400_1500, overwrite=TRUE)

dir_1500_1600 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1500_1600.xlsx'
df_1500_1600_r = data.frame(x=df_1500_1600[df_1500_1600$surround<bright_threshold,]$distance, 
                     y=df_1500_1600[df_1500_1600$surround<bright_threshold,]$surround)
write.xlsx(df_1500_1600_r, dir_1500_1600, overwrite=TRUE)

dir_1600_1700 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1600_1700.xlsx'
df_1600_1700_r = data.frame(x=df_1600_1700[df_1600_1700$surround<bright_threshold,]$distance, 
                     y=df_1600_1700[df_1600_1700$surround<bright_threshold,]$surround)
write.xlsx(df_1600_1700_r, dir_1600_1700, overwrite=TRUE)

dir_1700_1800 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1700_1800.xlsx'
df_1700_1800_r = data.frame(x=df_1700_1800[df_1700_1800$surround<bright_threshold,]$distance, 
                     y=df_1700_1800[df_1700_1800$surround<bright_threshold,]$surround)
write.xlsx(df_1700_1800_r, dir_1700_1800, overwrite=TRUE)

dir_1800_1900 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1800_1900.xlsx'
df_1800_1900_r = data.frame(x=df_1800_1900[df_1800_1900$surround<bright_threshold,]$distance, 
                     y=df_1800_1900[df_1800_1900$surround<bright_threshold,]$surround)
write.xlsx(df_1800_1900_r, dir_1800_1900, overwrite=TRUE)

dir_1900_2000 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_1900_2000.xlsx'
df_1900_2000_r = data.frame(x=df_1900_2000[df_1900_2000$surround<bright_threshold,]$distance, 
                     y=df_1900_2000[df_1900_2000$surround<bright_threshold,]$surround)
write.xlsx(df_1900_2000_r, dir_1900_2000, overwrite=TRUE)

dir_2000_2500 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_2000_2500.xlsx'
df_2000_2500_r = data.frame(x=df_2000_2500[df_2000_2500$surround<bright_threshold,]$distance, 
                     y=df_2000_2500[df_2000_2500$surround<bright_threshold,]$surround)
write.xlsx(df_2000_2500_r, dir_2000_2500, overwrite=TRUE)

dir_2500_3000 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_2500_3000.xlsx'
df_2500_3000_r = data.frame(x=df_2500_3000[df_2500_3000$surround<bright_threshold,]$distance, 
                     y=df_2500_3000[df_2500_3000$surround<bright_threshold,]$surround)
write.xlsx(df_2500_3000_r, dir_2500_3000, overwrite=TRUE)

dir_3000_4000 = 'C:/Users/jesse/OneDrive/Documents/MATLAB/spillover/spillover_3000_4000.xlsx'
df_3000_4000_r = data.frame(x=df_3000_4000[df_3000_4000$surround<bright_threshold,]$distance, 
                     y=df_3000_4000[df_3000_4000$surround<bright_threshold,]$surround)
write.xlsx(df_3000_4000_r, dir_3000_4000, overwrite=TRUE)
#######################



#make density plot to show the number of points and their intensities per pixel distance




# subtract control from equation, then use eq to find spillover per bead using distance
eq_dir = 'C:/Users/jesse/OneDrive/Documents/MATLAB/parabolic.xlsx'
eq = data.frame(read.xlsx(eq_dir, colNames = F))


list_of_df = list(df_1300_1400, df_1400_1500, df_1500_1600, df_1600_1700, df_1700_1800,
                  df_1800_1900, df_1900_2000, df_2000_2500, df_2500_3000, df_3000_4000)
n_eq = length(list_of_df)
empty = list()

for (x in 1:n_eq) {
  n_surround = nrow(data.frame((list_of_df[x])))
  sub = integer()
  print(x)
  for (y in 1:n_surround) {
    d = data.frame(list_of_df[x])$distance[y]  # x=distance, y=intensity
    intensity = data.frame(list_of_df[x])$surround[y]
    spillover = eq[x,1]*abs(d)^4 + eq[x,2]*abs(d)^3 + eq[x,3]*abs(d)^2 + eq[x,4]*abs(d) + eq[x,5] - control
    fin = intensity - spillover
    sub = append(sub, fin)
  }
  empty = append(empty, list(data.frame(sub)))
}
full = empty
rm(empty)

######### make final data frame and then plot
fin_1300_1400 = data.frame(distance = data.frame(list_of_df[1])$distance, fin = data.frame(full[1]))
fin_1400_1500 = data.frame(distance = data.frame(list_of_df[2])$distance, fin = data.frame(full[2]))
fin_1500_1600 = data.frame(distance = data.frame(list_of_df[3])$distance, fin = data.frame(full[3]))
fin_1600_1700 = data.frame(distance = data.frame(list_of_df[4])$distance, fin = data.frame(full[4]))
fin_1700_1800 = data.frame(distance = data.frame(list_of_df[5])$distance, fin = data.frame(full[5]))
fin_1800_1900 = data.frame(distance = data.frame(list_of_df[6])$distance, fin = data.frame(full[6]))
fin_1900_2000 = data.frame(distance = data.frame(list_of_df[7])$distance, fin = data.frame(full[7]))
fin_2000_2500 = data.frame(distance = data.frame(list_of_df[8])$distance, fin = data.frame(full[8]))
fin_2500_3000 = data.frame(distance = data.frame(list_of_df[9])$distance, fin = data.frame(full[9]))
fin_3000_4000 = data.frame(distance = data.frame(list_of_df[10])$distance, fin = data.frame(full[10]))


plot(fin_1300_1400$distance, fin_1300_1400$sub)

