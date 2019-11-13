## R Basics for the lab meeting

getwd() #pwd equivalent 
setwd('~/../../Volumes/Ainslie_USB/VibData/PreProcessedData')

# Doing basic maths
1+1
easyqn <- 3+4 #this way is prefered as it always means you are delaring something to the workspace
easyqn = 4-10
worse <- median(x=1:10)
better <- median(x<-1:10)
easymistake <- isTRUE(easyqn <- 1)
easyqn = 4-10
corrected <- isTRUE(easyqn < -1)
hardqn <- sqrt(abs((3*sin(4))-1))


# Variables
my_age <- 27 # Numeric variable
my_name <- 'Ainslie' # Character variable
# Am I a cat lover?: (yes/no) <=> (TRUE/FALSE)
is_catlover <- FALSE # logical variable


# Vectors
fav_numbs <- c(17, 9, 42, 3)
mean(fav_numbs)
max(fav_numbs) 
atluae <- fav_numbs[3]


# Matricies
other_numbs <- c(1,2,3,4)
m1 <- rbind(fav_numbs, other_numbs)
is.matrix(m1)
m2 <- t(m1)
m3 <- matrix(data=c(fav_numbs, other_numbs),  nrow=4, ncol=2, byrow=TRUE)
colnames(m3) <- c("fav_numbs","other_numbs")


# Data Frames
groups <- (c("grp1", "grp2", "grp1", "grp2"))
df <- data.frame(groups,m3)
df$groups <- factor(df$groups)
tapply(df[,2], df$groups, mean)



