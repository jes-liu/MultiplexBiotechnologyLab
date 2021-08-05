rm(list=ls())
# what is r? ----
# r is a programming language mainly used for statistical analysis and graphing

# what is r-studio? ----
# when you download r, it only downloads the console where all the function
# are executed. console on bottom
# r-studio, however, is an IDE (integrated development environment) to be able
# to create scripts and longer functions for r

# naming objects ----
# assignment convention & executing that line
a <- 2
a <- 1 + 1
b <- 5*a

# case-sensitive (a =/= A) ----
A <- 3

# variable names ----
12obj <- 5
_obj <- 6
(obj) <- 9
obj.1 <- 1
obj_2 <- 2

# data classes ----
84.9  # numeric
"char"  # character
T/F/TRUE/FALSE  #boolean

# data structures ----
c()  # combines inputted values into vector or list
?c 
x <- c(1, 2, 3, 5, 10)
y <- c(1:10)
z <- c("red", 12, TRUE, a)
x + 2
y*2
y^y

# vector ----
# 1 column/row of data of 1 data type
vector <- c("red", "blue", "green")
vector

# list ----
# 1 column/row of data of >1 data type
l <- list("red", 12, TRUE, a)
l
num <- list(1, 2, 3, 4)
num*2
is.vector()
is.list()
list(vector, a)
c(list, 36)
c(vector,36)

# subsetting ----
vector <- c("red", "blue", "green", "yellow", "purple", "grey")
vector[1]
vector[1:4]
vector[c(1, 3, 6)]
vector[-5]

# matrix ----
# 2+ columns/rows of data of 1 data type
mat <- matrix(1:9, nrow=3, ncol=3)
mat[1, 2]

# dataframe ----
# 2+ columns/rows of data of >1 data type
df <- data.frame("Age" = c(18:20),
                 "Gender" = c("Male", "Female", "Male"))
df$Age
df$Gender
df[,1]
df[,1:2]
df[1,]
df[1:2,]
df$Eyes <- c("Blue", "Brown", "Green")
df





