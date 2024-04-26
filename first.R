### Packages needed
library(reshape2)
library(ggplot2)

data = read.csv("./Data_Cortex_Nuclear.csv")
data = data[sample(1:nrow(data)),] # shuffle the data

### remove Mice ID
data = data[,-1]

summary(data)
cat("With 'summary(data)', we can see the 33 first features do not contain too much NA's.
We can also see the class are balanced.")

### Transform class as factor
classes = list("Genotype", "Behavior", "Treatment", "class")

for (text in classes){
  data[,text] = factor(data[,text])
}

### Transform protein class as numeric

for (i in 1:77){
  data[,i] = as.numeric(data[,i])
}

### Check for redundancies

cat("\nNumber of redudant rows :", length(data) - length(unique(data)))


### Handling Missing values

summary(data[1:43]) # the 43 first column have max 18 missing values

missinglist = c()
for (i in 1:43){
  missinglist = union(missinglist, which(is.na(data[i])))
}
cat("\nNumber of instances with missing values for the 43 first proteins :", length((missinglist)))

# Remove this missing values
data = data[-missinglist,]

# Other columns : replace by the mean
missingColumns <- colnames(data)[apply(data, 2, anyNA)]

cat("\nOther category with missing values :", missingColumns, "\n")

for (text in missingColumns){
  data[is.na(data[text]), text] = mean(data[,text], na.rm = TRUE)
}

### One class at a time
dataX = data[1:77]
dataclass = data[78:81]


### Data Analysis

#extract high correlated features to show only them
threshold = 0.75

highcor = which(abs(cor(dataX)) >= threshold, arr.ind = TRUE) # element with abs correlation >0.8
highcor = highcor[! highcor[,1] == highcor[,2],] # remove diagonal elements
cormap = cor(data[unique(highcor[,1])]) # correlation of the high correlated elements only
cormap[abs(cormap) < threshold] = NA # disable display of the 

x11()
ggplot(melt(cormap), aes(Var1, Var2, fill = value))+  # plot this map
  geom_tile(color = "grey")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", # use bilateral colors
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()

#x11()
#pairs(dataBehavior[unique(highcor[,1])], pch = 21, bg = rainbow(2)[dataBehavior[,78]])

pca = prcomp(dataX)
summary(pca)
for (i in 1:(length(dataclass)-1)){
  x11()
  plot(pca$x, main = colnames(dataclass)[i], pch = 21, bg = rainbow(2)[dataclass[,i]])
}
biplot(pca)

