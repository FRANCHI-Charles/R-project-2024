### Packages needed
library(reshape2)
library(ggplot2) #for heatmap
library(factoextra) #for PCA visualization
library(MASS) #for LDA
library(caTools) #for train test split
library(e1071) #for svm

set.seed(689)

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
p = ggplot(melt(cormap), aes(Var1, Var2, fill = value))+  # plot this map
  geom_tile(color = "grey")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", # use bilateral colors
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()
print(p)

# remove redundant features
todel = c()
for (i in 1:dim(highcor)[1]) {
  if (highcor[i, 1] > highcor[i,2]) { # to not do [2, 1] and [1, 2] for example
    todel = union(todel, highcor[i,1])
  }
}
dataX = dataX[-todel]


### Visualisation
data.pca = prcomp(dataX)
summary(data.pca)
for (i in 1:(length(dataclass)-1)){
  x11()
  p = fviz_pca_biplot(data.pca, geom.ind= c("point"), col.var = "black", col.ind = rainbow(2)[dataclass[,i]], select.var = list(cos2 = 0.3),
                  title = colnames(dataclass)[i]) +guides(shape="none", fill="none") +
    scale_color_manual(name = "Class", labels = unique(dataclass[i])[,1], values = rainbow(2))
  print(p)
  #plot(pca$x, main = colnames(dataclass)[i], pch = 21, bg = rainbow(2)[dataclass[,i]])
  #legend("bottomleft", legend = unique(dataclass[i])[,1], pch = 21, pt.bg=rainbow(2), cex = 1.3)
}

x11()
plot(data.pca$x, main = colnames(dataclass)[4], pch = 21, bg = rainbow(8)[dataclass[,4]])
legend("topleft", legend = unique(dataclass[4])[,1], pch = 21, pt.bg=rainbow(8), cex = 1)


data.lda = lda(dataX, grouping = dataclass[,4])
projection = as.matrix(dataX)%*%as.matrix(cbind(data.lda$scaling[,1],data.lda$scaling[,2]))
x11()
plot(projection, pch = 21, bg = rainbow(8)[dataclass[,4]], main= "LDA")
legend("top", legend = unique(dataclass[4])[,1], pch = 21, pt.bg=rainbow(8), cex = 1)


### SVM

# train test split
split = sample.split(dataX, SplitRatio = 0.85)
trainX = subset(dataX, split == TRUE)
trainy = subset(dataclass, split == TRUE)
summary(trainy) # check if we don"t have a too bad distribution of labels

testX = subset(dataX, split == FALSE)
testy = subset(dataclass, split == FALSE)
summary(testy)  # check if we don"t have a too bad distribution of labels

# all categories
clf = svm(trainX, y=trainy[,4], kernel="linear", type = 'C-classification')

predy = predict(clf, testX)
confmatrix = table(testy[,4], predy)
cat("\nAll categories\n")
print(confmatrix)

# Memantine
clf = svm(trainX, y=trainy[,2], kernel="linear", type = 'C-classification')

predy = predict(clf, testX)
confmatrix = table(testy[,2], predy)
cat("\nTreatment\n")
print(confmatrix)
cat("Important features : ", names(sort(abs(coef(clf)), decreasing = TRUE)[1:10]))


