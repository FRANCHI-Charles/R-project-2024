data = read.csv("./Data_Cortex_Nuclear.csv")

### remove Mice ID
data = data[,-1]

summary(data)

### Transform class as factor
classes = list("Genotype", "Behavior", "Treatment", "class")

for (text in classes){
  data[,text] = factor(data[,text])
}

### Transform protein class as numeric

for (i in 1:77){
  data[,i] = as.numeric(data[,i])
}

### Handling Missing values

summary(data[1:43]) # the 43 first column have max 18 missing values

missinglist = c()
for (i in 1:43){
  missinglist = union(missinglist, which(is.na(data[i])))
}
cat("Number of instances with missing values for the 43 first proteins :", length((missinglist)))

# Remove this missing values
data = data[-missinglist,]

# Other columns : replace by the mean
missingColumns <- colnames(data)[apply(data, 2, anyNA)]

for (text in missingColumns){
  data[is.na(data[text]), text] = mean(data[,text], na.rm = TRUE)
}


