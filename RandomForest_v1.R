
library(tidyverse)
library(randomForest)
#install.packages("randomForest")
#install.packages("caret")
library(caret)



##datasets
g1.meta <- read_csv("G1_MetaData.csv")
g1.metab <- read_csv("G1_Metab_Data.csv")
peridice <- read_csv("PERIDICE_metabolite_data_v2.csv")


#Wrangle Data

#Make RF friendly metabolite names (remove spaces, dashes, and numbers from start of metab names so that random forest function works)
metab.names <- peridice %>%
  select(metabolite) %>%
  unique() %>%
  mutate(metab.short = str_remove(metabolite, " ")) %>%
  mutate(metab.short = str_remove_all(metab.short, "-")) %>%
  mutate(metab.short = str_remove_all(metab.short, "[[:digit:]]")) %>%
  mutate(metab.short = str_remove_all(metab.short, "\\(")) %>%
  mutate(metab.short = str_remove_all(metab.short, "\\)")) %>%
  mutate(metab.short = str_remove_all(metab.short, ","))


#Wrangle PERI-DICE Data
peridice.metab <- peridice %>%
  select(filename, metabolite, nmol_per_pc) %>%
  group_by(filename, metabolite) %>%
  summarize(nmol_per_pc = mean(nmol_per_pc)) %>%    ###summarize with mean to handle duplicated metabolite values with different PC values
  ungroup() %>%
  unique() %>%
  left_join(., metab.names) %>%
  select(-metabolite) %>%
  pivot_wider(id_cols = filename, names_from = metab.short, values_from = nmol_per_pc)

peridice.N <- peridice %>%
  select(filename, added_N_uM) %>%
  unique() %>%
  filter(!is.na(added_N_uM)) 

peridice.dat <- left_join(peridice.N, peridice.metab) %>%
  select(-filename) %>%
  rename("N_uM" = added_N_uM) 
  
peridice.dat[is.na(peridice.dat)] <- 0




#Wrangle Gradients 1 Data
g1.metab.long <- g1.metab %>%
  select(-Compound_name_in_figures) %>%
  rename("metabolite" = Complete_compound_name) %>%
  pivot_longer(cols = S11C1_A:U9_B, names_to = "Sample_ID", values_to = "nmol") %>%
  left_join(., g1.meta) %>%
  mutate(nmol_per_pc = nmol/(PC_nM/1000))

g1.dat <- g1.metab.long %>%
  select(Sample_ID, metabolite, nmol_per_pc, NO3_NO2) %>%
  replace_na(list(nmol = 0)) %>%
  left_join(., metab.names) %>%
  select(-metabolite) %>%
  filter(!is.na(metab.short)) %>%
  pivot_wider(names_from = metab.short, values_from = nmol_per_pc) %>%
  rename("N_uM" = NO3_NO2) %>%
  select(-Sample_ID)


g1.dat[is.na(g1.dat)] <- 0



###########Run Random Forest Models 


###_________________PERI-DICE___________________

##split data set into test and train
set.seed(222)
peri.ind <- sample(2, nrow(peridice.dat), replace = TRUE, prob = c(0.7, 0.3))

peri.train <- peridice.dat[peri.ind==1,]

peri.test<- peridice.dat[peri.ind==2,]

##train RF 
peri.train.rf <- randomForest(N_uM ~ ., data = peri.train, ntree = 500, mtry = 16)
peri.train.rf

##test RF
peri.predict <- predict(peri.train.rf, peri.test)
peri.predict

#calculate mean sequared error
peri.predict.mse <- mean((peri.predict-peri.test$N_uM)^2)
peri.predict.mse



## train on full data set
peri.full.rf <- randomForest(N_uM ~ ., data = peridice.dat, ntree = 500)
peri.full.rf

#compare error to # of trees
plot(peri.full.rf)

#tune model using caret:
peri.control <- trainControl(method="repeatedcv", number=5, repeats=3, search="grid")
peri.tunegrid <- expand.grid(.mtry=c(1:40))
metric <- "RMSE"

##Runing grid search takes a long time so only do if you have to:
#peri.rf_gridsearch <- train(N_uM~., data=peridice.dat, method="rf", metric=metric, tuneGrid=peri.tunegrid, trControl=peri.control)
print(peri.rf_gridsearch)
plot(peri.rf_gridsearch)


##rerun model using tuned mtry parameters
peri.rf.tuned <- randomForest(N_uM ~ ., data = peridice.dat, ntree = 500, mtry = 16)
peri.rf.tuned


#pull out variable importance 
peri.imp <- as.data.frame(importance(peri.rf.tuned)) 

###___Plot Peridice Full Tuned Model Results
peri.pred <- as.data.frame(predict(peri.rf.tuned)) %>%
  rename("Predicted_N" = `predict(peri.rf.tuned)`)

peri.comparison <- peridice.dat %>%
  select(N_uM) %>%
  rename("Actual_N" = N_uM) %>%
  cbind(., peri.pred)

ggplot(peri.comparison, aes(x = Actual_N, y = Predicted_N)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 6) +
  ylim(0, 6) +
  xlab("Measured N (uM)") +
  ylab("Predicted N (uM)")



###___________Gradients 1_______________

##split data set into test and train
g1.ind <- sample(2, nrow(g1.dat), replace = TRUE, prob = c(0.7, 0.3))

g1.train <- g1.dat[g1.ind==1,]

g1.test<- g1.dat[g1.ind==2,]

##train RF 
g1.train.rf <- randomForest(N_uM ~ ., data = g1.train, ntree = 500, mtry = 32)
g1.train.rf

##test RF
g1.predict <- predict(g1.train.rf, g1.test)
g1.predict

#calculate mean sequared error
g1.predict.mse <- mean((g1.predict-g1.test$N_uM)^2)
g1.predict.mse


##full data set
g1.full.rf <- randomForest(N_uM ~ ., data = g1.dat, ntree = 500)
g1.full.rf


#tune model using caret:
g1.control <- trainControl(method="repeatedcv", number=5, repeats=3, search="grid")
g1.tunegrid <- expand.grid(.mtry=c(1:40))
metric <- "RMSE"

g1.rf_gridsearch <- train(N_uM~., data=g1.dat, method="rf", metric=metric, tuneGrid=g1.tunegrid, trControl=g1.control)
print(g1.rf_gridsearch)
plot(g1.rf_gridsearch)

##rerun model using tuned mtry parameters
g1.rf.tuned <- randomForest(N_uM ~ ., data = g1.dat, ntree = 500, mtry = 32)
g1.rf.tuned

#pull out variable importance 
g1.imp <- as.data.frame(importance(g1.rf.tuned)) 



###___Plot Gradients 1 Full Tuned Model Results
g1.pred <- as.data.frame(predict(g1.rf.tuned)) %>%
  rename("Predicted_N" = `predict(g1.rf.tuned)`)

g1.comparison <- g1.dat %>%
  select(N_uM) %>%
  rename("Actual_N" = N_uM) %>%
  cbind(., g1.pred)

ggplot(g1.comparison, aes(x = Actual_N, y = Predicted_N)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 6) +
  ylim(0, 6) +
  xlab("Measured N (uM)") +
  ylab("Predicted N (uM)")









####_____Test PERI-DICE and Gradients 1 on eachother to evaluate performance:

##Test PERI RF on G1 data:
xtest.peri <- predict(peri.rf.tuned, g1.dat)
xtest.peri

#calculate mean sequared error
xtest.peri.mse <- mean((xtest.peri-g1.dat$N_uM)^2)
xtest.peri.mse

###Plot Results:
xtest.peri.pred <- as.data.frame(xtest.peri) %>%
  rename("Predicted_N" = `xtest.peri`)

xtest.peri.comparison <- g1.dat %>%
  select(N_uM) %>%
  rename("Actual_N" = N_uM) %>%
  cbind(., xtest.peri.pred)

ggplot(xtest.peri.comparison, aes(x = Actual_N, y = Predicted_N)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 6) +
  ylim(0, 6) +
  xlab("Measured N (uM)") +
  ylab("Predicted N (uM)")



##_____________________________________________
##Test G1 RF on PERI data:
xtest.g1 <- predict(g1.rf.tuned, peridice.dat)
xtest.g1

#calculate mean sequared error
xtest.g1.mse <- mean((xtest.g1-peridice.dat$N_uM)^2)
xtest.g1.mse

###Plot Results:
xtest.g1.pred <- as.data.frame(xtest.g1) %>%
  rename("Predicted_N" = `xtest.g1`)

xtest.g1.comparison <- peridice.dat %>%
  select(N_uM) %>%
  rename("Actual_N" = N_uM) %>%
  cbind(., xtest.g1.pred)

ggplot(xtest.g1.comparison, aes(x = Actual_N, y = Predicted_N)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 6) +
  ylim(0, 6) +
  xlab("Measured N (uM)") +
  ylab("Predicted N (uM)")



####____Plot results of Variable Importance

#Combine into one dataset:
peri.var.imp <- peri.imp %>%
  rownames_to_column(var = "metab.short") %>%
  rename("PERI.import" = IncNodePurity) %>%
  left_join(., metab.names)

g1.var.imp <- g1.imp %>%
  rownames_to_column(var = "metab.short") %>%
  rename("g1.import" = IncNodePurity) %>%
  left_join(., metab.names)

full.var.imp <- left_join(peri.var.imp, g1.var.imp)

#scatter plot
ggplot(full.var.imp, aes(x = g1.import, y = PERI.import)) +
  geom_point()


#Trying a different approach: 
peri.var.imp.2 <- peri.imp %>%
  rownames_to_column(var = "metab.short") %>%
  rename("import" = IncNodePurity) %>%
  left_join(., metab.names) %>%
  mutate(dataset = "PERI")

g1.var.imp.2 <- g1.imp %>%
  rownames_to_column(var = "metab.short") %>%
  rename("import" = IncNodePurity) %>%
  left_join(., metab.names) %>%
  mutate(dataset = "G1")

full.var.imp.2 <- full_join(peri.var.imp.2, g1.var.imp.2) %>%
  group_by(dataset) %>%
  arrange(desc(import)) %>%
  slice_max(import, n = 10)
  

ggplot(full.var.imp.2, aes(x = import, y = reorder(metabolite, import))) +
  geom_col() +
  facet_grid(.~dataset, scales = "free") +
  ylab("10 most important metbolites for each RFR model") +
  xlab("metabolite importance") +
  theme_bw()


###Plotting NO3_NO2
ggplot(g1.meta, aes(x = Binned_latitude, y = NO3_NO2)) +
  geom_point(size = 3, alpha = 0.8, color = "darkblue") +
  geom_line(color = "darkblue") +
  theme_bw() +
  xlab("Latitude (Degrees N)") +
  ylab("Inorganic nitrogen concentration (uM)")

















