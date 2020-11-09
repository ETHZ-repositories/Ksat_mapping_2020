Global prediction of soil saturated hydraulic conductivity using random
forest in a Covariate-based Geo Transfer Functions (CoGTF) framework
================
Surya Gupta, Tom Hengl, Peter Lehmann, Sara Bonetti, Dani Or

    library(sp)
    library(rgdal)
    library(raster)
    library(foreign)
    library(hydroGOF)
    
    
    tt2<-read.csv("E:/Ksat_dataset_mapping.csv")
    
    ## Unique IDs for 5 degrees by 5 degrees
    
    unique(tt2$FID_Fish_n)
    
    
    source("E:/OpenLandMap/R/saveRDS_functions.R")
    source("E:/OpenLandMap/R/LandGIS_functions.R")
    
    ## saveRDS_functions.R and LandGIS_functions.R available at https://github.com/Envirometrix/LandGISmaps/tree/477460d1d0099646c508f65e68769b9edf050ce8/functions
    
    ## 3D modeling (see Hengl, T., & MacMillan, R. A. (2019). Predictive soil mapping with R. Lulu. com.)
    
    dfs <- hor2xyd(tt2, U="hzn_top", L="hzn_bot")
    
    
    I.vars = make.names(unique(unlist(sapply(c("s.no","clm_", "dtm_", "lcv", "veg_", "olm_c", "olm_s", "olm_bd", "FID_Fish_", "DEPTH"), function(i){names(dfs)[grep(i, names(dfs))]}))))
    
    t.vars = c("log_ksat")
    sel.n <- c(t.vars,I.vars)
    sel.r <- complete.cases(dfs[,sel.n])
    PTF_temp2 <- dfs[sel.r,sel.n]
    
    ##Set1
    set.seed(16)
    chosen <- sample(unique(PTF_temp2$FID_Fish_n), 30)
    
    ff<-subset(PTF_temp2, FID_Fish_n %in% chosen)
    
    final<-PTF_temp2[!(PTF_temp2$FID_Fish_n %in% ff$FID_Fish_n),]
    
    set.seed(6)
    chosen <- sample(unique(final$FID_Fish_n), 30)
    
    ff1<-subset(final, FID_Fish_n %in% chosen)
    
    final1<-final[!(final$FID_Fish_n %in% ff1$FID_Fish_n),]
    
    
    set.seed(11)
    chosen <- sample(unique(final1$FID_Fish_n), 30)
    
    ff2<-subset(final1, FID_Fish_n %in% chosen)
    
    final2<-final1[!(final1$FID_Fish_n %in% ff2$FID_Fish_n),]
    
    
    set.seed(14)
    chosen <- sample(unique(final2$FID_Fish_n), 30)
    
    ff3<-subset(final2, FID_Fish_n %in% chosen)
    
    final3<-final2[!(final2$FID_Fish_n %in% ff3$FID_Fish_n),]
    
    ##Set2
    
    set.seed(33)
    chosen <- sample(unique(PTF_temp2$FID_Fish_n), 30)
    
    ff<-subset(PTF_temp2, FID_Fish_n %in% chosen)
    
    final<-PTF_temp2[!(PTF_temp2$FID_Fish_n %in% ff$FID_Fish_n),]
    
    set.seed(39)
    chosen <- sample(unique(final$FID_Fish_n), 30)
    
    ff1<-subset(final, FID_Fish_n %in% chosen)
    
    final1<-final[!(final$FID_Fish_n %in% ff1$FID_Fish_n),]
    
    
    set.seed(57)
    chosen <- sample(unique(final1$FID_Fish_n), 30)
    
    ff2<-subset(final1, FID_Fish_n %in% chosen)
    
    final2<-final1[!(final1$FID_Fish_n %in% ff2$FID_Fish_n),]
    
    
    set.seed(67)
    chosen <- sample(unique(final2$FID_Fish_n), 30)
    
    ff3<-subset(final2, FID_Fish_n %in% chosen)
    
    final3<-final2[!(final2$FID_Fish_n %in% ff3$FID_Fish_n),]
    
    ##Set3
    
    set.seed(77)
    chosen <- sample(unique(PTF_temp2$FID_Fish_n), 30)
    
    ff<-subset(PTF_temp2, FID_Fish_n %in% chosen)
    
    final<-PTF_temp2[!(PTF_temp2$FID_Fish_n %in% ff$FID_Fish_n),]
    
    set.seed(85)
    chosen <- sample(unique(final$FID_Fish_n), 30)
    
    ff1<-subset(final, FID_Fish_n %in% chosen)
    
    final1<-final[!(final$FID_Fish_n %in% ff1$FID_Fish_n),]
    
    
    set.seed(96)
    chosen <- sample(unique(final1$FID_Fish_n), 30)
    
    ff2<-subset(final1, FID_Fish_n %in% chosen)
    
    final2<-final1[!(final1$FID_Fish_n %in% ff2$FID_Fish_n),]
    
    
    set.seed(115)
    chosen <- sample(unique(final2$FID_Fish_n), 30)
    
    ff3<-subset(final2, FID_Fish_n %in% chosen)
    
    final3<-final2[!(final2$FID_Fish_n %in% ff3$FID_Fish_n),]
    
    
    
    df1<-ff
    df2<-ff1
    df3<-ff2
    df4<-ff3
    df5<-final3
    
    
    Train1<- rbind(ff, ff1, ff2, ff3)
    
    Train2<- rbind (ff1, ff2, ff3,final3)
    
    Train3<- rbind(ff2, ff3,final3, ff)
    
    Train4<- rbind(ff3,final3, ff,ff1)
    
    Train5<- rbind(final3, ff,ff1, ff2)
    
    
    grid <- list.files("E:/maps_tests/new_layers/layers_RS/" , pattern = "*.tif$")
    All_cov <- raster::stack(paste0("E:/maps_tests/new_layers/layers_RS/", grid))
    
    set.seed(2) 
    fm.ksat <- as.formula(paste("log_ksat~ ",paste(names(All_cov), collapse = "+")))
    fm.ksat
    
    set.seed(2) 
    rm.ksat <- Train1[complete.cases(Train1[,all.vars(fm.ksat)]),]
    m.ksat <- ranger(fm.ksat, rm.ksat, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat
    
    df5$prediction<- predict(m.ksat,df5)$predictions
    
    ## Ist_part is computed
    
    rm.ksat1 <- Train2[complete.cases(Train2[,all.vars(fm.ksat)]),]
    m.ksat1 <- ranger(fm.ksat, rm.ksat1, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat1
    
    df1$prediction<- predict(m.ksat1,df1)$predictions
    
    ## 2nd_part is computed
    rm.ksat2 <- Train3[complete.cases(Train3[,all.vars(fm.ksat)]),]
    m.ksat2 <- ranger(fm.ksat, rm.ksat2, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat2
    
    df2$prediction<- predict(m.ksat2,df2)$predictions
    
    ## 3rd_part is computed
    
    rm.ksat3 <- Train4[complete.cases(Train4[,all.vars(fm.ksat)]),]
    m.ksat3 <- ranger(fm.ksat, rm.ksat3, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat3
    df3$prediction<- predict(m.ksat3,df3)$predictions
    
    ## 4th_part is computed
    rm.ksat4 <- Train5[complete.cases(Train5[,all.vars(fm.ksat)]),]
    m.ksat4 <- ranger(fm.ksat, rm.ksat4, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat4
    
    df4$prediction<- predict(m.ksat4,df4)$predictions
    
    Final_data<- rbind(df1,df2,df3,df4,df5)
    
    dd<- aggregate(Final_data[, 1:33], list(Final_data$s_no), mean)
    
    
    
    ll<- lm(dd$prediction~ dd$log_ksat)
    
    RMSE(dd$prediction,dd$log_ksat)
    
    summary(ll)
    
    
    ccc = DescTools::CCC(dd$prediction,dd$log_ksat, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
    ccc
    
    bias(dd$prediction,dd$log_ksat)
    
    colnames(dd)
    
    dd$log_ksat1<- 10^dd$log_ksat
    dd$prediction1<- 10^dd$prediction
    
    hexbinplot(log_ksat1~ prediction1, 
               panel = function(x, y, ...){
                 panel.hexbinplot(x, y, ...)
                 panel.loess(x, y,span = 2/3, col.line = "blue",type="l", lty=2, lwd = 4)
                 panel.abline(c(0, 1),lwd = 2)
               },
               data = dd,xlab = "Predicted Ksat [cm/day]", ylab = "Measured Ksat [cm/day]",cex.axis = 4, aspect="1", xbins=30, colramp = function(n) {viridis (8,  alpha = 1, begin = 0, end = 1, direction = -1,option = "C")},xlim=c(0.00000000001,10000000), ylim=c(0.00000000001,10000000),
               scales=list(
                 x = list(log = 10, equispaced.log = FALSE), 
                 y = list(log = 10, equispaced.log = FALSE)
               ),
               font.lab= 6, cex.labels = 1.2,font.axis = 2,colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1) )
    
    
    ##Fitting final model
    
    rm.ksat <- PTF_temp2[complete.cases(PTF_temp2[,all.vars(fm.ksat)]),]
    m.ksat <- ranger(fm.ksat, rm.ksat, num.trees=200, mtry=6, quantreg = TRUE)
    m.ksat
    
    
    ##Importance variable
    m.ksat <- randomForest(fm.ksat, rm.ksat, num.trees=200, mtry=6, quantreg = TRUE)
    
    varImpPlot(m.ksat, sort=TRUE, n.var=min(30, nrow(m.ksat$importance)))
    
    ## Then Produced the final CoGTF map
    
    #p2 = predict(All_cov,m.ksat, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)
    
    #writeRaster(p2, "/home/step/data/OpenLandMap/Final_selected_covariates/New_raster_layer/Final_RF_0cm.tif")
