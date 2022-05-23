# 1. Load Packages -----------------------------------------------------------
require(jtools) #provides the effect_plot and cat_plot functions, which allow you to plot the partial residuals of your models

# 2. Models ---------------------------------------------------------------
NameOfYourPlotWithInteractions<-cat_plot(NameOfYourModel,pred=YourPredictorFromInteraction,modx=YourOtherPredictorFromInteraction,x.label="NAMEX", y.label="NAMEY
",interval =TRUE,int.type = "confidence", partial.residuals = TRUE, data = ______)

NameOfYourPlotWithOutInteractions<-effect_plot(NameOfYourModel, pred=YourPredictorOfInterest, x.label="Experimental Phase", y.label="Weighted 
Eigenvector Centrality", interval =TRUE, int.type = "confidence", partial.residuals = TRUE, data = ________)