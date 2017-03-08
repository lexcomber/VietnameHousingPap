######### Lex Comber
######### a.comber@leeds.ac.uk
######### 08 March 2017
######### Please let me know if there any glitches! 
######### Enjoy!
######### for more fun see https://goo.gl/Lia9cI

# some packages 
# if you have not installed them before
# use install.packages("package.name", dep = T) 
library(GISTools)
library(foreign)
library(GWmodel)
library(car)
library(OpenStreetMap)
library(RgoogleMaps)
library(PBSmapping)
library(spgwr)
library(scales)
library(classInt)
library(tidyverse)
library(repmis)

##### Part 1: Data prep and initial investigation
# load data from GITHUB site
source_data("https://github.com/lexcomber/VietnameHousingPap/blob/master/Veitpap.RData?raw=True")

## 1. Exploratory Regression
# set up a regression model with ALL variables
terms <- names(df)[c(4:32)]
reg.mod <- paste(terms[1], "~")
for (i in 2:(length(terms)-1)) {
	reg.mod <- paste(reg.mod, terms[i], "+")
}
reg.mod <- paste(reg.mod, terms[length(terms)])
reg.mod <- as.formula(reg.mod)
# use for OLS regression
lm.1 <- lm(reg.mod, data = na.omit(df))
summary(lm.1)
## Stepwise to determine best set of variables
stepAIC(lm.1,direction = "both")
all.reg.mod <- as.formula(HPRICVND ~ HPRICPSM + GFA + PlotArea + DisCent + 
    TotNum + InnCY + HousePerm + HousePriv + ShopFYes + EdUni + 
    CarYes + StreetBus + WasteGood + SecGood)
# now select variables and clean for NAs
vars = gsub('[+]','',as.character(c(all.reg.mod)))
vars = gsub('  ',' ',vars)
vars <- unlist(strsplit(vars," "))
vars <- vars[-2]
index <- as.vector(apply(df[, vars], 1, function(x) sum(is.na(x))) == 0)
df.tmp <- df[index, vars]
df.tmp <- cbind(df[index, 1:3], df.tmp)

# overwrite data.sp with just selected variables
data.sp <- SpatialPointsDataFrame(df.tmp[,2:3], data = data.frame(df.tmp), 
	proj4string = crs.val)
# reproject data for distance based analysis such as GW models 
data.sp <- spTransform(data.sp, utm.proj)

# Map the data
library(OpenStreetMap)
fac = 0.01 # in case a bit is needed
tmp <- spTransform(data.sp, crs.val)
ul <- as.vector(cbind(bbox(tmp)[2,2]+fac, bbox(tmp)[1,1]-fac))
lr <- as.vector(cbind(bbox(tmp)[2,1]-fac, bbox(tmp)[1,2]+fac))
MyMap <- openmap(ul,lr,zoom = NULL)
#### Figure 1
# print to PDF
# setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
# pdf(file = "f1.pdf")
plot(MyMap, removeMargin=T) 
plot(spTransform(data.sp, osm()),add = T, col="#08519C80", pch = 19)
# dev.off()
##### Part 2: Initial analyses
## Exploratory Regression
# set up a regression model with SELECTED variables
terms <- names(data.sp)[4:18]
reg.mod <- paste(terms[1], "~")
for (i in 2:(length(terms)-1)) {
	reg.mod <- paste(reg.mod, terms[i], "+")
}
reg.mod <- paste(reg.mod, terms[length(terms)])
reg.mod <- as.formula(reg.mod)
# use for OLS regression
lm.1 <- lm(reg.mod, data = data.sp@data)
summary(lm.1)
tab <- summary(lm.1)[4]
tab <- round(tab[[1]], 3)
##### Table 1
## Global VIF and CN
# VIFs > 10, VDP > 0.5 suggestive of collinearity
vif(lm.1) # all low
tab <- cbind(tab, append(NA, vif(lm.1)))
colnames(tab)[5] <- "VIF"
tab <- round(tab, 3)
# setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
# write.csv(tab, "Table1.csv")
##### Condition number of the design matrix
X <- as.matrix(cbind(1,data.sp@data[, 4:18]))
BKWcn <- function(X) {
	p <- dim(X)[2]
	Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
	Xsvd <- svd(Xscale)$d
	Xsvd[1] / Xsvd[p]
}
BKWcn(X)
# The BKW condition number (Belsley et al. 1980), is found to be 27.4
# indicative that collinearity is at least, a global problem for this data. 

##### Part 3: Exploratory GWR collinearity 
# GWR euclidean and netowrk distance model calibration
bw.ed <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV") 
bw.nd <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV", dMat = dMat.ndp2) 
# GW collinearity diagnostic
gw.col <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.nd, 
  kernel = "bisquare", adaptive = TRUE, dMat = dMat.ndp2)
summary(gw.col$SDF[,15])
##### Figure 2
# now plot the data - note the use of the spTransform function
cols <- add.alpha((brewer.pal(5, "Greys")), 0.7)
index.CN30 <-  gw.col$SDF@data[,15] > 30
sum(index.CN30)
cex.i <- rep(0.7, length(index.CN30))
cex.i[index.CN30] <- 2.4
# setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
# pdf(file = "f2.pdf")
plot(MyMap, removeMargin=T) 
plot(spTransform(data.sp, osm()),add = T, col = cols[4], pch = 16, cex = cex.i)
# dev.off()

##### Part 4: Prediction Comparison and evalution
# ED-GWR ND-GWR ED-LCR ND-LCR 
# ED-GWR
gwr.mod.ed <- gwr.robust(reg.mod, data = data.sp, 
	kernel = "bisquare", bw = bw.ed, adaptive = TRUE) 
ed.gwr <- gwr.mod.ed$SDF
# ND-GWR
gwr.mod.nd <- gwr.robust(reg.mod, data = data.sp, kernel = "bisquare", 
	bw = bw.nd, adaptive = TRUE, dMat = dMat.ndp2) 
nd.gwr <- gwr.mod.nd$SDF
# ED-LCR  
lcr.bw.ed <- bw.gwr.lcr(reg.mod, data = data.sp, kernel = "bisquare", 
	adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
lcr.mod.ed <- gwr.lcr(reg.mod, data = data.sp, bw = lcr.bw.ed, kernel = "bisquare", 
	adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
ed.lcr <- lcr.mod.ed$SDF
# ND-LCR 
lcr.bw.nd <- bw.gwr.lcr(reg.mod, data = data.sp, kernel = "bisquare", 
	adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30, dMat = dMat.ndp2)
lcr.mod.nd <- gwr.lcr(reg.mod, data = data.sp, bw = lcr.bw.nd, kernel = "bisquare",
	adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30, dMat = dMat.ndp2)
nd.lcr <- lcr.mod.nd$SDF

# A global regression - GWR with a massive bandwidth = OLS
gwr.mod <- gwr.basic(reg.mod, data = data.sp, kernel = "bisquare", 
  bw = 10000000000000000000000, adaptive = FALSE) 
# the RSS scores
gwr.mod$GW.diagnostic$RSS
gwr.mod.ed$GW.diagnostic$RSS
gwr.mod.nd$GW.diagnostic$RSS
lcr.mod.ed$GW.diagnostic$RSS
lcr.mod.nd$GW.diagnostic$RSS


# extract data
vars = gsub('[+]','',as.character(reg.mod))
vars = gsub('  ',' ',vars)
vars <- unlist(strsplit(vars," "))[-c(1,2)]
sp.i <- data.sp[,c(vars)]
# now predict using function
predict.func <- function(sp.i, gwr.mod) {
  gwr.pred <- gwr.mod[, 2:15] 
  gwr.pred.intersect <- gwr.mod[, 1]
  res <- (sp.i@data * gwr.pred@data) 
  res <- cbind(gwr.pred.intersect@data, res)
  res <- rowSums(res)
  res <- cbind(res, data.sp@data[, c("HPRICVND")])
  colnames(res) <- c("res", "HPRICVND")
  return(res)
}
ed.gwr.res <- predict.func(sp.i, ed.gwr)
nd.gwr.res <- predict.func(sp.i, nd.gwr)
ed.lcr.res <- predict.func(sp.i, ed.lcr)
nd.lcr.res <- predict.func(sp.i, nd.lcr)

# map function with ggplot
ggplot.func.lex <- function(res) {
  # GET EQUATION AND R-SQUARED AS STRING
  # SOURCE: http://goo.gl/K4yh
  # from http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
  lm_eqn <- function(df){
    y <- df[,1]
    x <- df[,2]
    m <- lm(y ~ x, df);
    a <- sum(m$residuals^2)
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 3), 
                          b = format(coef(m)[2], digits = 3), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    return(list(as.character(as.expression(eq)), a)) }
    
    res <- data.frame(res)
    rss1 <- round(lm_eqn(res[, c(2, 1)])[[2]], 0)
    p <- ggplot(data = res[!index.CN30,], aes(x = res, y = HPRICVND)) +
    geom_point(color = "#25252580") +
    scale_y_continuous(name = "Survey Price") + 
    scale_x_continuous(name = "Predicted Price") +
    geom_point(data = res[index.CN30,], aes(x = res, y = HPRICVND), color= "#DE2D26CC") +
    geom_text(x = range(res[,1])[2]%/% 2, y = 50000, label = paste("RSS:", rss1), 
      parse = TRUE, cex = 6) +
    geom_text(x = range(res[,1])[2]%/% 2, y = 55000, label = lm_eqn(res[, c(2, 1)])[[1]], 
      parse = TRUE, cex = 6) +   
    geom_smooth(method = lm) +
    theme_minimal() +
    theme_bw()
  return(p)
}
#AIC(lm(res[, c(2,1)]))
p1 <- ggplot.func.lex(ed.gwr.res)
p2 <- ggplot.func.lex(nd.gwr.res)
p3 <- ggplot.func.lex(ed.lcr.res)
p4 <- ggplot.func.lex(nd.lcr.res)
# plot 
# you may wish to write these PDF etc
# eg pdf(file = "f3d.PDF"); p4; dev.off()
setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
pdf(file = "f3_a.PDF"); p1; dev.off()
pdf(file = "f3_b.PDF"); p2; dev.off()
pdf(file = "f3_c.PDF"); p3; dev.off()
pdf(file = "f3_d.PDF"); p4; dev.off()

## 5. FINAL LCR-GW  
# Now gwr.lcr to adjust for the ridge term
lcr.bw <- bw.gwr.lcr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
lcr.bw
# this will find the bw to fit a ridge term in every location 
# where the CN was greater than 30 
lcr.mod <- gwr.lcr(reg.mod, data = data.sp, regression.points = gwr.grid, bw = lcr.bw, 
  kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30, dMat = dMat.ndp2g2)

##### Figure 4: plot coeffienct estimates
lcr.res <- lcr.mod$SDF
tab <- t(round(apply(lcr.mod$SDF@data[, 2:15], 2, fivenum), 1))
colnames(tab) <- c("Min.", "1st Quartile", "Median", "3rd Quartile", "Max.")
tab[, c(2,4)]
setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/RPlay/")
write.csv(tab, "table2.csv")

var.names <- c("Intercept", "Price per m2", "Ground floor area (m2)", "Plot area (m2)", 
              "Travel time to centre (mins)",  "No. in household", "Inner Courtyard (T)", 
              "Permanent house (T)", "Private house (T)", "Shop front (T)", 
              "University education (T)", "Car ownership (T)", 
              "Business / commercial street (T)", "Good waste collection (T)", 
              "Good security (T)")
length(var.names)
# setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
# pdf("f4.pdf", w = 12)
par(mfrow = c(2,5))
for (i in c(5, 7:15)) {
	v <- v <- as.vector(lcr.res@data[,i])
	if (median(v) > 0) {
		cols <- add.alpha("#CD5C5C", 0.5)
		#cols = "indianred"
		plot(MyMap, removeMargin=T) 
		plot(spTransform(gwr.grid, osm()),add = T, col=cols[1], 
			cex = (v/max(v))-0.15, pch = 23, bg = cols[1])
		title(var.names[i])
	}
	if (median(v) < 0) {
		cols <- add.alpha(brewer.pal(5, "Blues")[4], 0.5)
		plot(MyMap, removeMargin=T) 
		plot(spTransform(gwr.grid, osm()),add = T, col=cols[1], 
			cex = (v/min(v))-0.15, pch = 23, bg = cols[1])
		title(var.names[i])
	}
}
#dev.off()

## END






