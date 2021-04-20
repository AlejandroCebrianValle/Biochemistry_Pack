# Extracted from Huitema Carly, Horsman Geoff; Analizing Kinetic Data Using the
# Powerful Ststistical Capabilities of R; bioRxiv 316588.

## Preparation to Script
rm(list = ls())

# Fitting the Michaellis-Menten equation
## Data
conc.uM <- c(0.5, 1, 2.5, 3.5, 5, 7.5, 10, 15, 25, 50, 70, 75, 100)
rate <- c(0.6, 1.1, 2.1, 2.3, 3.7, 3.0, 4.3, 4.8, 5.3, 6.0, 5.1, 5.7, 5.8)

### creation of data frame
expl.df <- data.frame(conc.uM, rate)
expl.df

### ploting data
plot(conc.uM, rate, 
     main = "Plot Title",
     xlab = "conc (uM)",
     ylab = "rate (uM/min)",
     lines(conc.uM, rate, lty = "dotted", col = "red"))
    # Upper line creates a line which joins all points with a red-dotted line

### Fit Michelis-Menten equation
mn.nls <- nls(rate ~ (Vmax * conc.uM / (Km + conc.uM)),
              data = expl.df,
              start = list(Km = 5, Vmax = 6))
summary(mn.nls)
#### Extraction coefficients
Km <- unname(coef(mn.nls) ["Km"])
Vmax <- unname(coef(mn.nls)["Vmax"])
#### Plotting data and line of best fit
x <- c(0:100)
y <- (Vmax * x / (Km + x))
lines(x,y, lty = "dotted", col = "blue")
#### Confidence intervals of parameters
confint(mn.nls)
#### Look at residuals and plot
mn.resid <- resid(mn.nls)
plot(expl.df$conc.uM, mn.resid)
#### add weighting to fit
expl.df$weight <- 1/expl.df$conc^2
mn.weight.nls <- nls(rate ~ (Vmax * conc.uM / (Km + conc.uM)),
                     data = expl.df,
                     start = list(Km = 5, Vmax = 6),
                     weight = expl.df$weight)
summary(mn.weight.nls)

# Determining IC_50
## Data
conc.uM <- c(300, 150, 75, 38, 19, 5, 2, 1, 0.6)
percent.activity <- c(2, 7, 12, 22, 36, 53, 67, 83, 85)
ic50.df <- data.frame(conc.uM, percent.activity)
ic50.df$conc.nM <- ic50.df$conc.uM * 1000
ic50.df$logconc.nM <- log10(ic50.df$conc.nM)
plot(ic50.df$logconc.nM, ic50.df$percent.activity)
### Estimation values of curve examining the plot
x <- c(2:12/2)
y <- 0 + (80 - 0)/(1 + (x / 4)^10)
lines(x, y, col="red")
### Fiting data using the nls()
rodbard.fit <- nls(formula(percent.activity ~
                               bot + (top - bot)/(1 +
                               (ic50.df$logconc.nM / logic50)^slope)),
                   data = ic50.df,
                   algorithm = "port",
                   start = list(bot = 0, top = 80, logic50 = 4, slope = 10),
                   lower = c(bot = -Inf, top = -Inf, logic50 = 0, slope = -Inf))
summary(rodbard.fit)
### extract coeficients
top <- unname(coef(rodbard.fit)["top"])
bot <- unname(coef(rodbard.fit)["bot"])
logic50 <- unname(coef(rodbard.fit)["logic50"])
slope <- unname(coef(rodbard.fit)["slope"])
### extract calculate line of best fit
y.fit <- bot + (top - bot)/(1 + (x / logic50)^slope)
lines(x, y.fit, col = "green")
### log scale to linear scale to get IC50 (results in nM)
10^logic50

