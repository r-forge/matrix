library(lme4)
lme(pixel ~ day + I(day^2), Pixel, list(Dog = ~ day, DS = ~ 1))
show(lme1(pixel ~ day + I(day^2), Pixel, list(Dog = ~ day, DS = ~ 1)))
q("no")
