library(lme4)
lme(mAch ~ meanses*cses + sector*cses,
    data = Hsb82, random = ~ cses|school)
lme1(mAch ~ meanses*cses + sector*cses,
    data = Hsb82, random = ~ cses|school)
q("no")
