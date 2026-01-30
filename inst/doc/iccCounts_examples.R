## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", fig.width=6, fig.height=4
)


## -----------------------------------------------------------------------------
library(iccCounts)

## -----------------------------------------------------------------------------
AF_P<-icc_counts(AF,y="y",id="id",met="met",type="con")

## -----------------------------------------------------------------------------
ICC(AF_P)

## -----------------------------------------------------------------------------
VarComp(AF_P)

## -----------------------------------------------------------------------------
set.seed(100)
AF_P.gof<-GOF_check(AF_P)

## ----fig.width=6, fig.height=4------------------------------------------------
plot(AF_P.gof,type="envelope")

## -----------------------------------------------------------------------------
plot(AF_P.gof,type="dispersion")
DispersionTest(AF_P.gof)

## -----------------------------------------------------------------------------
plot(AF_P.gof,type="zeros")
ZeroTest(AF_P.gof)

## -----------------------------------------------------------------------------
AF_NB1<-icc_counts(AF,y="y",id="id",met="met",type="con", fam="nbinom1")
AF_NB2<-icc_counts(AF,y="y",id="id",met="met",type="con", fam="nbinom2")
AIC(AF_NB1$model)
AIC(AF_NB2$model)

## -----------------------------------------------------------------------------
set.seed(100)
AF_NB2.gof<-GOF_check(AF_NB2)

## -----------------------------------------------------------------------------
plot(AF_NB2.gof,type="envelope")

## -----------------------------------------------------------------------------
plot(AF_NB2.gof,type="dispersion")
DispersionTest(AF_NB2.gof)

## -----------------------------------------------------------------------------
ICC(AF_NB2)

## -----------------------------------------------------------------------------
VarComp(AF_NB2)

## -----------------------------------------------------------------------------
G_P<-icc_counts(Grimso,y="Tot",id="TransectID")
ICC(G_P)
VarComp(G_P)

## ----fig.height=12, fig.width=8-----------------------------------------------
set.seed(100)
G_P.gof<-GOF_check(G_P)
plot(G_P.gof)

## -----------------------------------------------------------------------------
DispersionTest(G_P.gof)
ZeroTest(G_P.gof)


## -----------------------------------------------------------------------------
EPP_P<-icc_counts(EPP,y="Social",id="id")
ICC(EPP_P)
VarComp(EPP_P)

## -----------------------------------------------------------------------------
set.seed(100)
EPP_P.gof<-GOF_check(EPP_P)
plot(EPP_P.gof,type="envelope")

## -----------------------------------------------------------------------------
plot(EPP_P.gof,type="dispersion")
DispersionTest(EPP_P.gof)


## -----------------------------------------------------------------------------
plot(EPP_P.gof,type="zeros")
ZeroTest(EPP_P.gof)

## -----------------------------------------------------------------------------
EPP_NB1<-icc_counts(EPP,y="Social",id="id",fam="nbinom1")
EPP_NB2<-icc_counts(EPP,y="Social",id="id",fam="nbinom2")
AIC(EPP_NB1$model)
AIC(EPP_NB2$model)

## -----------------------------------------------------------------------------
set.seed(100)
EPP_NB1.gof<-GOF_check(EPP_NB1)
plot(EPP_NB1.gof,type="envelope")


## -----------------------------------------------------------------------------
plot(EPP_NB1.gof,type="dispersion")
DispersionTest(EPP_NB1.gof)


## -----------------------------------------------------------------------------
plot(EPP_NB1.gof,type="zeros")
ZeroTest(EPP_NB1.gof)

## -----------------------------------------------------------------------------
EPP_ZIP<-icc_counts(EPP,y="Social",id="id",fam="zip")
ICC(EPP_ZIP)
VarComp(EPP_ZIP)

## -----------------------------------------------------------------------------
set.seed(100)
EPP_ZIP.gof<-GOF_check(EPP_ZIP)
plot(EPP_ZIP.gof,type="envelope")

## -----------------------------------------------------------------------------
plot(EPP_ZIP.gof,type="dispersion")
DispersionTest(EPP_ZIP.gof)
plot(EPP_ZIP.gof,type="zeros")
ZeroTest(EPP_ZIP.gof)

