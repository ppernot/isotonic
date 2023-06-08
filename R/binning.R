figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

simul = FALSE # Read saved results of MC simulations

# BUS2022 ####
D = read.table(
  '../Data/BUS2022/qm9_U0_test_Orig.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
S  = D[, "formula"]
R  = D[, "U0"]
C  = D[, "prediction"]
uC = D[, "uncertainty_total"] ^ 0.5
E  = R - C
M = length(E)
aux = 1:M
uE = uC

# Uncertainty stratification from isotonic regression
dist = sort(table(uE),decreasing = TRUE)
print(length(dist))
print(sum(dist==1))

png(
  file = file.path(figDir, paste0('Fig01.png')),
  width  = 2*gPars$reso,
  height = gPars$reso
)
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))
par(
  mfrow = c(1, 2),
  mar = mar,
  mgp = mgp,
  pty = 'm',
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  cex.main = 1
)
x = as.numeric(names(dist))
y = as.numeric(dist)
plot(x, y, type = 'h',
     col = cols[5], lwd = 2*lwd,
     log = 'xy',
     xlab = 'Uncertainty [eV]',
     ylab = 'Frequency', yaxs = 'i', ylim =c(1,1700)
     )
grid(equilogs = FALSE)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)
plot(
  as.vector(dist), log = 'y',
  type = 'h', col = cols[5], lwd = 2*lwd,
  xlab = 'Order #',  xaxs = 'i', xlim = c(0,55),
  ylab = 'Frequency', yaxs = 'i', ylim =c(1,1700)
)
grid(equilogs = FALSE)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)
dev.off()

nBin = 50
png(
  file = file.path(figDir, paste0('Fig02a.png')),
  width  = gPars$reso,
  height = gPars$reso
)
res = ErrViewLib::plotRelDiag(
  uE, E,
  title = 'Data unordered',
  nBin = nBin,
  logX = TRUE,
  xlim = c(0.001, 0.2),
  ylim = c(0.001, 0.2),
  score = TRUE,
  method = 'cho',
  plot = TRUE,
  label = 1,
  gPars = gPars)
text(0.005,0.15, paste0('ENCE = ',signif(res$ENCE,2)))
dev.off()

png(
  file = file.path(figDir, paste0('Fig02b.png')),
  width  = gPars$reso,
  height = gPars$reso
)
res = ErrViewLib::plotRelDiag(
  uE, E, aux =abs(E),
  title = 'Data ordered by |E|',
  nBin = nBin,
  logX = TRUE,
  xlim = c(0.001, 0.2),
  ylim = c(0.001, 0.2),
  score = TRUE,
  method = 'cho',
  plot = TRUE,
  label = 2,
  gPars = gPars)
text(0.005,0.15, paste0('ENCE = ',signif(res$ENCE,2)))
dev.off()

# Dependence of ENCE & ZVE stats on aux ####

nBins = c(2,5,seq(10,160,by=10))
sel = M / nBins > 30
nBins = nBins[sel]

nMC = 250

if(simul) {
  e1 = e2 = z1 = z2 = matrix(NA,ncol = length(nBins), nrow = nMC)
  n1 = n2 = c()
  for(j in 1:nMC) {
    print(j)
    aux = sample(M,M)

    for(i in seq_along(nBins)) {
      nBin = nBins[i]

      rese1 = ErrViewLib::plotRelDiag(
        uE, E, aux = aux,
        nBin = nBin,
        score = TRUE,
        method = 'cho',
        plot = FALSE)

      resz1 = ErrViewLib::plotLZV(
        uE, E / uE, aux = aux,
        nBin = nBin,
        score = TRUE,
        plot = FALSE,
        method = 'cho')

      n1[i] = length(resz1$pc)
      e1[j,i] = rese1$ENCE
      z1[j,i] = resz1$ZVE
    }
  }

  # Original order
  aux = 1:M
  n1_0 = e1_0 = z1_0 = c()
  for(i in seq_along(nBins)) {
    nBin = nBins[i]

    rese1 = ErrViewLib::plotRelDiag(
      uE, E, aux = aux,
      nBin = nBin,
      score = TRUE,
      method = 'cho',
      plot = FALSE)

    resz1 = ErrViewLib::plotLZV(
      uE, E / uE, aux = aux,
      nBin = nBin,
      score = TRUE,
      plot = FALSE,
      method = 'cho')

    n1_0[i] = length(resz1$pc)
    e1_0[i] = rese1$ENCE
    z1_0[i] = resz1$ZVE
  }

  save(n1,e1,z1,n1_0,e1_0,z1_0,
       file='tmp.Rd')

} else {

  load(file='tmp.Rd')
}

png(
  file = file.path(figDir, paste0('Fig03.png')),
  width  = 2*gPars$reso,
  height = gPars$reso
)
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))
par(
  mfrow = c(1, 2),
  mar = mar,
  mgp = mgp,
  pty = 'm',
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  cex.main = 1
)
wf = 10 # Width factor for error bars
xlim = 6
x = sqrt(n1)
sel = x > xlim
xs = x[sel]
matplot(
  x, t(e1),
  pch = 16,
  col = gPars$cols_tr[5],
  xlim = c(0, max(x)),
  xlab = expression(N^{1/2}),
  xaxs = 'i',
  ylim = c(-0.01, 0.1),
  ylab = 'ENCE'
)
grid()
abline(v = xlim, lty = 3)
abline(h = 0, lty = 2)
valide = c()
for (j in 1:nMC) {
  zs = e1[j, sel]
  reg = lm(zs ~ xs)
  abline(reg = reg, col = gPars$cols_tr2[5])
  i1 = summary(reg)$coefficients[1,1]
  Ui1 = 2*summary(reg)$coefficients[1,2]
  y = i1; Uy = Ui1
  segments(0,y-Uy,0,y+Uy, lwd=wf*lwd, col = gPars$cols[5])
  valide[j] = (i1-Ui1)*(i1+Ui1) < 0
}
print(mean(valide))
points(x,e1_0,pch = 16, col = gPars$cols[2])
zs = e1_0[sel]
reg = lm(zs ~ xs)
abline(reg = reg, col = gPars$cols[2], lwd=6)
i1 = summary(reg)$coefficients[1,1]
Ui1 = 2*summary(reg)$coefficients[1,2]
y = i1; Uy = Ui1
segments(0,y-Uy,0,y+Uy, lwd=wf*lwd, col = gPars$cols[2])
box()

x = sqrt(n1)
sel = x > 6
xs = x[sel]
matplot(
  x, t(z1),
  pch = 16,
  col = gPars$cols_tr[5],
  xlim = c(0, max(x)),
  xlab = expression(N^{1/2}),
  xaxs = 'i',
  ylim = c(0.96, 1.25),
  ylab = 'ZVE'
)
grid()
abline(v = xlim, lty = 3)
abline(h = 1, lty = 2)
validz = c()
for (j in 1:nMC) {
  zs = z1[j, sel]
  reg = lm(zs ~ xs)
  abline(reg = reg, col = gPars$cols_tr2[5])
  i1 = summary(reg)$coefficients[1,1]
  Ui1 = 2*summary(reg)$coefficients[1,2]
  y = i1; Uy = Ui1
  segments(0,y-Uy,0,y+Uy, lwd=wf*lwd, col = gPars$cols[5])
  validz[j] = (i1-Ui1 - 1)*(i1+Ui1 - 1) < 0
}
print(mean(validz))
points(x,z1_0,pch = 16, col = gPars$cols[2])
zs = z1_0[sel]
reg = lm(zs ~ xs)
abline(reg = reg, col = gPars$cols[2], lwd=6)
i1 = summary(reg)$coefficients[1,1]
Ui1 = 2*summary(reg)$coefficients[1,2]
y = i1; Uy = Ui1
segments(0,y-Uy,0,y+Uy, lwd=wf*lwd, col = gPars$cols[2])
box()

dev.off()

