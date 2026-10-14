# Quantile regression-based polygenic scores (QR-PGS) and prediction intervals
# Train quantile-specific PGS models and construct prediction intervals
# Reference (TBA)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(fsQRPPA)
library(RcppParallel)
options(scipen = 999)

phe = "Field_30710" # phenotype identifier
adj = 1 # inclusion of covariate predictors: 1 => covariate-adjusted; 0 => unadjusted
tau = 0.5 # quantile level: e.g., 0.025, 0.1, ..., 0.5, ..., 0.9, 0.975
nmz = "norm"
phe.nmz = paste0("Trait.", str_to_title(nmz))
nthread = 8


## 1. Phenotype and genotype data
dir = paste0("./example")
pfile = paste0(dir, "/example.phenotypes.tsv")
pheno = fread(file = pfile,
              header = T, sep = "\t", data.table = F, stringsAsFactors = F)
gfile = paste0(dir, "/example.genotypes.tsv")
geno = fread(file = gfile,
             header = T, sep = "\t", data.table = F, stringsAsFactors = F)
pheno = geno %>% 
  select(IID) %>% 
  left_join(y = pheno, by = join_by(IID))
cname = pheno %>% select(SEX:PC10) %>% colnames()
if (adj == 1) {
  geno = cbind(pheno %>% select(any_of(cname)), geno)  
}

## 2. QR-PGS training
#### 2.1 Model parameters
mu = 1e-5
c_eta <- 0.001
tol_abs <- 1e-4
tol_rel = 1e-3
lambda_min_ratio <- 0.01 
lambda_max_ratio <- 1.1
n_lambda = 9
step_LLA <- 1 
max_iter <- 32768
min_iter <- 16
n.group = min(ncol(geno), nthread)
RcppParallel::setThreadOptions(numThreads = nthread)

pen = "ELASTIC_NET"
a_param <- ifelse(pen == "SCAD", 3.7, 3)
lambda.pool = c(0) # lambda = 0 => ridge-penalized quantile regression
zeta.pool = c(1e-2, 1e-4, 1e-6, 0)
if (pen != "ELASTIC_NET") {
  zeta.pool = 0
}
param.grid = expand.grid(lambda = lambda.pool,
                         zeta = zeta.pool) %>% 
  mutate(warmup = lambda + zeta) %>%
  arrange(desc(warmup), desc(lambda), desc(zeta)) %>%
  filter(warmup != 0)

#### 2.2 Training
fit = fsQRPPA(X = as.matrix(geno[pheno$Set == "tr", ]),
              y = (pheno %>% filter(Set == "tr") %>% pull(as.symbol(phe.nmz))),
              tau = tau,
              lambda = param.grid$lambda,
              zeta = param.grid$zeta,
              penalty = pen,
              a = a_param,
              c_eta = c_eta,
              mu = mu,
              max_iter = max_iter,
              min_iter = min_iter,
              tol_abs = tol_abs,
              tol_rel = tol_rel,
              step_LLA = step_LLA,
              G = n.group,
              standardized = F,
              verbose = F)
  
loss.pinball = pinball_loss_path(X = as.matrix(geno[pheno$Set == "va", ]),
                                 y = (pheno %>% filter(Set == "va") %>% pull(as.symbol(phe.nmz))),
                                 tau = tau,
                                 G = n.group,
                                 Beta = fit$beta,
                                 standardized = F)

opt.idx = which.min(loss.pinball)
opt.lambda = param.grid$lambda[opt.idx]
opt.zeta = param.grid$zeta[opt.idx]
opt.beta = (fit$beta)[, opt.idx]
opt.iter = (fit$n_iter)[opt.idx]
opt.pred = as.matrix(geno %>% mutate(Ones = 1, .before = 1)) %*% as.matrix(opt.beta)
df.beta = data.frame(Tau = tau, 
                     Predictor = c("Intercept", colnames(geno)), Beta = opt.beta,  
                     N.Predictor = length(opt.beta) - 1, N.Feat = sum(opt.beta != 0) - 1, 
                     Lambda = opt.lambda, Zeta = opt.zeta, N.Iter = opt.iter)
df.pred = as.data.frame(opt.pred) %>%
  rename(PGS = V1) %>%
  mutate(Tau = tau) %>%
  mutate(pheno %>% select(IID, Set))

## 3. Model and predictions
#### 3.1 QR-PGS model coefficients at quantile level tau
ofile = paste0(dir, "/qr.beta.", nmz, ".adj", adj, ".tau", tau, ".tsv.gz")
fwrite(x = df.beta,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T, compress = "gzip")

#### 3.2 QR-PGS predictions at quantile level tau for all samples
ofile = paste0(dir, "/qr.pred.", nmz, ".adj", adj, ".tau", tau, ".tsv.gz")
fwrite(x = df.pred,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T, compress = "gzip")


