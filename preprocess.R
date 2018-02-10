library(io);
library(dplyr);
library(rstan);

x <- qread("lung-bm-incidence-meta.tsv", stringsAsFactors=FALSE);

# set all blank values to NA
x[x == ""] <- NA;

d <- select(x, -other_primary_cancers, -other_metastases, -comments);

message("Number of studies: ", length(unique(d$pmid)))
message("Number of populations: ", nrow(d))

# exclude cohorts that specifically enrich or deplete BM
d <- filter(d,
	!grepl("B?M+", cancer_type),
	!grepl("BM- at dx", cancer_type),
	!grepl("neurologic", cancer_type)
);

message("Number of studies: ", length(unique(d$pmid)))
message("Number of populations: ", nrow(d))

# focus on studies counting all BM occurrences
d <- filter(d,
	n_bm_pos > 0
);

message("Number of studies: ", length(unique(d$pmid)))
message("Number of populations: ", nrow(d))

####

logit <- function(x) {
	log(x) - log(1 - x)
}

logistic <- function(x) {
	1 / (1 + exp(-x))
}

####

a <- d$n_bm_pos;
b <- d$n_bm_neg;
n <- a + b;
p <- a / n;
pvar <- p * (1 - p) / n;

data <- list(
	I = nrow(d),
	y = p,
	w = sqrt(pvar)
);

fit1 <- stan("random-effects_normal.stan", data = data);

####

data <- list(
	I = nrow(d),
	m = a,
	n = n
);

fit2 <- stan("random-effects_binomial.stan", data = data);
mu <- extract(fit2, "mu")[[1]];
logistic(mean(mu))
lmu <- logistic(mu);
quantile(lmu, c(0.025, 0.25, 0.5, 0.75, 0.975))

####

d <- mutate(d,
	p_luad = logit( (n_histology_luad + 1) / (n_histology_luad + n_histology_sclc + n_histology_other + 3) ),
	p_sclc = logit( (n_histology_sclc + 1) / (n_histology_luad + n_histology_sclc + n_histology_other + 3) ),
);

d2 <- filter(d, !is.na(p_luad_sclc));
X <- cbind(d2$p_luad, d2$p_sclc);

data <- list(
	I = nrow(X),
	J = ncol(X),
	m = d2$n_bm_pos,
	n = d2$n_bm_pos + d2$n_bm_neg,
	x = X
);

fit3 <- stan("random-effects_binomial_regression.stan", data = data);

mu <- extract(fit3, "mu")[[1]];
logistic(mean(mu))
lmu <- logistic(mu);
quantile(lmu, c(0.025, 0.25, 0.5, 0.75, 0.975))

plot(fit3, pars=c("mu", "tau", "beta"))

stan_scat(fit3, pars=c("mu", "beta[1]"))
stan_scat(fit3, pars=c("mu", "beta[2]"))

thetas <- extract(fit3, "theta")[[1]];
lthetas <- logistic(thetas);
summary(lthetas)
boxplot(lthetas)
d2$theta <- colMeans(lthetas);

transmute(d2,
	study, cancer_type,
	p_luad = format(logistic(p_luad), digits=2),
	p_sclc = format(logistic(p_sclc), digits=2),
	theta
)

####

X_type <- as.matrix(select(d,
	n_histology_luad, n_histology_sclc, n_histology_other
));
X_type[is.na(X_type)] <- 0;

X_stage <- select(d,
	n_stage_i_ii, n_stage_iii, n_stage_i_ii_iii, n_stage_iv
);

idx <- with(X_stage, is.na(n_stage_i_ii) & !is.na(n_stage_i_ii_iii));
X_stage[idx, ]

# distribute counts for stage I-III evenly among stage I-II and III
X_stage$n_stage_i_ii[idx] <- round(X_stage$n_stage_i_ii_iii[idx] / 2);
X_stage$n_stage_iii[idx] <- round(X_stage$n_stage_i_ii_iii[idx] / 2);
X_stage[idx, ]
X_stage <- as.matrix(select(X_stage, -n_stage_i_ii_iii));
X_stage[is.na(X_stage)] <- 0;

X_age <- as.matrix(select(d,
	n_age_lt_50, n_age_ge_50_lt_70, n_age_ge_70
));
X_age[is.na(X_age)] <- 0;

X_egfr <- as.matrix(select(d,
	n_egfr_wildtype, n_egfr_mutant
));
X_egfr[is.na(X_egfr)] <- 0;

X_race <- as.matrix(select(d,
	n_race_nfe_fin, n_race_afr, n_race_eas, n_race_sas, n_race_oth_amr
));
X_race[is.na(X_race)] <- 0;

X_smoker <- as.matrix(select(d,
	n_smoker_never, n_smoker_ever
));
X_smoker[is.na(X_smoker)] <- 0;

X_sex <- as.matrix(select(d,
	n_female, n_male
));
X_sex[is.na(X_sex)] <- 0;


data <- list(
	I = nrow(X_type),
	m = d$n_bm_pos,
	n = d$n_bm_pos + d$n_bm_neg,
	K_type = ncol(X_type),
	K_stage = ncol(X_stage),
	K_age = ncol(X_age),
	K_egfr = ncol(X_egfr),
	K_race = ncol(X_race),
	K_smoker = ncol(X_smoker),
	K_sex = ncol(X_sex),
	alpha_type = rep(1, ncol(X_type)),
	alpha_stage = rep(1, ncol(X_stage)),
	alpha_age = rep(1, ncol(X_age)),
	alpha_egfr = rep(1, ncol(X_egfr)),
	alpha_race = rep(1, ncol(X_race)),
	alpha_smoker = rep(1, ncol(X_smoker)),
	alpha_sex = rep(1, ncol(X_sex)),
	X_type = X_type,
	X_stage = X_stage,
	X_age = X_stage,
	X_egfr = X_egfr,
	X_race = X_race,
	X_smoker = X_smoker,
	X_sex = X_sex
);

str(data)

fit4 <- stan("random-effects_binomial_regression_evm.stan", data = data);

mu <- extract(fit4, "mu")[[1]];
logistic(mean(mu))
lmu <- logistic(mu);
quantile(lmu, c(0.025, 0.25, 0.5, 0.75, 0.975))

plot(fit4, pars=c("mu", "tau"));

covariates <- c("type", "stage", "age", "egfr", "race", "smoker", "sex");
plot(fit4, pars=c(paste0("beta_", covariates)));

stan_scat(fit4, pars=c("mu", "beta_type[1]"))
stan_scat(fit4, pars=c("mu", "beta_type[2]"))
stan_scat(fit4, pars=c("mu", "beta_stage[1]"))
stan_scat(fit4, pars=c("mu", "beta_stage[2]"))

stan_scat(fit4, pars=c("tau", "beta_type[1]"))
stan_scat(fit4, pars=c("tau", "beta_type[2]"))
stan_scat(fit4, pars=c("tau", "beta_stage[1]"))
stan_scat(fit4, pars=c("tau", "beta_sage[2]"))

thetas <- extract(fit4, "theta")[[1]];
lthetas <- logistic(thetas);
summary(lthetas)
boxplot(lthetas)
d$theta <- colMeans(lthetas);

transmute(d,
	study, cancer_type,
	n_luad = n_histology_luad,
	n_sclc = n_histology_sclc,
	#n_other = n_histology_other,
	theta
)

pis <- extract(fit4, "pi1")[[1]];
cbind(
	select(d, n_histology_luad, n_histology_sclc, n_histology_other),
	format(apply(pis, c(2, 3), mean), digits=2)
)

