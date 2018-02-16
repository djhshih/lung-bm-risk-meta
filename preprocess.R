library(io);
library(dplyr);
library(rstan);
library(MCMCpack);   # rdirichlet

out.fname <- filename("lung-bm-risk");
pdf.fname <- insert(out.fname, ext="pdf");

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
	X_age = X_age,
	X_egfr = X_egfr,
	X_race = X_race,
	X_smoker = X_smoker,
	X_sex = X_sex
);

str(data)

#fit4 <- qread("lung-bm-risk-meta_fit4.rds");
fit4 <- stan("random-effects_binomial_regression_evm.stan", data = data);

mu <- extract(fit4, "mu")[[1]];
logistic(mean(mu))
lmu <- logistic(mu);
quantile(lmu, c(0.025, 0.25, 0.5, 0.75, 0.975))

plot(fit4, pars=c("mu", "tau"));

covariates <- c("type", "stage", "age", "egfr", "race", "smoker", "sex");

# between-study effect dominates all fixed effects
qdraw(
	{
		plot(fit4, pars=c("tau", paste0("beta_", covariates)))
	},
	file = insert(pdf.fname, "covars")
);

stan_scat(fit4, pars=c("mu", "beta_type[1]"))
stan_scat(fit4, pars=c("mu", "beta_type[2]"))
stan_scat(fit4, pars=c("mu", "beta_stage[1]"))
stan_scat(fit4, pars=c("mu", "beta_stage[2]"))
stan_scat(fit4, pars=c("mu", "beta_age[1]"))
stan_scat(fit4, pars=c("mu", "beta_age[2]"))
stan_scat(fit4, pars=c("mu", "beta_egfr[1]"))
stan_scat(fit4, pars=c("mu", "beta_race[1]"))
stan_scat(fit4, pars=c("mu", "beta_race[2]"))
stan_scat(fit4, pars=c("mu", "beta_race[3]"))
stan_scat(fit4, pars=c("mu", "beta_smoker[1]"))
stan_scat(fit4, pars=c("mu", "beta_sex[1]"))

stan_scat(fit4, pars=c("tau", "beta_type[1]"))
stan_scat(fit4, pars=c("tau", "beta_type[2]"))
stan_scat(fit4, pars=c("tau", "beta_stage[1]"))
stan_scat(fit4, pars=c("tau", "beta_stage[2]"))
stan_scat(fit4, pars=c("tau", "beta_age[1]"))
stan_scat(fit4, pars=c("tau", "beta_age[2]"))
stan_scat(fit4, pars=c("tau", "beta_egfr[1]"))
stan_scat(fit4, pars=c("tau", "beta_race[1]"))
stan_scat(fit4, pars=c("tau", "beta_race[2]"))
stan_scat(fit4, pars=c("tau", "beta_race[3]"))
stan_scat(fit4, pars=c("tau", "beta_smoker[1]"))
stan_scat(fit4, pars=c("tau", "beta_sex[1]"))


stan_scat(fit4, pars=c("beta_stage[1]", "beta_type[1]"))

stan_scat(fit4, pars=c("beta_stage[1]", "beta_age[1]"))
stan_scat(fit4, pars=c("beta_stage[2]", "beta_age[2]"))
stan_scat(fit4, pars=c("beta_stage[1]", "beta_age[2]"))
stan_scat(fit4, pars=c("beta_stage[2]", "beta_age[1]"))

stan_scat(fit4, pars=c("beta_stage[1]", "beta_egfr[1]"))
stan_scat(fit4, pars=c("beta_stage[1]", "beta_race[1]"))
stan_scat(fit4, pars=c("beta_stage[1]", "beta_smoker[1]"))
stan_scat(fit4, pars=c("beta_stage[1]", "beta_sex[1]"))

stan_scat(fit4, pars=c("beta_stage[1]", "beta_stage[2]"))
stan_scat(fit4, pars=c("beta_age[1]", "beta_age[2]"))
stan_scat(fit4, pars=c("beta_race[1]", "beta_race[2]"))
stan_scat(fit4, pars=c("beta_race[1]", "beta_race[3]"))

param.names <- c("mu", "tau", paste0("beta_", covariates));

# mu +/- tau
# beta_type * logit(pi_type)

params <- extract(fit4, pars=param.names);

# only take covariate values from first study as an example
inputs <- list(
	X_type = data$X_type[1, ],
	X_stage = data$X_stage[1, ],
	X_age = data$X_age[1, ],
	X_egfr = data$X_egfr[1, ],
	X_race = data$X_race[1, ],
	X_smoker = data$X_smoker[1, ],
	X_sex = data$X_sex[1, ]
);

hparams <- list(
	alpha_type = data$alpha_type,
	alpha_stage = data$alpha_stage,
	alpha_age = data$alpha_age,
	alpha_egfr = data$alpha_egfr,
	alpha_race = data$alpha_race,
	alpha_smoker = data$alpha_smoker,
	alpha_sex = data$alpha_sex
);

# convert vector into a column vector (i.e. matrix with 1 column)
column_vector <- function(x) {
	matrix(x, nrow=length(x))
}

#' @return a p x 1 column vector where p is the dimension of x
transform_covariates <- function(x, alpha) {
	z <- x + alpha;
	# first beta is assumed to be fixed at 0
	column_vector(logit(z[-1] / sum(z)))
}

#' @return a n by p matrix where n is the number of random samples
#'         and p is dimension of x
transform_uncertain_covariates <- function(x, alpha, n) {
	alpha2 <- x + alpha;
	Z <- rdirichlet(n, alpha2);
	# first beta is assumed to be fixed at 0
	logit(Z[, -1] / rowSums(Z))
}

n.samples <- length(params$mu);

# inputs are based on first study
# so covars and ucovars are both based on the first study as well
covars <- mapply(transform_covariates, inputs, hparams);
ucovars <- mapply(transform_uncertain_covariates, inputs, hparams, n.samples);

qs <- c(0.025, 0.5, 0.975);
z <- extract(fit4, "z")[[1]];

betas <- with(params, list(
	beta_type,
	beta_stage,
	beta_age,
	beta_egfr,
	beta_race,
	beta_smoker,
	beta_sex
));

quantile(
	logistic(
		params$mu
	), qs
)

# calculate theta0 for a new data point from the first study
quantile(
	logistic(
		params$mu + z[, 1] * params$tau
	), qs
)

beta_covar_product <- function(betas, covars) {
	rowSums(mapply("%*%", betas, covars))
}

# calculate theta for new data point with some setting of covars
# (e.g. from first study) based on population mean,
# but ignoring uncertainty in covariates
quantile(
	logistic(
		params$mu + beta_covar_product(betas, covars)
	), qs
)

beta_ucovar_product <- function(betas, ucovars) {
	rowSums(mapply(
		function(beta, x) {
			rowSums(beta * x)
		}, betas, ucovars
	))
}

# calculate theta for new data point with some setting of covars
# (e.g. from first study) based on population mean
quantile(
	logistic(
		params$mu + beta_ucovar_product(betas, ucovars)
	), qs
)

# new data point for first study with some setting some of covars,
# ignoring uncertainty in covariates
quantile(
	logistic(
		params$mu + z[,1] * params$tau + beta_covar_product(betas, covars)
	), qs
)

# new data point for first study with some setting some of covars
quantile(
	logistic(
		params$mu + z[,1] * params$tau + beta_ucovar_product(betas, ucovars)
	), qs
)

# calculate theta0 indirectly, ignoring uncertainty in covariates
quantile(
	logistic(
		theta1 - beta_covar_product(betas, covars)
	), qs
)

theta_1 <- extract(fit4, "theta[1]")[[1]];
z_1 <- extract(fit4, "z[1]")[[1]];

# theta of first study at the observed covariates
quantile(
	logistic(theta_1), qs
)

# calculate theta0 indirectly for first study
quantile(
	logistic(
		theta_1 - beta_ucovar_product(betas, ucovars)
	), qs
)

# calculate theta0 more directly for first study
quantile(
	logistic(
		params$mu + z_1 * params$tau
	), qs
)

stan_scat(fit4, c("tau", "z[1]"))

####

thetas <- extract(fit4, "theta")[[1]];

lthetas <- logistic(thetas);
summary(lthetas)
boxplot(lthetas)
d$theta <- colMeans(lthetas);

plot(d$period_start, d$theta)
plot(d$period_end, d$theta)
plot((d$period_start + d$period_end) / 2, d$theta)

summary(lm(theta ~ period_start, data=d))


theta0s <- as.numeric(params$mu) + as.numeric(params$tau) * z;
ltheta0s <- logistic(theta0s);
boxplot(ltheta0s)
d$theta0 <- colMeans(ltheta0s);
d$theta0_lower = apply(ltheta0s, 2, quantile, probs=0.1);
d$theta0_upper = apply(ltheta0s, 2, quantile, probs=0.9);

d$period_mid <- 0.5 * (d$period_start + d$period_end);

d.main <- filter(d, !duplicated(pmid));

with(d.main, plot(period_start, theta0))
with(d.main, plot(period_end, theta0))
with(d.main, plot((period_start + period_end) / 2, theta0))

library(ggplot2);

qdraw(
	ggplot(d.main, aes(x = period_start, xend = period_end, y = theta0, yend = theta0)) +
		geom_segment() + theme_bw() +
		geom_pointrange(aes(x = period_mid, y = theta0, ymin = theta0_lower, ymax=theta0_upper))
	,
	width = 7,
	file = insert(pdf.fname, c("period-vs-theta0", "cross"))
);

qdraw(
	ggplot(d.main, aes(x = period_mid, y = theta0)) +
		geom_point() +
		geom_smooth(method="lm") + theme_bw()
	,
	width = 7,
	file = insert(pdf.fname, c("period-vs-theta0", "lm"))
);

filter(d.main, theta0 > 0.4)

filter(d, pmid == 22733536)

####

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

