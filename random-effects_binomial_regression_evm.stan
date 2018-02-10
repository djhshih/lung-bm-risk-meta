// Mixed meta-regression with errors-in-variable model extension

// NB ragged arrays (not available yet) could greatly simplify this code

data {
	int<lower=1> I;       // number of studies

	int<lower=0> m[I];    // numbers of successes
	int<lower=1> n[I];    // numbers of trials

	// number of groups for each categorical covariate
	int<lower=1> K_type;
	int<lower=1> K_stage;
	int<lower=1> K_age;
	int<lower=1> K_egfr;
	int<lower=1> K_race;
	int<lower=1> K_smoker;
	int<lower=1> K_sex;

	// prior counts for each group
	vector<lower=0>[K_type] alpha_type;
	vector<lower=0>[K_stage] alpha_stage;
	vector<lower=0>[K_age] alpha_age;
	vector<lower=0>[K_egfr] alpha_egfr;
	vector<lower=0>[K_race] alpha_race;
	vector<lower=0>[K_smoker] alpha_smoker;
	vector<lower=0>[K_sex] alpha_sex;

	// group counts
	int<lower=0> X_type[I, K_type];
	int<lower=0> X_stage[I, K_stage];
	int<lower=0> X_age[I, K_age];
	int<lower=0> X_egfr[I, K_egfr];
	int<lower=0> X_race[I, K_race];
	int<lower=0> X_smoker[I, K_smoker];
	int<lower=0> X_sex[I, K_sex];
}

parameters {
	real mu;              // overall effect
	real<lower=0> tau;    // between-study sd

	// regression coefficients (relative to the first group)
	real beta_type[K_type - 1];
	real beta_stage[K_stage - 1];
	real beta_age[K_age - 1];
	real beta_egfr[K_egfr - 1];
	real beta_race[K_race - 1];
	real beta_smoker[K_smoker - 1];
	real beta_sex[K_sex - 1];

	// group probabilities
	simplex[K_type] pi_type[I];
	simplex[K_stage] pi_stage[I];
	simplex[K_age] pi_age[I];
	simplex[K_egfr] pi_egfr[I];
	simplex[K_race] pi_race[I];
	simplex[K_smoker] pi_smoker[I];
	simplex[K_sex] pi_sex[I];

	real z[I];
}

transformed parameters {
	real theta[I];          // study mean effect

	for (i in 1:I) {
		real t;
		t = mu + tau*z[i];

		for (k in 1:(K_type-1)) {
			t += beta_type[k] * logit(pi_type[i, k+1]);
		}
		for (k in 1:(K_stage-1)) {
			t += beta_stage[k] * logit(pi_stage[i, k+1]);
		}
		for (k in 1:(K_age-1)) {
			t += beta_age[k] * logit(pi_age[i, k+1]);
		}
		for (k in 1:(K_egfr-1)) {
			t += beta_egfr[k] * logit(pi_egfr[i, k+1]);
		}
		for (k in 1:(K_race-1)) {
			t += beta_race[k] * logit(pi_race[i, k+1]);
		}
		for (k in 1:(K_smoker-1)) {
			t += beta_smoker[k] * logit(pi_smoker[i, k+1]);
		}
		for (k in 1:(K_sex-1)) {
			t += beta_sex[k] * logit(pi_sex[i, k+1]);
		}

		theta[i] = t;
	}
	
}

model {
	mu ~ normal(0, 5);
	tau ~ cauchy(0, 5);

	beta_type ~ normal(0, 5);
	beta_stage ~ normal(0, 5);
	beta_age ~ normal(0, 5);
	beta_egfr ~ normal(0, 5);
	beta_race ~ normal(0, 5);
	beta_smoker ~ normal(0, 5);
	beta_sex ~ normal(0, 5);

	for (i in 1:I) {
		pi_type[i] ~ dirichlet(alpha_type);
		pi_stage[i] ~ dirichlet(alpha_stage);
		pi_age[i] ~ dirichlet(alpha_age);
		pi_egfr[i] ~ dirichlet(alpha_egfr);
		pi_race[i] ~ dirichlet(alpha_race);
		pi_smoker[i] ~ dirichlet(alpha_smoker);
		pi_sex[i] ~ dirichlet(alpha_sex);

		X_type[i] ~ multinomial(pi_type[i]);
		X_stage[i] ~ multinomial(pi_stage[i]);
		X_age[i] ~ multinomial(pi_age[i]);
		X_egfr[i] ~ multinomial(pi_egfr[i]);
		X_race[i] ~ multinomial(pi_race[i]);
		X_smoker[i] ~ multinomial(pi_smoker[i]);
		X_sex[i] ~ multinomial(pi_sex[i]);
	}

	z ~ normal(0, 1);

	m ~ binomial_logit(n, theta);
}

