# PURPOSE: Functions to compute metrics on measurement or structure

# Function to apply over manifests of a construct and return VIF values
compute_vif <- function(target, predictors, model_data) {
  independents_regr <- stats::lm(paste("`",target,"` ~.",sep = ""),
                                 data = as.data.frame(model_data[,predictors]))

  r_squared <- summary(independents_regr)$r.squared
  1/(1 - r_squared)
}

# Function to calculate SRMR for an implied and observed cor matrix
compute_SRMR <- function(observed, implied) {

  # calculate residuals for only upper triangle
  res <- implied - observed

  # Calculate model SRMR
  ## BUT this is different from Adanco, which is different from SmarTPLS!!
  ## Henseler (2014) mentions ignoring within block residuals - is this the difference?
  sqrt(mean(res^2))
}

# BIC function using rsq, SST, n pk
BIC_func <- function(rsq, pk, N, construct_score){
  SSerrk <- (1-rsq)*(stats::var(construct_score)*(N-1))
  N*log(SSerrk/N) + (pk+1)*log(N)
}

# AIC function using rsq, SST, n pk
AIC_func <- function(rsq, pk, N, construct_score){
  SSerrk <- (1-rsq)*(stats::var(construct_score)*(N-1))
  2*(pk+1)+N*log(SSerrk/N)
}

# function to compute Henseler's rhoA
compute_construct_rhoA <- function(weights,mmMatrix,construct, obsData) {
  # get the weights for the construct
  w <- as.matrix(weights[mmMatrix[mmMatrix[,"construct"]==construct,"measurement"],construct])

  # Get empirical covariance matrix of lv indicators (S)
  indicators <- scale(obsData[,mmMatrix[mmMatrix[,"construct"]==construct,"measurement"]],TRUE,TRUE)
  S <- stats::cov(indicators,indicators)
  diag(S) <- 0

  # Get AA matrix without diagonal
  AAnondiag <- w %*% t(w)
  diag(AAnondiag) <- 0

  # Calculate rhoA
  return((t(w) %*% w)^2 * ((t(w) %*% (S) %*% w)/(t(w) %*% AAnondiag %*% w)))
}

compute_model_implied <- function(seminr_model) {

  ## Collect matrices
  S <- cor(seminr_model$data[,seminr_model$mmVariables])
  Lambda <- t(seminr_model$outer_loadings[seminr_model$mmVariables,])
  Theta  <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(Theta) <- dimnames(S)

  saturated <- cor(seminr_model$construct_scores)

  m         <- t(seminr_model$path_coef)
  Cons_endo <- rownames(m)[rowSums(m) != 0]
  Cons_exo  <- setdiff(colnames(m), Cons_endo)

  B      <- t(seminr_model$path_coef[Cons_endo, Cons_endo, drop = FALSE])
  Gamma  <- t(seminr_model$path_coef[Cons_exo,Cons_endo, drop = FALSE])
  Phi    <- saturated[Cons_exo, Cons_exo, drop = FALSE]
  I      <- diag(length(Cons_endo))


  ## Calculate variance of the zetas
  # Note: this is not yet fully correct, athough it does not currently affect
  # the results. This may have to be fixed in the future to avoid potential
  # problems that might arise in setups we have not considered yet.
  vec_zeta <- 1 - rowSums(m * saturated)
  names(vec_zeta) <- rownames(saturated)

  vcv_zeta <- matrix(0, nrow = nrow(I), ncol = ncol(I))
  diag(vcv_zeta) <- vec_zeta[Cons_endo]

  ## Correlations between exogenous and endogenous constructs
  Corr_exo_endo <- Phi %*% t(Gamma) %*% t(solve(I-B))
  ## Correlations between endogenous constructs
  Cor_endo <- solve(I-B) %*% (Gamma %*% Phi %*% t(Gamma) + vcv_zeta) %*% t(solve(I-B))
  diag(Cor_endo) <- 1

  vcv_construct <- rbind(cbind(Phi, Corr_exo_endo),
                         cbind(t(Corr_exo_endo), Cor_endo))
  ## Make symmetric
  vcv_construct[lower.tri(vcv_construct)] <- t(vcv_construct)[lower.tri(vcv_construct)]

  ## Calculate model-implied VCV of the indicators
  vcv_ind <- t(Lambda) %*% vcv_construct %*% Lambda

  Sigma <- vcv_ind + Theta

  ## Make symmetric
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

  ## Replace indicators connected to a composite by their correponding elements of S.

  mod <- seminr_model$mmMatrix
  composites <- unique(mod[mod[,"type"] == "A" | mod[,"type"] == "B","construct" ])
  index <- seminr_model$outer_loadings[,composites, drop = FALSE] %*% t(seminr_model$outer_loadings[,composites , drop = FALSE])

  Sigma[which(index != 0)] <- S[which(index != 0 )]

  # Replace indicators whose measurement errors are allowed to be correlated by s_ij
  compute_SRMR(S, Sigma)

  Sigma[.object$Information$Model$error_cor == 1] = S[.object$Information$Model$error_cor == 1]

  return(Sigma)
}



# cSEM
require(cSEM)
data(mobi)

## Note: The operator "<~" tells cSEM that the construct to its left is modelled
##       as a composite.
##       The operator "=~" tells cSEM that the construct to its left is modelled
##       as a common factor.
##       The operator "~" tells cSEM which are the dependent (left-hand side) and
##       independent variables (right-hand side).

model <- "
# Structural model
EXPE ~ IMAG
QUAL ~ EXPE
VAL  ~ EXPE + QUAL
SAT  ~ IMAG + EXPE + QUAL + VAL
COMP ~ SAT
LOY  ~ IMAG + SAT + COMP

# Composite model
IMAG <~ IMAG1 + IMAG2 + IMAG3 +IMAG4 +IMAG5
EXPE <~ CUEX1 + CUEX2 + CUEX3
QUAL <~ PERQ1 + PERQ2 + PERQ3 + PERQ4 + PERQ5 + PERQ6 +PERQ7
VAL  <~ PERV1  + PERV2
SAT  <~ CUSA1  + CUSA2  + CUSA3
COMP <~ CUSCO
LOY  <~ CUSL1 + CUSL2 +CUSL3
"
res <- csem(.data = mobi, .model = model)
res$Estimates$Loading_estimates
res$Information$Model

