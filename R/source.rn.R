#----------------------------------------------------------------
#	QGglmm R package (reaction norms)
#	Functions to estimate QG parameters from non-linear reaction norms
#----------------------------------------------------------------
#	Pierre de Villemereuil (2023)
#----------------------------------------------------------------

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

##  ---------------------------- Helper functions ------------------------- ##

## Non-Bessel-corrected variance
var_nocorrect <- function(vec) {
    N <- length(vec)

    ((N - 1) / N) * var(vec)
}

## Compute the average phenotype conditional to the environment
# Args: - e: the environmental value to condition to
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
# Value: The value for E_g_e (numeric)
rn_avg_e <- function(e, shape, theta, G_theta, width = 10, fixed = NA) {
    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Average
    avg <- cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            shape(e, full_x) * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral

    return(avg)
}

## Compute the genetic variance conditionnally to the environment
# Args: - e: the environmental value to condition to
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
# Value: The value for V_g_e (numeric)
rn_vg_e <- function(e, shape, theta, G_theta, width = 10, fixed = NA) {

    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Average
    avg <- rn_avg_e(e, shape, theta, G_theta, width = 10, fixed = fixed)

    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            shape(e, full_x)^2 * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral - avg^2
}

## Compute the additive genetic variance conditionnally to the environment
# Args: - e: the environmental value to condition to
#       - d_shape: the derivative of the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
# Value: The value for V_A_e (numeric)
rn_psi <- function(e, d_shape, theta, G_theta, width = 10, fixed = NA) {
    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Computing the integral for Psi
    psi <- cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            d_shape(e, full_x) * matrix(rep(vec_mvnorm(x, var_theta, G_theta, logdet), d),
                                      nrow = d,
                                      byrow = TRUE)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = d,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral

    # Naming the vector if names are available
    names(psi) <- names(theta[-fixed])

    # Now, computing V_A_e
    return(psi)
}

##  ---------------------------- General functions ------------------------- ##

## Utilitary function (exposed to the user) to compute a derivative automatically
# Args: - expr: expression from which the derivative needs to the computed (expression)
#       - dpars: the parameters wrt which we need to compute the derivative (character)
#       - allpars: list of all the parameters (character)
# Value: A function that can be used to compute V_A using QGrn_va()
QGrn_generate_derivative <- function(expr, dpars, allpars) {
    # Generating the new parameter "names" for the function to provide QGglmm
    matpars <- str_c("pars[",1:length(allpars),", ]")
    names(matpars) <- allpars

    # Computing the derivative
    deriv <-
        lapply(dpars, \(p) D(expr, p)) |>
        as.character() |>
        str_replace_all(matpars)

    # Constructing the body of the function
    body <- parse(text = str_c(c("matrix(c(",
                                          str_c(deriv, collapse = ","),
                                          "), nrow = 2, byrow = TRUE)"),
                                        collapse = ""))

    # Generating the list of arguments without defaults
    args <- list()
    args[["x"]] <- alist(x=)$x
    args[["pars"]] <- alist(pars=)$pars

    # Return the function to provide QGglmm
    eval(call("function", as.pairlist(args), body[[1]]), parent.frame())
}


## Compute the mean phenotype by environment
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
# Value: The value for V_plas (numeric)
QGrn_mean_by_env <- function(env, shape, theta, G_theta, width = 10, fixed = NA) {
    sapply(env,
            \(e) rn_avg_e(e = e, shape = shape, theta = theta, G_theta = G_theta,
                          width = width, fixed = fixed))
}

## Compute the plastic variance (V_plas)
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - correction: should the Bessel correction be used or not?
# Value: The value for V_plas (numeric)
QGrn_vplas <- function(env, shape, theta, G_theta, width = 10, fixed = NA, correction = FALSE) {
    # Should we apply Bessel's correction (R default) or not (QGrn default)?
    if (correction) {
        var_func <- var
    } else {
        var_func <- var_nocorrect
    }

    # Compute the average for each environment, then take the variance
    sapply(env,
            \(e) rn_avg_e(e = e, shape = shape, theta = theta, G_theta = G_theta,
                          width = width, fixed = fixed)) |>
        var_func()
}

## Compute the genetic variance (V_gen)
# Args: - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - env: the environmental values over which the model has been estimated
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
# Value: The value for V_gen (numeric)
QGrn_vgen <- function(shape, theta, G_theta, env, width = 10, fixed = NA, average = TRUE) {

    # Computing V_gen for each environment
    out <-
        sapply(env,
        \(e) rn_vg_e(e = e, shape = shape, theta = theta, G_theta = G_theta,
                            width = width, fixed = fixed))

    # Averaging if requested (default)
    if (average) { out <- mean(out) }

    return(out)
}

## Compute the additive genetic variance (V_A)
# Args: - env: the environmental values over which the model has been estimated
#       - d_shape: the derivative of function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
# Value: The value for V_A and Gamma-decomposition as a data.frame for each environment or averaged
QGrn_va <- function(env, d_shape, theta, G_theta, width = 10, fixed = NA, average = TRUE) {

    # Computing psi for each environment
    psi <-
        lapply(env,
               \(e) {rn_psi(e = e, d_shape = d_shape, theta = theta, G_theta = G_theta,
                           width = width, fixed = fixed)})
    psi <- do.call("rbind", psi)

    # Computing V_A for each environment
    v_a <- apply(psi, 1, \(psi_) t(psi_) %*% G_theta %*% psi_)

    # Computing the gamma-decomposition of V_A
    gamma_i <- t(apply(psi, 1, \(psi_) psi_^2 * diag(G_theta)))
    colnames(gamma_i) <- paste0("Gamma:",colnames(gamma_i))

    gamma_ij <- apply(psi, 1, simplify = FALSE, FUN = \(psi_) {
        out <- numeric(sum(upper.tri(G_theta)))
        names_out <- character(sum(upper.tri(G_theta)))
        k   <- 1
        for (i in 1:(length(psi_) - 1)) {
            for (j in (i+1):length(psi_)) {
                out[k] <- 2 * psi_[i] * psi_[2] * G_theta[i, j]
                names_out[k] <- paste(names(psi_)[i], names(psi_)[j], sep = "-")
                k <- k + 1
            }
        }
        names(out) <- names_out
        return(out)
    })
    gamma_ij <- do.call("rbind", gamma_ij)
    colnames(gamma_ij) <- paste0("Gamma:",colnames(gamma_ij))

    out <-
        cbind(
            gamma_i,
            gamma_ij
        )

    # Averaging if requested (default)
    if (average) {
        v_a <- mean(v_a)
        out <- colMeans(out) / v_a
        out <- cbind(data.frame(V_A = v_a), t(as.data.frame(out)))
        rownames(out) <- NULL
    } else {
        out <- out / v_a
        out <- cbind(data.frame(Env = env, V_A = v_a), out)
    }

    return(out)
}
