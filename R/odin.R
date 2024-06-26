## Automatically generated by odin 1.5.10 - do not edit
odin_model_2sp_ <- R6::R6Class(
  "odin_model",
  cloneable = FALSE,

  private = list(
    ptr = NULL,
    use_dde = NULL,

    odin = NULL,
    variable_order = NULL,
    output_order = NULL,
    n_out = NULL,
    ynames = NULL,
    interpolate_t = NULL,
    cfuns = list(
      rhs_dde = "odin_model_2sp_rhs_dde",
      rhs_desolve = "odin_model_2sp_rhs_desolve",
      initmod_desolve = "odin_model_2sp_initmod_desolve",
      output_dde = "odin_model_2sp_output_dde"),
    dll = "ICDMM",
    user = c("aD", "age", "age_20_factor", "age_rate", "age05", "age20l",
             "age20u", "age59", "b0", "b1", "betaL", "bites_Bed_species1",
             "bites_Bed_species2", "bites_Indoors_species1",
             "bites_Indoors_species2", "cD", "chi_species1", "chi_species2",
             "cT", "cU", "custom_seasonality", "d_IRS0_1", "d_IRS0_2",
             "d_ITN0_1", "d_ITN0_2", "d1", "dB", "dCA", "dCM", "dE", "dEL",
             "delayGam", "delayMos", "den", "density_vec", "dID", "dLL",
             "dPL", "eta", "fD0", "foi_age", "ft", "gamma1", "gammaD",
             "gammaL", "het_wt", "IB0", "IC0", "ID0", "init_A", "init_D",
             "init_EL", "init_Ev", "init_IB", "init_ICA", "init_ICM",
             "init_ID", "init_Iv", "init_LL", "init_P", "init_PL", "init_S",
             "init_Sv", "init_T", "init_U", "irs_cov", "IRS_interval",
             "irs_loss", "itn_cov", "ITN_interval", "ITN_IRS_on", "itn_loss",
             "kB", "kC", "kD", "mu0", "muEL", "muLL", "muPL", "mv0", "na",
             "nh", "num_int", "omega", "p10", "p2", "phi0", "phi1", "pi",
             "PM", "Q0_species1", "Q0_species2", "r_IRS0_1", "r_IRS0_2",
             "r_ITN0_1", "r_ITN0_2", "r_ITN1", "r_ITN2", "rA", "rD",
             "rel_foi", "rP", "rT", "rU", "ssa0", "ssa1", "ssa2", "ssa3",
             "ssb1", "ssb2", "ssb3", "tau1", "tau2", "theta_c",
             "time_length", "uB", "uCA", "uD", "x_I"),

    ## This is never called, but is used to ensure that R finds our
    ## symbols that we will use from the package; without this they
    ## cannot be found by dynamic lookup now that we use the package
    ## FFI registration system.
    registration = function() {
      if (FALSE) {
        .C("odin_model_2sp_rhs_dde", package = "ICDMM")
        .C("odin_model_2sp_rhs_desolve", package = "ICDMM")
        .C("odin_model_2sp_initmod_desolve", package = "ICDMM")
        .C("odin_model_2sp_output_dde", package = "ICDMM")
      }
    },

    ## This only does something in delay models
    set_initial = function(t, y, use_dde) {
      .Call("odin_model_2sp_set_initial", private$ptr, t, y, use_dde,
            PACKAGE= "ICDMM")
    },

    update_metadata = function() {
      meta <- .Call("odin_model_2sp_metadata", private$ptr,
                    PACKAGE = "ICDMM")
      private$variable_order <- meta$variable_order
      private$output_order <- meta$output_order
      private$n_out <- meta$n_out
      private$ynames <- private$odin$make_names(
        private$variable_order, private$output_order, FALSE)
      private$interpolate_t <- meta$interpolate_t
    }
  ),

  public = list(
    initialize = function(..., user = list(...), use_dde = FALSE,
                          unused_user_action = NULL) {
      private$odin <- asNamespace("odin")
      private$ptr <- .Call("odin_model_2sp_create", user, PACKAGE = "ICDMM")
      self$set_user(user = user, unused_user_action = unused_user_action)
      private$use_dde <- use_dde
      private$update_metadata()
    },

    ir = function() {
      path_ir <- system.file("odin/odin_model_2sp.json", mustWork = TRUE,
                             package = "ICDMM")
      json <- readLines(path_ir)
      class(json) <- "json"
      json
    },

    ## Do we need to have the user-settable args here? It would be
    ## nice, but that's not super straightforward to do.
    set_user = function(..., user = list(...), unused_user_action = NULL) {
      private$odin$support_check_user(user, private$user, unused_user_action)
      .Call("odin_model_2sp_set_user", private$ptr, user, PACKAGE = "ICDMM")
      private$update_metadata()
    },

    ## This might be time sensitive and, so we can avoid computing
    ## it. I wonder if that's an optimisation we should drop for now
    ## as it does not seem generally useful. This would bring us
    ## closer to the js version which requires that we always pass the
    ## time in.
    initial = function(t) {
      .Call("odin_model_2sp_initial_conditions", private$ptr, t, PACKAGE = "ICDMM")
    },

    rhs = function(t, y) {
      .Call("odin_model_2sp_rhs_r", private$ptr, t, y, PACKAGE = "ICDMM")
    },

    deriv = function(t, y) {
      self$rhs(t, y)
    },

    contents = function() {
      .Call("odin_model_2sp_contents", private$ptr, PACKAGE = "ICDMM")
    },

    transform_variables = function(y) {
      private$odin$support_transform_variables(y, private)
    },

    engine = function() {
      "c"
    },

    run = function(t, y = NULL, ..., use_names = TRUE) {
      private$odin$wrapper_run_delay(
        self, private, t, y, ..., use_names = use_names)
    }
  ))


odin_model_2sp <- function(..., user = list(...), use_dde = FALSE,
                     unused_user_action = NULL) {
  asNamespace("odin")$deprecated_constructor_call("odin_model_2sp")
  odin_model_2sp_$new(user = user, use_dde = use_dde,
                unused_user_action = unused_user_action)
}
class(odin_model_2sp) <- "odin_generator"
attr(odin_model_2sp, "generator") <- odin_model_2sp_
odin_model_2sp_tv1_ <- R6::R6Class(
  "odin_model",
  cloneable = FALSE,

  private = list(
    ptr = NULL,
    use_dde = NULL,

    odin = NULL,
    variable_order = NULL,
    output_order = NULL,
    n_out = NULL,
    ynames = NULL,
    interpolate_t = NULL,
    cfuns = list(
      rhs_dde = "odin_model_2sp_tv1_rhs_dde",
      rhs_desolve = "odin_model_2sp_tv1_rhs_desolve",
      initmod_desolve = "odin_model_2sp_tv1_initmod_desolve",
      output_dde = "odin_model_2sp_tv1_output_dde"),
    dll = "ICDMM",
    user = c("aD", "age", "age_20_factor", "age_rate", "age05", "age20l",
             "age20u", "age59", "b0", "b1", "betaL", "bites_Bed_species1",
             "bites_Bed_species2", "bites_Indoors_species1",
             "bites_Indoors_species2", "cD", "chi_species1", "chi_species2",
             "cT", "cU", "custom_seasonality", "d_IRS0_1", "d_IRS0_2",
             "d_ITN0_1", "d_ITN0_2", "d1", "dB", "dCA", "dCM", "dE", "dEL",
             "delayGam", "delayMos", "den", "density_vec", "density_vec_sp1",
             "dID", "dLL", "dPL", "eta", "fD0", "foi_age", "ft", "gamma1",
             "gammaD", "gammaL", "het_wt", "IB0", "IC0", "ID0", "init_A",
             "init_D", "init_EL", "init_Ev", "init_IB", "init_ICA",
             "init_ICM", "init_ID", "init_Iv", "init_LL", "init_P",
             "init_PL", "init_S", "init_Sv", "init_T", "init_U", "irs_cov",
             "IRS_interval", "irs_loss", "itn_cov", "ITN_interval",
             "ITN_IRS_on", "itn_loss", "kB", "kC", "kD", "mu0", "muEL",
             "muLL", "muPL", "mv0", "na", "nh", "num_int", "omega", "p10",
             "p2", "phi0", "phi1", "pi", "PM", "Q0_species1", "Q0_species2",
             "r_IRS0_1", "r_IRS0_2", "r_ITN0_1", "r_ITN0_2", "r_ITN1",
             "r_ITN2", "rA", "rD", "rel_foi", "rP", "rT", "rU", "ssa0",
             "ssa1", "ssa2", "ssa3", "ssb1", "ssb2", "ssb3", "tau1", "tau2",
             "theta_c", "time_length", "uB", "uCA", "uD", "x_I"),

    ## This is never called, but is used to ensure that R finds our
    ## symbols that we will use from the package; without this they
    ## cannot be found by dynamic lookup now that we use the package
    ## FFI registration system.
    registration = function() {
      if (FALSE) {
        .C("odin_model_2sp_tv1_rhs_dde", package = "ICDMM")
        .C("odin_model_2sp_tv1_rhs_desolve", package = "ICDMM")
        .C("odin_model_2sp_tv1_initmod_desolve", package = "ICDMM")
        .C("odin_model_2sp_tv1_output_dde", package = "ICDMM")
      }
    },

    ## This only does something in delay models
    set_initial = function(t, y, use_dde) {
      .Call("odin_model_2sp_tv1_set_initial", private$ptr, t, y, use_dde,
            PACKAGE= "ICDMM")
    },

    update_metadata = function() {
      meta <- .Call("odin_model_2sp_tv1_metadata", private$ptr,
                    PACKAGE = "ICDMM")
      private$variable_order <- meta$variable_order
      private$output_order <- meta$output_order
      private$n_out <- meta$n_out
      private$ynames <- private$odin$make_names(
        private$variable_order, private$output_order, FALSE)
      private$interpolate_t <- meta$interpolate_t
    }
  ),

  public = list(
    initialize = function(..., user = list(...), use_dde = FALSE,
                          unused_user_action = NULL) {
      private$odin <- asNamespace("odin")
      private$ptr <- .Call("odin_model_2sp_tv1_create", user, PACKAGE = "ICDMM")
      self$set_user(user = user, unused_user_action = unused_user_action)
      private$use_dde <- use_dde
      private$update_metadata()
    },

    ir = function() {
      path_ir <- system.file("odin/odin_model_2sp_tv1.json", mustWork = TRUE,
                             package = "ICDMM")
      json <- readLines(path_ir)
      class(json) <- "json"
      json
    },

    ## Do we need to have the user-settable args here? It would be
    ## nice, but that's not super straightforward to do.
    set_user = function(..., user = list(...), unused_user_action = NULL) {
      private$odin$support_check_user(user, private$user, unused_user_action)
      .Call("odin_model_2sp_tv1_set_user", private$ptr, user, PACKAGE = "ICDMM")
      private$update_metadata()
    },

    ## This might be time sensitive and, so we can avoid computing
    ## it. I wonder if that's an optimisation we should drop for now
    ## as it does not seem generally useful. This would bring us
    ## closer to the js version which requires that we always pass the
    ## time in.
    initial = function(t) {
      .Call("odin_model_2sp_tv1_initial_conditions", private$ptr, t, PACKAGE = "ICDMM")
    },

    rhs = function(t, y) {
      .Call("odin_model_2sp_tv1_rhs_r", private$ptr, t, y, PACKAGE = "ICDMM")
    },

    deriv = function(t, y) {
      self$rhs(t, y)
    },

    contents = function() {
      .Call("odin_model_2sp_tv1_contents", private$ptr, PACKAGE = "ICDMM")
    },

    transform_variables = function(y) {
      private$odin$support_transform_variables(y, private)
    },

    engine = function() {
      "c"
    },

    run = function(t, y = NULL, ..., use_names = TRUE) {
      private$odin$wrapper_run_delay(
        self, private, t, y, ..., use_names = use_names)
    }
  ))


odin_model_2sp_tv1 <- function(..., user = list(...), use_dde = FALSE,
                     unused_user_action = NULL) {
  asNamespace("odin")$deprecated_constructor_call("odin_model_2sp_tv1")
  odin_model_2sp_tv1_$new(user = user, use_dde = use_dde,
                unused_user_action = unused_user_action)
}
class(odin_model_2sp_tv1) <- "odin_generator"
attr(odin_model_2sp_tv1, "generator") <- odin_model_2sp_tv1_
