
bttestISClass <- if (requireNamespace('jmvcore')) R6::R6Class(
  "bttestISClass",
  inherit = bttestISBase,
  private = list(
    .init = function() {

      ci <- paste0(self$options$ciWidth, '% Credible Interval')
      self$results$ttest$getColumn('cil')$setSuperTitle(ci)
      self$results$ttest$getColumn('ciu')$setSuperTitle(ci)

    },
    .run = function() {

      # `self$data` contains the data
      # `self$options` contains the options
      # `self$results` contains the results object (to populate)

    })
)
