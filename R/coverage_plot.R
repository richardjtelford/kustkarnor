#'coverage_plot
#'
#'@description  A simple diagnostic plot showing the coverage of fossil taxa in modern calibration set
#'
#'@param spp data.frame of modern species abundances
#'@param fos data.frame of fossil species abundances
#'@param n2_rare numeric value of Hill's N2 below which species are highlighted as rare
#'@param label numeric label taxa where maximum fossil abundance - maximum modern abundance > label. Defaults to NULL which does not add labels
#'
#' @details   Finds the maximum abundance of fossil taxa and plots this against the maximum abundance the taxa in the modern calibration set. Taxa with a Hill's N2 less than \code{rare} in the calibration set are highlighted in blue. Taxa absent from the calibration set are highlighed in red. If there are many taxa above the 1:1 line, or important fossil taxa have a low N2 in the calibration set, reconstructions should be interpreted with caution.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' data("SWAP", package = "rioja")
#' data("RLGH", package = "rioja")
#' coverage_plot(spp = SWAP$spec, fos=RLGH$spec, n2_rare = 5, label = 0)
#'
#' @importFrom dplyr bind_rows filter select mutate if_else group_by inner_join summarise_all
#' @importFrom tidyr spread gather replace_na
#' @importFrom tibble enframe
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs scale_colour_brewer
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_explicit_na fct_relevel
#' @importFrom utils data
#' @importFrom stats coef
#' @export

coverage_plot <- function(spp, fos, n2_rare = 5, label = NULL){
  mod_fos <- bind_rows(fos = fos, spp = spp, .id = "data")

  #calculate N2
  N2 <- mod_fos %>%
    filter(.data$data == "spp") %>%
    ungroup() %>%
    select_if(is.numeric) %>%
    Hill.N2() %>%
    enframe(name = "Taxon", value = "n2") %>%
    #missing spp get n2 = Inf. replace with NA
    mutate(n2 = if_else(is.infinite(.data$n2), NA_real_, .data$n2),
           n2_cut = cut(.data$n2,
                        breaks = c(0, n2_rare, Inf),
                        labels = paste("N2", c("<=", ">"), n2_rare)),
           n2_cut = fct_explicit_na(.data$n2_cut, na_level = "-"),
           n2_cut = fct_relevel(.data$n2_cut, "-"))

  #find max and join to N2
  max_n2 <- mod_fos %>%
    group_by(data, add = TRUE) %>%
    summarise_all(max) %>%
    group_by_if(is.character) %>%
    pivot_longer(cols = -one_of(group_vars(.)), names_to = "Taxon", values_to = "max") %>%
    ungroup()

  max_n2 <- left_join(
    max_n2 %>% filter(data == "fos") %>% select(-data, fos = max),
    max_n2 %>% filter(data == "spp") %>% select(Taxon, spp = max),
    by = "Taxon") %>%
    replace_na(replace = list(spp = 0)) %>%
    #filter out double missing taxa
    filter(!(fos == 0 & spp == 0)) %>%
    inner_join(N2, by = "Taxon")


  #plot
  g <- ggplot(max_n2, aes(x = .data$spp, y = .data$fos,
                          colour = .data$n2_cut)) +
    geom_abline() +
    geom_point() +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Maximum modern %", y = "Maximum fossil %", colour = "N2")

  if(!is.null(label)){
    to_label <- max_n2 %>%
      filter(.data$fos - .data$spp > label)
    g <- g +
      geom_text_repel(data = to_label,
                      mapping = aes(label = .data$Taxon),
                      show.legend = FALSE,
                      colour = "black", size = 3)
  }

  return(g)
}
