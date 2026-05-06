#' @export
  prepare_trait_data <- function(df, species_col = "Species") {

    df <- na.omit(df)

    trait_cols <- setdiff(names(df), species_col)

    data_scaled <- as.data.frame(scale(df[, trait_cols]))
    data_scaled[[species_col]] <- df[[species_col]]
    data_scaled <- data_scaled[, c(species_col, trait_cols)]

    return(data_scaled)
  }

#' @export

select_aci_variables <- function(data_scaled,
                                     species_col = "Species",
                                     cor_threshold = 0.30,
                                     loading_threshold = 0.70,
                                     var_threshold = 0.10,
                                     nfactors = NULL) {

      trait_data <- data_scaled[, setdiff(names(data_scaled), species_col)]

      cor_matrix <- cor(trait_data)
      diag(cor_matrix) <- NA

      selected_cor_vars <- na.omit(
        colnames(cor_matrix)[apply(abs(cor_matrix) >= cor_threshold, 2, any)]
      )

      if(is.null(nfactors)) {
        nfactors <- ncol(trait_data)
      }

      fit <- psych::principal(trait_data, rotate = "varimax", nfactors = nfactors)

      loadings_matrix <- as.matrix(fit$loadings)
      contribution <- abs(loadings_matrix)

      selected_variables <- apply(
        contribution,
        2,
        function(x) rownames(contribution)[x > loading_threshold]
      )

      prop_var <- fit$Vaccounted["Proportion Var", ]
      keep <- which(prop_var >= var_threshold)

      selected_variables_unique <- unique(unlist(selected_variables[keep]))

      out <- list(
        pca_model = fit,
        correlation_matrix = cor_matrix,
        selected_by_correlation = selected_cor_vars,
        selected_by_pca = selected_variables_unique,
        selected_table = data.frame(Variable = selected_variables_unique)
      )

      return(out)
    }


#' @export
calculate_species_aci <- function(data_scaled,
                                        variables,
                                        signs,
                                        species_col = "Species",
                                        rescale = TRUE) {

        if(length(variables) != length(signs)) {
          stop("The number of variables and signs must be the same.")
        }

        if(!all(variables %in% names(data_scaled))) {
          stop("Some variables are not present in data_scaled.")
        }

        var <- data_scaled[, variables]

        signs <- as.numeric(signs)

        ACI_raw <- rowSums(sweep(var, 2, signs, `*`)) / length(variables)

        Species_ACI <- data.frame(
          Species = data_scaled[[species_col]],
          var,
          ACI_raw = ACI_raw
        )

        if(rescale) {
          val <- abs(min(Species_ACI$ACI_raw, na.rm = TRUE))
          Species_ACI$ACI_rescaled <- Species_ACI$ACI_raw + val
        }

        return(Species_ACI)
      }


#' @export
calculate_community_aci <- function(abundance_df,
                                            species_aci,
                                            species_col = "Species",
                                            abundance_col = "Abundance",
                                            transect_col = "Site_sample",
                                            locality_col = "Locality",
                                            aci_col = "ACI_rescaled") {

          required_cols <- c(species_col, abundance_col, transect_col, locality_col)

          if(!all(required_cols %in% names(abundance_df))) {
            stop("Some required columns are missing from abundance_df.")
          }

          if(!all(c(species_col, aci_col) %in% names(species_aci))) {
            stop("species_aci must contain species_col and aci_col.")
          }

          df_joined <- abundance_df %>%
            dplyr::left_join(species_aci[, c(species_col, aci_col)], by = species_col) %>%
            na.omit()

          df_species_weighted <- df_joined %>%
            dplyr::group_by(.data[[transect_col]]) %>%
            dplyr::mutate(
              RA = .data[[abundance_col]] / sum(.data[[abundance_col]], na.rm = TRUE),
              ACIw = .data[[aci_col]] * RA,
              ACIc = sum(ACIw, na.rm = TRUE)
            ) %>%
            dplyr::ungroup()

          df_comm <- df_species_weighted %>%
            dplyr::group_by(.data[[locality_col]], .data[[transect_col]]) %>%
            dplyr::summarise(
              ACIc = dplyr::first(ACIc),
              total_abundance = sum(.data[[abundance_col]], na.rm = TRUE),
              richness = dplyr::n_distinct(.data[[species_col]]),
              .groups = "drop"
            )

          return(list(
            species_weighted = df_species_weighted,
            community_transect = df_comm
          ))
        }

