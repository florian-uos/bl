#' @importFrom rlang .data .env
#' @importFrom dplyr filter pull select %>%
#' @export
bl_community_index <- function(species_occurences,
                               trait_table,
                               trait,
                               ignore_missing_species = FALSE,
                               ignore_abundance = FALSE) {

    # Artenspalte finden
    species_col_name <- intersect(names(trait_table),
                                  c("Species", "Art", "species", "art"))[1]
    if (is.na(species_col_name)) {
        stop("'trait_table' must contain species column named 'Art' or 'Species'")
    }

    # Trait-Spalte vorhanden?
    if (!trait %in% names(trait_table)) {
        stop("'trait' is not a column name of 'trait_table'.")
    }

    species_names_trait      <- trait_table[[species_col_name]]
    species_names_occurences <- names(species_occurences)

    # Arten ohne Match in trait_table
    missing_species <- setdiff(species_names_occurences, species_names_trait)
    if (length(missing_species) > 0) {
        msg <- paste0("The following species in 'species_occurences' have no match in 'trait_table': ",
                      paste(missing_species, collapse = ", "))
        if (!ignore_missing_species) stop(msg)
        warning(msg)
    } else {
        message("All species in 'species_occurences' have a match in 'trait_table'")
    }

    # Vorkommen vorbereiten: NA -> 0, ggf. Presence/Absence
    occ <- as.matrix(species_occurences)
    occ[is.na(occ)] <- 0
    if (ignore_abundance) occ[] <- as.numeric(occ > 0)

    # Trait-Werte numerisch machen und an die Reihenfolge der Spalten in occ anpassen
    trait_values <- suppressWarnings(
        as.numeric(gsub(",", ".", as.character(trait_table[[trait]])))
    )
    trait_values <- trait_values[match(species_names_occurences, species_names_trait)]

    # Arten, die im trait_table stehen, aber keinen Wert haben
    no_trait <- setdiff(species_names_occurences[is.na(trait_values)], missing_species)
    if (length(no_trait) > 0) {
        warning("The following species have no trait value: ",
                paste(no_trait, collapse = ", "),
                ". These are omitted from the trait summary.")
    }

    # gewichteter Mittelwert pro Zeile (nur Arten mit Trait-Wert)
    has_trait <- !is.na(trait_values)
    occ_t     <- occ[, has_trait, drop = FALSE]
    total     <- rowSums(occ_t)
    index     <- as.numeric(occ_t %*% trait_values[has_trait]) / total
    index[total == 0] <- NA   # leere Standorte
    names(index) <- rownames(species_occurences)

    index
}
