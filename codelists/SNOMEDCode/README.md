# SNOMEDCode

Codelists in here were created searching the GP_snomed sheet in the linkage variable catalogue provided by ALSPAC, mapped to BNF codes to identify additional searchterms, including brand names of drugs.

The underlying code takes the DATA and keeps only rows where in the specified COLUMN at least one of the SEARCHTERMS is found and none of the EXCLUSIONTERMS is present. Using R and dplyr it looks something like this:

```
termsearch <- function(lookup, searchterms) {
    stringr::str_detect(lookup, stringr::regex(paste(searchterms, collapse = '|'), ignore_case = TRUE))
}

initial <- dplyr::filter(DATA, termsearch(COLUMN, SEARCHTERMS))
excluded <- dplyr::filter(initial, termsearch(COLUMN, EXCLUSIONTERMS))
final <- dplyr::setdiff(inital, excluded)
```

For an interactive interface, see this [Shiny app](https://julian-matthewman.shinyapps.io/codelisttools/).