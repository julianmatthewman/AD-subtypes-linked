make_dummy_data <- function(codelists) {
    # Specifications ----------------------------------------------------------

    # Set size
    n <- 3000 # Number of patients and rows per file

    # Set patids & pracids
    patids <- 1:344

    # Set daterange
    daterange <- as.Date(as.Date("1990-01-01"):as.Date("2007-01-01"), origin = "1970-01-01")

    # Process codes -----------------------------------------------------------

    # Make lists of all readcodes and SNOMEDcodes
    readcodes <- codelists$codes[which(codelists$codevar == "readcode")]
    SNOMEDcodes <- codelists$codes[which(codelists$codevar == "SNOMEDCode")]

    # Make readcodes that sample evenly (make a sample of each code that is the length of the longest codelist)
    longest_codelist_length <- max(map_int(c(readcodes, SNOMEDcodes), length))
    readcodes_even <- unlist(map(readcodes, sample, size = longest_codelist_length, replace = TRUE), use.names = FALSE)
    SNOMEDcodes_even <- unlist(map(SNOMEDcodes, sample, size = longest_codelist_length, replace = TRUE), use.names = FALSE)



    # Make dummy data ---------------------------------------------------------

    # Make sample function where replace=TRUE
    rsample <- function(x, size) {
        sample(x, size, replace = TRUE)
    }

    # Make dummy data
    bind_rows(
    	readcode_records = tibble(
    		sid2606 = rsample(patids, n),
    		eventdate = rsample(daterange, n),
    		READCODE = rsample(readcodes_even, n),
    		readcode_match = 1
    	),
    	SNOMEDCode_records = tibble(
    		sid2606 = rsample(patids, n),
    		eventdate = rsample(daterange, n),
    		SNOMEDCT_CODE = rsample(SNOMEDcodes_even, n),
    		snomed_match=1
    	)
    ) |> mutate(event_yr=lubridate::year(eventdate),
    						event_mnth=lubridate::month(eventdate),
    						birth_yr=1991,
    						birth_mnth=1) |> 
    	arrange(sid2606)
}

make_dummy_patients <- function(flexderm_vars, n=5000) {
	tib <- tibble(sid2606=1:n,
								sid2606preg=1:n,
								gp_start_age=0.01,
								gp_end_age=30,
								in_gp_data=rbinom(n,1, 0.9),
								ad_severity_subtype=factor(sample(c(1:5, -9999), n, replace = TRUE)),
								dob = as.Date("1990-01-01")) |> 
		mutate(ad_severity_subtype=factor(ad_severity_subtype),
					 ad_severity_subtype_new = fct_collapse(ad_severity_subtype,
					 													 onetwothree = c("1", "2", "3"),
					 													 fourfive = c("4", "5")))
	pres <- flexderm_vars$flexderm
	sev <- flexderm_vars$severity_flexderm
	tib[, c(pres, sev)] <- NA

	tib <- tib |> 
		mutate(across(all_of(pres), ~sample(c(1, 2), n, replace = TRUE)),
					 across(all_of(sev), ~factor(sample(1:4, n, replace = TRUE))))

	tib
}
