#' Make ever/never variables
#' @return A tibble with an id column and ever/never variables for each codelist
make_ever_never_vars <- function(patients, codelists, eventdata) {
	
    map(eventdata, \(x) 
    		as.integer(patients$id %in% x$id)) |> 
    	set_names(paste("ever",codelists$name, sep="_")) |> 
    	bind_rows() |> bind_cols(tibble(id=patients$id)) |> 
		mutate(across(everything(), \(x) replace_na(x, 0)))
}


#' Make count variables
#' @return A tibble with an id column and count variables for each codelist
make_count_vars <- function(patients, codelists, eventdata) {
	
	map2(eventdata, codelists$name, \(x,y)
			 count(x, id) |> 
			 	right_join(patients["id"], by="id") |> 
			 	arrange(id) |> 
			 	mutate(n=replace_na(n,0)) |> 
			 	select(n) |> 
			 	set_names(paste("count", y, sep = "_"))
	) |> bind_cols(tibble(id=patients$id)) |> 
		mutate(across(everything(), \(x) replace_na(x, 0)))
}


#' Make age of onset variables
#' @return A tibble with an id column and age_of_onset date variables (for only those ids where there is a code)
make_age_of_onset <- function(patients, eventdata, codelists) {
	
	map2(eventdata, codelists$name, \(x,y) {
			 	temp <- x |> 
			 		group_by(id) |> 
			 		summarise(date_of_onset=first(eventdate),
			 							dob=first(dob))
			 	
			 	patients |> 
			 		select(id) |> 
			 		left_join(temp, by = "id") |> 
			 		mutate(age_of_onset_cont=as.integer(date_of_onset-dob)/365.25,
			 					 age_of_onset=case_when(age_of_onset_cont<2 ~ "<2",
			 					 											 age_of_onset_cont<6 ~ "<6",
			 					 											 age_of_onset_cont<10 ~ "<10",
			 					 											 age_of_onset_cont<14 ~ "<14",
			 					 											 age_of_onset_cont>=14 | is.na(age_of_onset_cont) ~ "never"),
			 					 age_of_onset=factor(age_of_onset, levels = c("never", "<2", "<6", "<10", "<14" ))) |> #People with AD after 14 years also get classified as "never", since we don't want to use future information
			 		select(age_of_onset) |> 
			 		set_names(paste("age_first", y, sep = "_"))
			 }) |> bind_cols(tibble(id=patients$id))

		
}


#' Make code in window variables
#' @return A tibble with an id column and multiple code_in_window variables
make_code_in_window <- function(patients, eventdata, codelists) {
	
	map2(eventdata, codelists$name, \(x,y) {
		
		temp <- patients |> 
			select(id) |> 
			left_join(x, by = "id") |> 
			select(id, eventdate, dob) |> 
			mutate(code_at_age=as.integer(eventdate-dob)/365.25) |> 
			select(id, code_at_age) |> 
			group_by(id) |> 
			summarise(in_window_0_1=any(code_at_age>0 & code_at_age<=1),
								in_window_1_2=any(code_at_age>1 & code_at_age<=2),
								in_window_2_3=any(code_at_age>2 & code_at_age<=3),
								in_window_3_4=any(code_at_age>3 & code_at_age<=4),
								in_window_4_5=any(code_at_age>4 & code_at_age<=5),
								in_window_5_6=any(code_at_age>5 & code_at_age<=6),
								in_window_6_7=any(code_at_age>6 & code_at_age<=7),
								in_window_7_8=any(code_at_age>7 & code_at_age<=8),
								in_window_8_9=any(code_at_age>8 & code_at_age<=9),
								in_window_9_10=any(code_at_age>9 & code_at_age<=10),
								in_window_10_11=any(code_at_age>10 & code_at_age<=11),
								in_window_11_12=any(code_at_age>11 & code_at_age<=12),
								in_window_12_13=any(code_at_age>12 & code_at_age<=13),
								in_window_13_14=any(code_at_age>13 & code_at_age<=14)) |> 
			mutate(across(c(everything(),-id), \(x) replace_na(x, FALSE))) |> 
			mutate(across(c(everything(),-id), \(x) as.integer(x))) |> 
			select(-id)
			temp <- set_names(temp, paste(y, names(temp)[1:length(temp)], sep = "_"))
			temp
	}) |> bind_cols(tibble(id=patients$id))
	
}

#' Make a variable indicating if someone had a diagnosis and a treatment code in a given time window
#' @return A tibble with an id column and a dx_and_rx_in_window for each timepoint
make_dx_and_rx_in_window <- function(patients, eventdata, codelists) {
	
	temp <- bind_rows(
	eventdata[[which(codelists$name=="eczema")]] |> mutate(type="diagnosis"),
	bind_rows(eventdata[[which(codelists$name=="emollients")]],
						eventdata[[which(codelists$name=="mild_topical_corticosteroids")]], 
						eventdata[[which(codelists$name=="moderate_topical_corticosteroids")]], 
						eventdata[[which(codelists$name=="potent_topical_corticosteroids")]], 
						eventdata[[which(codelists$name=="very_potent_topical_corticosteroids")]]) |> 
		mutate(type="treatment")) |> 
		select(id, eventdate, type, dob) |> 
		arrange(id)
	
	patients |> 
		select(id) |> 
		left_join(temp, by="id") |> 
		mutate(code_at_age=as.integer(eventdate-dob)/365.25) |> 
		select(id, code_at_age, type) |> 
		group_by(id) |> 
		summarise(dx_and_rx_in_window_0_1=any(code_at_age>0 & code_at_age<=1 & type=="treatment") & any(code_at_age>0 & code_at_age<=1 & type=="diagnosis"),
							dx_and_rx_in_window_1_2=any(code_at_age>1 & code_at_age<=2 & type=="treatment") & any(code_at_age>1 & code_at_age<=2 & type=="diagnosis"),
							dx_and_rx_in_window_2_3=any(code_at_age>2 & code_at_age<=3 & type=="treatment") & any(code_at_age>2 & code_at_age<=3 & type=="diagnosis"),
							dx_and_rx_in_window_3_4=any(code_at_age>3 & code_at_age<=4 & type=="treatment") & any(code_at_age>3 & code_at_age<=4 & type=="diagnosis"),
							dx_and_rx_in_window_4_5=any(code_at_age>4 & code_at_age<=5 & type=="treatment") & any(code_at_age>4 & code_at_age<=5 & type=="diagnosis"),
							dx_and_rx_in_window_5_6=any(code_at_age>5 & code_at_age<=6 & type=="treatment") & any(code_at_age>5 & code_at_age<=6 & type=="diagnosis"),
							dx_and_rx_in_window_6_7=any(code_at_age>6 & code_at_age<=7 & type=="treatment") & any(code_at_age>6 & code_at_age<=7 & type=="diagnosis"),
							dx_and_rx_in_window_7_8=any(code_at_age>7 & code_at_age<=8 & type=="treatment") & any(code_at_age>7 & code_at_age<=8 & type=="diagnosis"),
							dx_and_rx_in_window_8_9=any(code_at_age>8 & code_at_age<=9 & type=="treatment") & any(code_at_age>8 & code_at_age<=9 & type=="diagnosis"),
							dx_and_rx_in_window_9_10=any(code_at_age>9 & code_at_age<=10 & type=="treatment") & any(code_at_age>9 & code_at_age<=10 & type=="diagnosis"),
							dx_and_rx_in_window_10_11=any(code_at_age>10 & code_at_age<=11 & type=="treatment") & any(code_at_age>10 & code_at_age<=11 & type=="diagnosis"),
							dx_and_rx_in_window_11_12=any(code_at_age>11 & code_at_age<=12 & type=="treatment") & any(code_at_age>11 & code_at_age<=12 & type=="diagnosis"),
							dx_and_rx_in_window_12_13=any(code_at_age>12 & code_at_age<=13 & type=="treatment") & any(code_at_age>12 & code_at_age<=13 & type=="diagnosis"),
							dx_and_rx_in_window_13_14=any(code_at_age>13 & code_at_age<=14 & type=="treatment") & any(code_at_age>13 & code_at_age<=14 & type=="diagnosis"),
							# Rx only
							rx_in_window_0_1=any(code_at_age>0 & code_at_age<=1 & type=="treatment"),
							rx_in_window_1_2=any(code_at_age>1 & code_at_age<=2 & type=="treatment"),
							rx_in_window_2_3=any(code_at_age>2 & code_at_age<=3 & type=="treatment"),
							rx_in_window_3_4=any(code_at_age>3 & code_at_age<=4 & type=="treatment"),
							rx_in_window_4_5=any(code_at_age>4 & code_at_age<=5 & type=="treatment"),
							rx_in_window_5_6=any(code_at_age>5 & code_at_age<=6 & type=="treatment"),
							rx_in_window_6_7=any(code_at_age>6 & code_at_age<=7 & type=="treatment"),
							rx_in_window_7_8=any(code_at_age>7 & code_at_age<=8 & type=="treatment"),
							rx_in_window_8_9=any(code_at_age>8 & code_at_age<=9 & type=="treatment"),
							rx_in_window_9_10=any(code_at_age>9 & code_at_age<=10 & type=="treatment"),
							rx_in_window_10_11=any(code_at_age>10 & code_at_age<=11 & type=="treatment"),
							rx_in_window_11_12=any(code_at_age>11 & code_at_age<=12 & type=="treatment"),
							rx_in_window_12_13=any(code_at_age>12 & code_at_age<=13 & type=="treatment"),
							rx_in_window_13_14=any(code_at_age>13 & code_at_age<=14 & type=="treatment"),
							# Dx only
							dx_in_window_0_1=any(code_at_age>0 & code_at_age<=1 & type=="diagnosis"),
							dx_in_window_1_2=any(code_at_age>1 & code_at_age<=2 & type=="diagnosis"),
							dx_in_window_2_3=any(code_at_age>2 & code_at_age<=3 & type=="diagnosis"),
							dx_in_window_3_4=any(code_at_age>3 & code_at_age<=4 & type=="diagnosis"),
							dx_in_window_4_5=any(code_at_age>4 & code_at_age<=5 & type=="diagnosis"),
							dx_in_window_5_6=any(code_at_age>5 & code_at_age<=6 & type=="diagnosis"),
							dx_in_window_6_7=any(code_at_age>6 & code_at_age<=7 & type=="diagnosis"),
							dx_in_window_7_8=any(code_at_age>7 & code_at_age<=8 & type=="diagnosis"),
							dx_in_window_8_9=any(code_at_age>8 & code_at_age<=9 & type=="diagnosis"),
							dx_in_window_9_10=any(code_at_age>9 & code_at_age<=10 & type=="diagnosis"),
							dx_in_window_10_11=any(code_at_age>10 & code_at_age<=11 & type=="diagnosis"),
							dx_in_window_11_12=any(code_at_age>11 & code_at_age<=12 & type=="diagnosis"),
							dx_in_window_12_13=any(code_at_age>12 & code_at_age<=13 & type=="diagnosis"),
							dx_in_window_13_14=any(code_at_age>13 & code_at_age<=14 & type=="diagnosis")) |> 
		mutate(across(c(everything(),-id), \(x) replace_na(x, FALSE))) |> 
		mutate(across(c(everything(),-id), \(x) as.integer(x)))
	
}


#' Make counts of codes in each chapter
#' @return
make_readchapter_vars <- function(gp, patients) {
	
	gp[!is.na(gp$readcode_match),] |> 
		mutate(readchapter=str_sub(readcode, start=1, end=3)) |> 
		group_by(id) |> 
		count(readchapter) |> 
		ungroup() |> 
		pivot_wider(names_from = readchapter, values_from = n) |> 
		right_join(tibble(id=patients$id), by="id") |> 
		mutate(across(everything(), \(x) replace_na(x, 0))) |> 
		arrange(desc(id))
}

#' Make eventdata
#' @return A dataframes containing id and timepoint of an event and a column with the variable name of all TRUE values
#' @details While events can occur multiple times at each timepoint, the returned dataframe only contains one row per id and timepoint
make_eventdata <- function(gp, timepoints, codelists) {
	
	gp |>
		filter(eval(sym(codelists$codevar)) %in% unlist(codelists$codes)) |> 
		mutate(age_at_event=as.integer(eventdate-dob),
					 timepoint=cut(age_at_event, c(0,timepoints*30), labels = timepoints), #The cut function creates a factor labelled with the timepoints
					 timepoint=as.numeric(as.character(timepoint)), #To get the factor label as a number, first convert to character, and then to numeric
					 variable=TRUE) |> 
		filter(!is.na(timepoint)) #Only keep eventdata up to the age of 14 (the last timepoint for ALSPAC subtypes)
}


#' Make eventdata with timepoints
#' @return A dataframes containing id and timepoint of an event and a column with the variable name of all TRUE values
#' @details While events can occur multiple times at each timepoint, the returned dataframe only contains one row per id and timepoint
make_eventdata_timepoints <- function(eventdata, timepoints, codelists) {
	
	eventdata |> 
		select(id, timepoint, variable) |> 
		slice(1, .by = c(id, timepoint)) |> #Keep only one row per patient and timepoint
		set_names(c("id", "timepoint", codelists$name))
}


#' Join all the eventdata to the patients
#' @return The long patient dataframe with joined variables
join_patients_timepoints <- function(patients_long, eventdata_timepoints) {
	
	reduce(c(list(patients_long), eventdata_timepoints),
				 left_join, by=c("id", "timepoint"))
}


