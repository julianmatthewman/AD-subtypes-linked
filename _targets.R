library(targets)
library(tarchetypes)
library(future)
library(future.callr)
plan(callr)

#Paths to the raw data need to be supplied (path_patients & path_gp)
source("paths/paths.R")

# Source all functions from the "R" folder
sapply(list.files("R", full.names = TRUE), source, .GlobalEnv)

tar_option_set(
	packages = c(
		"haven",
		"labelled",
		"skimr",
		"readr",
		"readxl",
		"tidyverse",
		"tidymodels",
		"vip",
		"glmnet"
	)
)

list(
	# Codelists ----------
	tar_target(
		# Codelists used to extract eventdata
		codelists,
		tribble(~ name, ~ codevar, ~ extract_from,
			# readcode
			"eczema","readcode","readcode_records",
			"atopic_only_eczema","readcode","readcode_records",
			"allergic_rhinitis","readcode","readcode_records",
			"asthma","readcode","readcode_records",
			"asthma_diagnosis","readcode","readcode_records",
			"eczema_infections","readcode","readcode_records",
			"eosinophilic_eosophagitis","readcode","readcode_records",
			"folliculitis","readcode","readcode_records",
			"foodallergy","readcode","readcode_records",
			"insomnia","readcode","readcode_records",
			"urticaria","readcode","readcode_records",
			"phototherapy","readcode","readcode_records",
			# SNOMED
			"adrenaline_pens","SNOMEDCode","SNOMEDCode_records",
			"antibiotics","SNOMEDCode","SNOMEDCode_records",
			"antihistamines","SNOMEDCode","SNOMEDCode_records",
			"asthma_inhalers","SNOMEDCode","SNOMEDCode_records",
			"insomnia_drugs","SNOMEDCode","SNOMEDCode_records",
			"emollients","SNOMEDCode","SNOMEDCode_records",
			"mild_topical_corticosteroids","SNOMEDCode","SNOMEDCode_records",
			"moderate_topical_corticosteroids","SNOMEDCode","SNOMEDCode_records",
			"potent_topical_corticosteroids","SNOMEDCode","SNOMEDCode_records",
			"very_potent_topical_corticosteroids","SNOMEDCode","SNOMEDCode_records",
			"oral_corticosteroids","SNOMEDCode","SNOMEDCode_records",
			"systemic_immunosupressants","SNOMEDCode","SNOMEDCode_records",
			"topical_calcineurin_inhibitors","SNOMEDCode","SNOMEDCode_records",
			"topical_antibiotics","SNOMEDCode","SNOMEDCode_records",
		) |>
			mutate(
				path = paste0("codelists/", codevar, "/", name, "_codelist.csv"),
				full = map(path, ~ read_csv(.x)),
				codes = map2(path, codevar, ~ read_csv(.x)[[.y]]),
				codes = ifelse(
					codevar == "SNOMEDCode",
					map(codes, ~ str_sub(.x, start = 2)),
					#Take ticks of SNOMEDCodes
					map(codes, ~ gsub("[00]*$", "", .x))
				)
			) #Take 00s of ends of readcodes
	),
	
	tar_target(
		flexderm_vars,
		tibble(
			timepoint = c(6, 18, 30, 42, 57, 69, 81, 103, 128, 140, 166),
			flexderm = c(
				"kb086",
				"kd085",
				"kf110",
				"kj100",
				"kl100",
				"kn1120",
				"kq090",
				"ks1280",
				"kv1111",
				"kw1280",
				"tb1111"
			),
			severity_flexderm = c(
				"kb087",
				"kd086",
				"kf111",
				"kj101",
				"kl101",
				"kn1121",
				"kq091",
				"ks1281",
				"kv1112",
				"kw1281",
				"tb1112"
			)
		)
	),
	tar_target(timepoints, flexderm_vars$timepoint),
	
	# Get raw data ----------
	tar_target(patients_raw, read_dta(path_patients)),
	tar_target(
		patients_dict,
		patients_raw |> labelled::generate_dictionary()
	),
	tar_target(gp_raw, read_dta(path_gp)),
	
	# Extract Pateint data ------------------
	tar_target(
		patients_all,
		patients_raw |>
			rename(id = sid2606) |>
			mutate(
				sex=factor(kz021),
				social_class_mat=as.numeric(fct_recode(factor(c755), NULL="65")),
				social_class_pat=as.numeric(fct_recode(factor(c765), NULL="65")),
				social_class=pmin(social_class_mat, social_class_pat, na.rm = TRUE) |> factor(),
				ethnicity=c800,
				asthma_mat=d153 |> factor() |> fct_recode("1"="2", "2"="3") |> as.numeric(),
				asthma_pat=pa173 |> factor() |> fct_recode("1"="2", "2"="3") |> as.numeric(),
				asthma_parental=pmin(asthma_mat, asthma_pat, na.rm = TRUE) |> factor(),
				eczema_mat=d154 |> factor() |> fct_recode("1"="2", "2"="3") |> as.numeric(),
				eczema_pat=pa174 |> factor() |> fct_recode("1"="2", "2"="3") |> as.numeric(),
				eczema_parental=pmin(eczema_mat, eczema_pat, na.rm = TRUE) |> factor(),
				eczema_doctors_diag_128=if_else(kv1070 %in% c(1,2,3,4), kv1070, NA), # 1=Yes, asthma; 2=Yes, eczema; 3=Yes, asthma & eczema; 4=No
				eczema_doctors_diag_128=if_else(kv1070 %in% c(2,3), TRUE, FALSE),
				eczema_doctors_diag_166=if_else(tb1070 %in% c(1,2,3,4), tb1070, NA),
				eczema_doctors_diag_166=if_else(tb1070 %in% c(2,3), TRUE, FALSE),
				ad_severity_subtype = factor(ad_severity_subtype),
				in_gp_data = factor(in_gp_data),
				ad_severity_subtype_binary = fct_collapse(
					ad_severity_subtype,
					"1_2_3_4" = c("1", "2", "3", "4"),
					"5" = c("5")
				),
				ad_severity_subtype_tert = fct_collapse(
					ad_severity_subtype,
					"1" = c("1"),
					"2_3_4" = c("2", "3", "4"),
					"5" = "5"
				),
				ad_severity_subtype_quart = fct_collapse(
					ad_severity_subtype,
					"1" = "1",
					"2" = "2",
					"3_4" = c("3", "4"),
					"5" = "5"
				),
				ad_severity_subtype_early_gp_start= if_else(gp_start_age < 2 & gp_end_age>13, ad_severity_subtype, NA),
				ad_severity_subtype_complete_alspac= if_else(if_any(flexderm_vars$flexderm, ~ is.na(.)), NA, ad_severity_subtype),
				ad_severity_subtype_excl_rare= fct_drop(if_else(ad_severity_subtype=="5", NA, ad_severity_subtype)),
				
			)
	),
	tar_target(
		patients_non_missing, 
		patients_all |> #Keep only people that ...
			filter(!is.na(ad_severity_subtype) & ad_severity_subtype!=-9999) |>  #... have their AD phenotype recorded ...
			# filter(!if_all(flexderm_vars$flexderm, ~ is.na(.))) |> #Remove participants that don't have any questionnaire responses recorded
			filter(in_gp_data==1) |>   #... have GP data recorded ...
			droplevels()
	),
	tar_target(
		patients,
		patients_non_missing #|> filter(gp_start_age<1 & gp_end_age>13)
	),
	
	# Extract GP data --------------------
	tar_target(
		gp,
		gp_raw |>
			rename(
				id = sid2606,
				readcode = READCODE,
				SNOMEDCode = SNOMEDCT_CODE
			) |>
			mutate(eventdate = as.Date(
				paste(event_yr, event_mnth, "15", sep = "-")
			), #Assume all events happen on 15th of the month since no information on day available
			dob = as.Date(
				paste(birth_yr, birth_mnth, "15", sep = "-")
			))
	), 
	# Dummy data
	#   tar_target(
	#     patients,
	# 		make_dummy_patients(flexderm_vars) |> rename(id=sid2606)
	#   ),
	#   tar_target(
	#   	gp,
	#   	make_dummy_data(codelists) |> rename(id=sid2606)
	#   ),
	
	# Labels:
	# value: label
	# -9999: Consent withdrawn by mother
	# -11: Triplet / quadruplet
	# -10: Not completed
	# -1: No response
	# 0: Other
	# 1: Yes
	# 2: No
	# 9: Don't know
	
	# Extract GP Eventdata ----------
	
	tar_target(
		# Get eventdata for every codelist
		eventdata,
		make_eventdata(gp, timepoints, codelists),
		pattern = map(codelists),
		iteration = "list"
	),
	
	tar_target(
		eventdata_timepoints,
		make_eventdata_timepoints(eventdata, timepoints, codelists),
		pattern = map(eventdata, codelists),
		iteration = "list"
	),
	
	# Transform ---------------------------------------------------------------
	
	tar_target(
		patients_long,
		patients |>
			select(
				id,
				in_gp_data,
				#sid2606preg, gp_start_age, gp_end_age,
				#generation, alsp_cnst, hlth_cnst, in_core, in_phase2, in_phase3, in_phase4,
				ad_severity_subtype,
				ad_severity_subtype_binary,
				ad_severity_subtype_tert,
				ad_severity_subtype_quart,
				ad_severity_subtype_complete_alspac,
				eczema_doctors_diag_128,
				eczema_doctors_diag_166,
				`6-ad_presence` = kb086,
				`6-ad_symptom_severity` = kb087,
				`18-ad_presence` = kd085,
				`18-ad_symptom_severity` = kd086,
				`30-ad_presence` = kf110,
				`30-ad_symptom_severity` = kf111,
				`42-ad_presence` = kj100,
				`42-ad_symptom_severity` = kj101,
				`57-ad_presence` = kl100,
				`57-ad_symptom_severity` = kl101,
				`69-ad_presence` = kn1120,
				`69-ad_symptom_severity` = kn1121,
				`81-ad_presence` = kq090,
				`81-ad_symptom_severity` = kq091,
				`103-ad_presence` = ks1280,
				`103-ad_symptom_severity` = ks1281,
				`128-ad_presence` = kv1111,
				`128-ad_symptom_severity` = kv1112,
				`140-ad_presence` = kw1280,
				`140-ad_symptom_severity` = kw1281,
				`166-ad_presence` = tb1111,
				`166-ad_symptom_severity` = tb1112
			) |>
			pivot_longer(
				cols = ends_with(c("ad_presence", "ad_symptom_severity")),
				names_to = c("timepoint", ".value"),
				names_sep = "-"
			) |>
			mutate(timepoint = as.numeric(timepoint),
						 ad_presence = case_when(ad_presence == 1 ~ TRUE,
						 												ad_presence == 2 ~ FALSE))
	),
	tar_target(
		ecz_treatments,
		c("emollients",
			"mild_topical_corticosteroids",
			"moderate_topical_corticosteroids",
			"potent_topical_corticosteroids",
			"very_potent_topical_corticosteroids")
	),
	tar_target(
		patients_long_gp,
		join_patients_timepoints(patients_long, eventdata_timepoints[which(codelists$name %in% c("eczema", ecz_treatments))]) |>
			mutate(across(all_of(
				c("eczema", ecz_treatments)
			), \(x) replace_na(x, FALSE))) |> # Replace NAs with FALSE
			mutate(
				dx=ifelse(
						emollients == TRUE |
						mild_topical_corticosteroids == TRUE |
						moderate_topical_corticosteroids == TRUE |
						potent_topical_corticosteroids == TRUE |
						very_potent_topical_corticosteroids == TRUE,
					TRUE,
					FALSE
				),
				dx_or_rx = ifelse(eczema == TRUE | dx == TRUE, TRUE, FALSE),
				dx_and_rx = ifelse(eczema == TRUE & dx == TRUE, TRUE, FALSE)
			) # Add different definitions
	),
	
	
	# All-follow-up Variables --------------
	
	tar_target(
		ever_never_vars,
		make_ever_never_vars(patients, codelists, eventdata)
	),
	tar_target(count_vars,
						 make_count_vars(patients, codelists, eventdata)),
	tar_target(
		age_of_onset_var,
		make_age_of_onset(patients, eventdata, codelists)
	),
	tar_target(
		code_in_window_var,
		make_code_in_window(patients, eventdata, codelists)
	),
	tar_target(
		dx_and_rx_in_window_var,
		make_dx_and_rx_in_window(patients, eventdata, codelists)
	),
	tar_target(#Only variables for which codelists were provided were extracted; but also fuzzy matched for readcode, e.g. with code "J1017" for Eosinophicic Eosophagitis, also codes "J10" for Eosophagitis were returned
		readchapter_vars,
		make_readchapter_vars(gp, patients)),
	# tar_target(
	# 	BNFchapter_vars,
	# 	make_BNFchapter_vars(..., ...) #Make when structure is clearer
	# ),
	
	
	# Sensitivity & Specificity ---------------------------------------------
	
	tar_target(
		sens_spec_by_timepoint,
		patients_long_gp |>
			filter(!is.na(ad_presence)) |>
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(timepoint, definition) |>
			summarise(
				tp = sum(ad_presence == TRUE & value == TRUE),
				tn = sum(ad_presence == FALSE & value == FALSE),
				fp = sum(ad_presence == FALSE & value == TRUE),
				fn = sum(ad_presence == TRUE & value == FALSE),
				sens = tp / (tp + fn),
				spec = tn / (tn + fp)
			) |>
			mutate(timepoint = factor(timepoint))
	),
	tar_target(
		prop_timepoint_very_bad,
		patients_long_gp |>
			filter(ad_symptom_severity==1) |>
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(timepoint, definition) |>
			summarise(
				yes = sum(value == TRUE),
				no = sum(value == FALSE),
				prop = yes / (yes + no),
				pct = round(prop * 100)
			) |>
			mutate(timepoint = factor(timepoint))
	),
	tar_target(
		sens_spec_overall,
		map(c(1,2), \(number_of_reports)
		patients_long_gp |>
			filter(!is.na(ad_presence)) |>
			select(id, ad_presence, eczema, dx_or_rx, dx_and_rx) |>
			group_by(id) |>
			summarise(across(everything(), \(x) sum(as.integer(x)))) |>
			mutate(ad_presence=ad_presence>=number_of_reports) |> 
			mutate(across(all_of(c("eczema", "dx_or_rx", "dx_and_rx")), \(x) x>=1)) |> 
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(definition) |>
			summarise(
				tp = sum(ad_presence == TRUE & value == TRUE),
				tn = sum(ad_presence == FALSE & value == FALSE),
				fp = sum(ad_presence == FALSE & value == TRUE),
				fn = sum(ad_presence == TRUE & value == FALSE),
				sens = tp / (tp + fn),
				spec = tn / (tn + fp)
			))
	),
	tar_target(
		sens_spec_overall_gp_as_gold_standard,
		map(c(1,2), \(number_of_reports)
				patients_long_gp |>
					filter(!is.na(ad_presence)) |>
					select(id, ad_presence, eczema, dx_or_rx, dx_and_rx) |>
					group_by(id) |>
					summarise(across(everything(), \(x) sum(as.integer(x)))) |>
					mutate(ad_presence=ad_presence>=number_of_reports) |> 
					mutate(across(all_of(c("eczema", "dx_or_rx", "dx_and_rx")), \(x) x>=1)) |> 
					pivot_longer(
						cols = c(eczema, dx_or_rx, dx_and_rx),
						names_to = "definition",
						values_to = "value"
					) |>
					group_by(definition) |>
					summarise(
						tp = sum(value == TRUE & ad_presence == TRUE),
						tn = sum(value == FALSE & ad_presence == FALSE),
						fp = sum(value == FALSE & ad_presence == TRUE),
						fn = sum(value == TRUE & ad_presence == FALSE),
						sens = tp / (tp + fn),
						spec = tn / (tn + fp)
					))
	),
	
	tar_target(
		prop_by_subtype,
		patients_long_gp |>
			select(id, ad_severity_subtype, eczema, dx_or_rx, dx_and_rx) |>
			group_by(ad_severity_subtype, id) |>
			summarise(across(everything(), any)) |>
			ungroup() |>
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(ad_severity_subtype, definition) |>
			summarise(
				yes = sum(value == TRUE),
				no = sum(value == FALSE),
				prop = yes / (yes + no),
				pct = round(prop * 100)
			)
	),
	
	tar_target(
		prop_by_doc_diag,
		patients_long_gp |>
			select(id, eczema_doctors_diag_166, eczema, dx_or_rx, dx_and_rx) |>
			group_by(eczema_doctors_diag_166, id) |>
			summarise(across(everything(), any)) |>
			ungroup() |>
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(eczema_doctors_diag_166, definition) |>
			summarise(
				yes = sum(value == TRUE),
				no = sum(value == FALSE),
				prop = yes / (yes + no),
				pct = round(prop * 100)
			)
	),
	
	tar_target(
		doc_diag_by_subtype,
		table(patients$eczema_doctors_diag_166, patients$ad_severity_subtype)
	),
	
	tar_target(
		prop_by_doc_diag_only_severe_with_self_report,
		patients_long_gp |>
			filter(ad_severity_subtype==1 & eczema_doctors_diag_166==TRUE) |> 
			select(id, eczema, dx_or_rx, dx_and_rx) |>
			group_by(id) |>
			summarise(across(everything(), any)) |>
			ungroup() |>
			pivot_longer(
				cols = c(eczema, dx_or_rx, dx_and_rx),
				names_to = "definition",
				values_to = "value"
			) |>
			group_by(definition) |>
			summarise(
				yes = sum(value == TRUE),
				no = sum(value == FALSE),
				prop = yes / (yes + no),
				pct = round(prop * 100)
			)
	),
	
	
	# Checks ------------------------------------------------------------------
	
	tar_target(
		codelist_counts,
		codelists$full[[1]] |>
			left_join(eventdata) |>
			count(eval(sym(codelists$codevar))) |>
			set_names(c(codelists$codevar, "n")) |>
			left_join(codelists$full[[1]]),
		pattern = map(codelists, eventdata),
		iteration = "list"
	),
	
	tar_target(
		ever_never_vars_colsums,
		ever_never_vars |> select(-id) |> colSums()
	),
	tar_target(
		ever_never_vars_colsums_ecz_only,
		ever_never_vars |> select(-id) |> filter(ever_eczema == 1) |>  colSums()
	),
	tar_target(
		patients_long_summarised,
		patients_long_gp |> 
			select(id, ad_severity_subtype, timepoint, ad_presence, ad_symptom_severity, dx, dx_or_rx, dx_and_rx) |> 
			mutate(ad_symptom_severity_recoded=case_when(ad_symptom_severity==4 ~ 1,
																									 ad_symptom_severity==3 ~ 2,
																									 ad_symptom_severity==2 ~ 3,
																									 ad_symptom_severity==1 ~ 4)) |> 
			group_by(id, ad_severity_subtype) |> 
			summarise(n_ad_present=sum(as.integer(ad_presence), na.rm = TRUE),
								n_dx=sum(as.integer(dx), na.rm = TRUE),
								n_dx_or_rx=sum(as.integer(dx_or_rx), na.rm = TRUE),
								n_dx_and_rx=sum(as.integer(dx_and_rx), na.rm = TRUE),
								sum_symptom_severity=sum(as.integer(ad_symptom_severity_recoded), na.rm = TRUE),
								mean_symptom_severity=mean(ad_symptom_severity, na.rm=TRUE),
								median_timepoint_ad_present=median(timepoint[ad_presence==TRUE], na.rm = TRUE),
								mean_timepoint_ad_present=mean(timepoint[ad_presence==TRUE], na.rm = TRUE),
								median_timepoint_dx=median(timepoint[dx==TRUE], na.rm = TRUE),
								mean_timepoint_dx=mean(timepoint[dx==TRUE], na.rm = TRUE))
	),
	tar_target(
		sum_of_codes_by_phenotype,
		map2(eventdata, codelists$name, 
				 \(x,y) mutate(x, name=y )) |> 
			bind_rows() |> 
			select(id, name) |> 
			group_by(id) |> 
			count(name) |>
			ungroup() |> 
			right_join(patients[c("id", "ad_severity_subtype")]) |> 
			group_by(ad_severity_subtype, name) |> 
			summarise(sum_of_codes=sum(n)) |> 
			ungroup()
	),
	
	tar_target(	
		mean_ad_presence,
		patients_long_summarised |> 
			ungroup() |> 
			group_by(ad_severity_subtype) |> 
			summarise(mean_n_ad_present=mean(n_ad_present), 
								min_n_ad_present=min(n_ad_present), 
								max_n_ad_present=max(n_ad_present), 
								mean_mean_symptom_severity=mean(mean_symptom_severity),
								min_mean_symptom_severity=min(mean_symptom_severity),
								max_mean_symptom_severity=max(mean_symptom_severity))
		
	),
	tar_target(
		excluded,
		nrow(patients_non_missing)-nrow(patients)
	),
	tar_target(
		n_lost_no_ad_presnece, #The number of people excluded for sensitivity and specificity calculations due to not having AD presence recorded anywhere; should be zero
		length(unique(patients_long_gp$id))- length(unique(patients_long_gp[!is.na(patients_long_gp$ad_presence),]$id))
	),
	
	
	
	
	# Modeling ----------------------------------------------------------------
	
	# Specifications
	tar_target(
		outcome_spec,
		c(
			"ad_severity_subtype",
			"ad_severity_subtype_binary",
			"ad_severity_subtype_tert",
			"ad_severity_subtype_quart",
			"ad_severity_subtype_early_gp_start",
			"ad_severity_subtype_complete_alspac",
			"ad_severity_subtype_excl_rare"
		),
	),
	tar_target(
		predictor_spec_label,
		c(
		"dx_and_rx_in_window_var",
		"ever_never_vars",
		"count_vars",
		"age_of_onset_var",
		"code_in_window_var",
		"dx_and_rx_in_window_var & count_vars",
		"dx_and_rx_in_window_var & count_vars & code_in_window_var")
	),
	tar_target(
		predictor_spec,
		list(dx_and_rx_in_window_var,
				 ever_never_vars,
				 count_vars,
				 age_of_onset_var,
				 code_in_window_var,
				 left_join(dx_and_rx_in_window_var, count_vars, by="id"),
				 left_join(left_join(dx_and_rx_in_window_var, count_vars, by="id"), code_in_window_var, by="id")
				 ),
		iteration = "list"
	),
	tar_target(
		outcome_predictor_spec,
		paste(outcome_spec, predictor_spec_label),
		pattern = cross(outcome_spec, predictor_spec_label),
	),
	tar_target(
		patients_wide,
		patients[c("id", outcome_spec)] |>
			left_join(predictor_spec, by = "id") |>
			mutate(outcome = eval(sym(outcome_spec))) |>
			select(-all_of(outcome_spec)),
		pattern = cross(outcome_spec, predictor_spec),
		iteration = "list"
	),
	tar_target(
		split_data,
		initial_split(patients_wide),
		pattern = patients_wide,
		iteration = "list"
	),
	tar_target(
		n_outcomes,
		count(patients_wide, outcome),
		pattern = patients_wide,
		iteration = "list"
	),
	tar_target(
		train_data,
		training(split_data),
		pattern = split_data,
		iteration = "list"
	),
	tar_target(
		test_data,
		testing(split_data),
		pattern = split_data,
		iteration = "list"
	),
	tar_target(
		recipe,
		define_recipe(train_data),
		pattern = train_data,
		iteration = "list"
	),
	tar_target(
		model_spec,
		multinom_reg(penalty = tune(), mixture = tune()) |> set_engine("glmnet")
	),
	tar_target(
		workflow,
		workflow() |>
			add_recipe(recipe) |>
			add_model(model_spec),
		pattern = recipe,
		iteration = "list"
	),
	tar_target(
		lambda_grid,
		grid_regular(penalty(), mixture(), levels = 5)
	),
	tar_target(
		lasso_grid,
		parallel_tune(workflow, train_data, lambda_grid),
		pattern = map(workflow, train_data),
		iteration = "list"
	),
	tar_target(
		hyperparameters,
		lasso_grid |> select_best("roc_auc"),
		pattern = lasso_grid,
		iteration = "list"
	),
	tar_target(
		final_lasso,
		finalize_workflow(workflow, hyperparameters),
		pattern = map(workflow, hyperparameters),
		iteration = "list"
	),
	tar_target(
		metrics,
		last_fit(final_lasso, split_data) |> collect_metrics(),
		pattern = map(final_lasso, split_data),
		iteration = "list"
	),
	tar_target(
		augmented,
		last_fit(final_lasso, split_data) |> augment(),
		pattern = map(final_lasso, split_data),
		iteration = "list"
	),
	tar_target(
		trained_model,
		final_lasso |> fit(train_data) |> extract_fit_parsnip(),
		pattern = map(final_lasso, train_data),
		iteration = "list"
	),
	tar_target(
		model_glance,
		glance(trained_model),
		pattern = map(trained_model),
		iteration = "list"
	),
	tar_target(
		var_importance,
		trained_model |> vip::vi(lambda = hyperparameters$penalty),
		pattern = map(trained_model, hyperparameters),
		iteration = "list"
	),
	tar_target(
		confusion_matrix,
		conf_mat(augmented, outcome, .pred_class),
		pattern = map(augmented),
		iteration = "list"
	),
	tar_target(
		# Consider removing
		youdens_index,
		j_index(augmented, outcome, .pred_class),
		pattern = map(augmented),
		iteration = "list"
	),
	tar_target(
		# Consider removing
		ppv,
		ppv(augmented, outcome, .pred_class),
		pattern = map(augmented),
		iteration = "list"
	),
	tar_target(
		sens,
		sens(augmented, outcome, .pred_class),
		pattern = map(augmented),
		iteration = "list"
	),
	tar_target(
		spec,
		spec(augmented, outcome, .pred_class),
		pattern = map(augmented),
		iteration = "list"
	),
	tar_target(
		metrics_by_outcome_spec,
		pmap(list(metrics, sens, spec, model_glance), \(x,y,z,g) bind_rows(x,y,z) |> bind_cols(g)) |> 
			map2(outcome_predictor_spec, \(x,y) x |> mutate(outcome_spec=y)) |> bind_rows()
	),
	tar_target(
		# Consider removing
		ppvs,
		map_df(
			seq(
				from = 0.01,
				to = 0.1,
				length.out = 10
			),
			~ ppv(augmented, outcome, .pred_class, prevalence = .x) |> mutate(.prev = .x)
		),
		pattern = map(augmented),
		iteration = "list"
	)
)
