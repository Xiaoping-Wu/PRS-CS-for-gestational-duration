# ===============================
# prepare target ID  
# ===============================
library(optparse)
library(dplyr)
library(data.table)


parsed_opts <- list(
  make_option("--phen",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--linkfile",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--genoID",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--phen_Mother",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--phen_Child",
              action = "store",
              default = NA,
              type = 'character')
)

if (interactive()) {
  opt <- list(
    phen="/mnt/archive/alspac/pheno/B4346/alspac-B4346-phenotype.tsv",
    linkfile="/mnt/archive/alspac/pheno/B4346/alspac-B4346-linkage.tsv",
    genoID="/mnt/archive/alspac/alspac/B4346/imputed/data/swapped.sample",
    phen_Mother="results/ALSPAC_PGS/Mother/phen.tsv",
    phen_Child="results/ALSPAC_PGS/Child/phen.tsv"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

infile <- fread(opt$phen)

phen <- infile %>%
  filter(stillbirths==2 & outcome_b_livebirth == 1 & live_births == 3 & alive_at_1_year==1) %>%            #remove stillbirth, keep live birth, keep children alove within 1 year
  filter(death_with_congenital_defect==-2) %>%         # -2 survivor
  filter(termination_for_fetal_problem==7) %>%    
  filter(pregnancy_size_summary == 1) %>%                      #remove twins
  filter(duplicate_pregnancy ==2) %>%                    # remove duplicate_pregnancy
  filter(b6s_in_vitro_fertilisation_conception != 1) %>%   # not ART
  filter(c4i_how_labour_started==1 |c4i_how_labour_started==-10) %>%                    # how labor start
  #  filter(c6a_method_of_delivery==0 ) %>%                    #delivery method, we care only "how the labors starts"
  filter(c6a_method_of_delivery != 3  & c6c_delivery_by_caesarean_section != 1 & c6c_delivery_by_caesarean_section !=2 ) %>%                                            #remove c-section
  filter(maternal_age_at_delivery>=15 & maternal_age_at_delivery<=44) %>%    #materal age
  #  filter(dv_gestation_in_days_based_on_final_clinical_estimate_of_edd>=154 & dv_gestation_in_days_based_on_final_clinical_estimate_of_edd<=310) %>%       # ga, this variable will lost many records
  filter(gestation_length>20 & gestation_length<44) %>% 
  mutate(PARITET=case_when(
    parity==0 ~ 0,
    parity>0  ~ 1,
    TRUE ~ NA)
  )%>%
  mutate(maternal_age2=maternal_age_at_delivery^2) %>%
  mutate(cidB4346= cidb4346) %>%
  select(cidB4346,PARITET,maternal_age_at_delivery,maternal_age2,sex_assigned_at_birth,dv_gestation_in_days_based_on_lmp,dv_gestation_in_days_based_on_final_clinical_estimate_of_edd,gestation_length) 

phen$c4i_how_labour_started
phen$c3c_membranes_ruptured_before_or_after_onset_of_regular_contractions



infile %>%
  filter(dv_gestation_in_days_based_on_lmp>140 & dv_gestation_in_days_based_on_lmp<310) %>%
  nrow()
#7761
infile %>%
  filter(dv_gestation_in_days_based_on_final_clinical_estimate_of_edd>140 & dv_gestation_in_days_based_on_final_clinical_estimate_of_edd<310) %>%
  nrow()
#8458
infile %>%
  filter(gestation_length>20 & gestation_length<44) %>%
  nrow()
#14063


phen$dv_gestation_in_days_based_on_final_clinical_estimate_of_edd_toweeks <- round(phen$dv_gestation_in_days_based_on_final_clinical_estimate_of_edd/7,2)
phen$gestation_length_todays <- phen$gestation_length*7


# > cor(phen$dv_gestation_in_days_based_on_final_clinical_estimate_of_edd_toweeks,phen$gestation_length)
# [1] 0.9354128
# > cor(phen$dv_gestation_in_days_based_on_final_clinical_estimate_of_edd,phen$gestation_length_todays)
# [1] 0.9354124

## 1. stillbirth: stillbirths==2
# b011	stillbirths	-9999	Consent withdrawn by mother
# b011	stillbirths	-7	HaB short
# b011	stillbirths	-1	Missing
# b011	stillbirths	1	Y
# b011	stillbirths	2	N
# b011	stillbirths	9	DK


## 2. congenital defect: death_with_congenital_defect==2
# kz017	death_with_congenital_defect	-9999	Consent withdrawn by mother
# kz017	death_with_congenital_defect	-11	Triplet / quadruplet
# kz017	death_with_congenital_defect	-7	Birth outcome not known
# kz017	death_with_congenital_defect	-3	Fetal death <20wk, not classified
# kz017	death_with_congenital_defect	-2	Survivor
# kz017	death_with_congenital_defect	-1	Not enrolled
# kz017	death_with_congenital_defect	1	Yes
# kz017	death_with_congenital_defect	2	No


## 3. fetal problem: termination_for_fetal_problem==7
# mz011b	termination_for_fetal_problem	-9999	Consent withdrawn by mother
# mz011b	termination_for_fetal_problem	-7	Unknown (new mum/pregnancy 2021)
# mz011b	termination_for_fetal_problem	-2	Outcome NK
# mz011b	termination_for_fetal_problem	1	yes, CM
# mz011b	termination_for_fetal_problem	3	yes, PRM
# mz011b	termination_for_fetal_problem	7	No

## 4. twins: pregnancy_size == 1
# mz010	pregnancy_size	-9999	Consent withdrawn by mother
# mz010	pregnancy_size	-7	Unknown (new mum/pregnancy 2021)
# mz010	pregnancy_size	-2	No delivery details
# mz010	pregnancy_size	1	Singleton
# mz010	pregnancy_size	2	Twin
# mz010	pregnancy_size	3	Triplet
# mz010	pregnancy_size	4	Quadruplet

## 5. maternal age: maternal_age_at_delivery>=15 & maternal_age_at_delivery<=44
# mz028b	maternal_age_at_delivery	-9999	Consent withdrawn by mother
# mz028b	maternal_age_at_delivery	-11	Triplet / quadruplet
# mz028b	maternal_age_at_delivery	-10	Not in core sample
# mz028b	maternal_age_at_delivery	-4	Outcome NK
# mz028b	maternal_age_at_delivery	-2	Miscarried
# mz028b	maternal_age_at_delivery	-1	Ma's DOB NK
# mz028b	maternal_age_at_delivery	15	< 16
# mz028b	maternal_age_at_delivery	44	>43

## 6. GA: dv_gestation_in_days_based_on_final_clinical_estimate_of_edd>=154 & dv_gestation_in_days_based_on_final_clinical_estimate_of_edd<=310
# DEL_P1008	dv_gestation_in_days_based_on_final_clinical_estimate_of_edd	-9999	Consent withdrawn by mother
# DEL_P1008	dv_gestation_in_days_based_on_final_clinical_estimate_of_edd	-10	Not completed
# DEL_P1008	dv_gestation_in_days_based_on_final_clinical_estimate_of_edd	-7	Unresolved
# DEL_P1008	dv_gestation_in_days_based_on_final_clinical_estimate_of_edd	-5	Exact EDD not known
# DEL_P1008	dv_gestation_in_days_based_on_final_clinical_estimate_of_edd	-1	Not answered
# DEL_P1006	dv_gestation_in_days_based_on_lmp	-9999	Consent withdrawn by mother
# DEL_P1006	dv_gestation_in_days_based_on_lmp	-10	Not completed
# DEL_P1006	dv_gestation_in_days_based_on_lmp	-7	Unresolved
# DEL_P1006	dv_gestation_in_days_based_on_lmp	-5	Exact LMP date not known
# DEL_P1006	dv_gestation_in_days_based_on_lmp	-1	Not answered
# bestgest	gestation_length	-9999	Consent withdrawn by mother
# bestgest	gestation_length	-11	Triplet/quadruplet
# bestgest	gestation_length	-10	Not in core sample
# bestgest	gestation_length	-3	Early miscarriage - Gestation not known
# bestgest	gestation_length	-2	No delivery details
# bestgest	gestation_length	-1	Missing



# 7. Spontaneously labor: c4i_how_labour_started==1
# DEL_P1160	c4i_how_labour_started	-9999	Consent withdrawn by mother
# DEL_P1160	c4i_how_labour_started	-10	Not completed
# DEL_P1160	c4i_how_labour_started	-7	Unresolved
# DEL_P1160	c4i_how_labour_started	-1	Not Answered
# DEL_P1160	c4i_how_labour_started	1	Spontaneously
# DEL_P1160	c4i_how_labour_started	2	After induction (incl failed induction)
# DEL_P1160	c4i_how_labour_started	3	No labour (e.g elective)
# DEL_P1160	c4i_how_labour_started	4	In other way

## 8. parity
# b032	parity	-9999	Consent withdrawn by mother
# b032	parity	-7	HaB short
# b032	parity	-2	Inconsistent data
# b032	parity	-1	Missing

## keep c4i_how_labour_started==1 |c4i_how_labour_started==-10
## caesarean_section related traits
# for (vv in c(
#   "c6c_delivery_by_caesarean_section",
#   "c6a_method_of_delivery",
#   "e041_had_caesarean_section"
# )) {
#   print(vv)
#   print(table(phen[[vv]], useNA = "always"))
#   print(table(phen$c4i_how_labour_started,phen[[vv]]))
#   
# }
# 
# ## These variables doesn't help  
# for (vv in c(
#   "c4iia_labour_induced_by_prostaglandin_gel_vaginal",
#   "c4iib_labour_induced_by_prostaglandin_pessaries",
#   "c4iic_labour_induced_by_extra_amniotic_prostaglandins",
#   "c4iid_labour_induced_by_oral_prostaglandins",
#   "c4iie_labour_induced_by_artificial_rupture_of_membranes_arm",
#   "c4iif_labour_induced_by_syntocinon_infusion",
#   "c4iig_labour_induced_by_other_means"
# )) {
#   print(vv)
#   # print(table(phen[[vv]], useNA = "always"))
#   print(table(phen$c4i_how_labour_started,phen[[vv]]))
#   
# }

genoID<- fread(opt$genoID) %>%
  select(ID_2) 
#all sex=1 in geno file

linkfile <- fread(opt$linkfile)
linkfile$geno_m <- paste0(linkfile$gi_hrc_g0m_g1,"M")
linkfile$geno_f <- paste0(linkfile$gi_hrc_g0m_g1,"A")


m <- phen %>%
  inner_join(linkfile,by="cidB4346")  %>%
  inner_join(genoID,by=c("geno_m"="ID_2")) %>%
  select(cidB4346,PARITET,maternal_age_at_delivery,maternal_age2,sex_assigned_at_birth,dv_gestation_in_days_based_on_lmp,dv_gestation_in_days_based_on_final_clinical_estimate_of_edd,gestation_length,geno_m)
nrow(m)
#[1] 3341

f <- phen %>%
  inner_join(linkfile,by="cidB4346")  %>%
  inner_join(genoID,by=c("geno_f"="ID_2")) %>%
  select(cidB4346,PARITET,maternal_age_at_delivery,maternal_age2,sex_assigned_at_birth,dv_gestation_in_days_based_on_lmp,dv_gestation_in_days_based_on_final_clinical_estimate_of_edd,gestation_length,geno_f)
nrow(f)
#[1] 3355

## 10. remove duplicated ID
## multiple??
# DEL_P50	mult_multiple_birth	-9999	Consent withdrawn by mother
# DEL_P50	mult_multiple_birth	-10	Not completed
# DEL_P50	mult_multiple_birth	-7	Unresolved
# DEL_P50	mult_multiple_birth	-1	Not Answered
# DEL_P50	mult_multiple_birth	1	Yes
# DEL_P50	mult_multiple_birth	2	No

# DEL_P1107	b8a_multiple_pregnancy	-9999	Consent withdrawn by mother
# DEL_P1107	b8a_multiple_pregnancy	-10	Not completed
# DEL_P1107	b8a_multiple_pregnancy	-7	Unresolved
# DEL_P1107	b8a_multiple_pregnancy	-1	Not Answered
# DEL_P1107	b8a_multiple_pregnancy	1	Yes
# DEL_P1107	b8a_multiple_pregnancy	2	No

write.table(m,opt$phen_Mother,col.names = T,row.names=F,quote=F,sep="\t")
write.table(f,opt$phen_Child,col.names = T,row.names=F,quote=F,sep="\t")
