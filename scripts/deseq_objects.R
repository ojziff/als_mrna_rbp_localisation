# run this script in Rstudio server OnDemand Terminal Tab (not with conda env r4.0.3) with:
# cd /camp/home/ziffo/home/projects/motor-neuron-mislocalisation/scripts
# sbatch -N 1 -c 12 --mem=0 -t 3:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/scripts/ml240_sample_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 12 --mem=72G -t 8:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/scripts/ml240_sample_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 6 --mem=50G -t 3:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/scripts/ml240_sample_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects 
# sbatch --cpus-per-task=10 --mem=1500G -N 1 --partition=hmem -t 3:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/scripts/ml240_sample_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects

# .libPaths("/camp/lab/luscomben/home/users/ziffo/R/x86_64-pc-linux-gnu-library/4.0")
# .libPaths("/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(here) # set here path https://cran.r-project.org/web/packages/here/vignettes/here.html
# camp_path = here("/Volumes/lab-luscomben/home/users/ziffo")
camp_path = here("/camp/lab/luscomben/home/users/ziffo")
# shared_path = here("/Volumes/lab-luscomben/home/users/ziffo/proj-luscombn-patani/working")
shared_path = here("/camp/project/proj-luscombn-patani/working")
# collab_path = here("/Volumes/lab-luscomben/home/users/ziffo/patani-collab")
collab_path = here("/camp/lab/luscomben/home/users/ziffo/patani-collab")
# proj_path = here("/Volumes/lab-luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation")
proj_path = here("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation")
# load(here(proj_path,"scripts/ml240.RData")) # overwrites with old functions
# load(here(proj_path,"sample-details/metadata.RData"))
# load(here(proj_path,"scripts/untreated_deseq_objects.RData"))
# load(here(proj_path,"scripts/ml240_sample_deseq_objects.RData"))
# source("/Volumes/lab-luscomben/home/users/ziffo/scripts/functions/R_workspace.R") 
source(here(camp_path,"scripts/functions/R_workspace.R"))


# # Metadata ----------------------------------------------------------------
all.metadata = read_csv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  filter(sample != "ctrl3_untreated_nucleus") %>% #failed nucleus QC - clusters with whole so remove
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-inhibitor-harley-2022"), mutation = case_when(cellline %in% c("ctrl1", "ctrl3", "ctrl4", "ctrl6") ~ "ctrl", cellline %in% c("glia", "glib", "cb1d", "cb1e", "ncrme6", "ncrmc2") ~ "vcp", cellline %in% c("tdpn1", "tdpn2") ~ "tardbp"), 
         genotype = case_when(cellline %in% c("cb1d", "cb1e") ~ "r155c", cellline %in% c("glia", "glib", "ncrme6", "ncrmc2") ~ "r191q", cellline %in% c("tdpn1", "tdpn2") ~ "g298s", TRUE ~ "ctrl"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), 
         isogenic = case_when(cellline %in% c("ncrme6", "ncrmc2") ~ "isogenic", mutation == "ctrl" ~ "ctrl", TRUE ~ "no"), parent_cellline = case_when(cellline %in% c("ncrmc2", "ncrme6") ~ "ctrl6", TRUE ~ cellline),
         treatment = factor(treatment, levels = c("untreated", "ml240")), fraction = gsub("nuclear","nucleus",fraction), fraction = factor(fraction, levels = c("cytoplasm", "nucleus", "whole"))) %>%
  select(sample, fraction, treatment, condition, mutation, cellline, genotype, isogenic, parent_cellline, database_dir)

metadata = read_csv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  filter(fraction != "whole", sample != "ctrl3_untreated_nucleus") %>% #failed nucleus QC - clusters with whole so remove
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-inhibitor-harley-2022"), mutation = case_when(cellline %in% c("ctrl1", "ctrl3", "ctrl4", "ctrl6") ~ "ctrl", cellline %in% c("glia", "glib", "cb1d", "cb1e", "ncrme6", "ncrmc2") ~ "vcp", cellline %in% c("tdpn1", "tdpn2") ~ "tardbp"), 
         genotype = case_when(cellline %in% c("cb1d", "cb1e") ~ "r155c", cellline %in% c("glia", "glib", "ncrme6", "ncrmc2") ~ "r191q", cellline %in% c("tdpn1", "tdpn2") ~ "g298s", TRUE ~ "ctrl"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), 
         isogenic = case_when(cellline %in% c("ncrme6", "ncrmc2") ~ "isogenic", mutation == "ctrl" ~ "ctrl", TRUE ~ "no"), parent_cellline = case_when(cellline %in% c("ncrmc2", "ncrme6") ~ "ctrl6", TRUE ~ cellline),
         treatment = factor(treatment, levels = c("untreated", "ml240")), fraction = gsub("nuclear","nucleus",fraction), fraction = factor(fraction, levels = c("cytoplasm", "nucleus"))) %>%
  select(sample, fraction, treatment, condition, mutation, cellline, genotype, isogenic, parent_cellline, database_dir)

ctrl.metadata = metadata %>% filter(mutation == "ctrl")
# vcp.metadata = metadata %>% filter(mutation == "vcp")
# tardbp.metadata = metadata %>% filter(mutation == "tardbp")
# als.metadata = metadata %>% filter(condition == "als")

untreated.metadata <- metadata %>% filter(treatment == "untreated")
ctrl.untreated.metadata <- untreated.metadata %>% filter(mutation == "ctrl")
# tardbp.untreated.metadata = untreated.metadata %>% filter(mutation == "tardbp")
# vcp.untreated.metadata = untreated.metadata %>% filter(mutation == "vcp")

# nucleus.untreated.metadata <- untreated.metadata %>% filter(fraction == "nucleus")
# cytoplasm.untreated.metadata <- untreated.metadata %>% filter(fraction == "cytoplasm")
# untreated.vcp.iso_noiso.metadata = untreated.metadata %>% filter(mutation %in% c("ctrl","vcp")) %>% group_by(mutation) %>% mutate(denom=-1/(n()-1), isogenic_pair = ifelse(parent_cellline=="ctrl6", 1, denom)) %>% ungroup() # account for isogenic pairing

ml240.metadata <- metadata %>% filter(treatment == "ml240")
# nucleus.ml240.metadata <- ml240.metadata %>% filter(fraction == "nucleus")
# cytoplasm.ml240.metadata <- ml240.metadata %>% filter(fraction == "cytoplasm")
# fraction.ml240.vcp.iso_noiso.metadata = ml240.metadata %>% filter(mutation %in% c("ctrl","vcp")) %>% group_by(mutation) %>% mutate(denom=-1/(n()-1), isogenic_pair = ifelse(parent_cellline=="ctrl6", 1, denom)) %>% ungroup() # account for isogenic pairing

whole.metadata = read_csv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  filter(fraction == "whole") %>% 
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-inhibitor-harley-2022"), mutation = case_when(cellline %in% c("ctrl1", "ctrl3", "ctrl4", "ctrl6") ~ "ctrl", cellline %in% c("glia", "glib", "cb1d", "cb1e", "ncrme6", "ncrmc2") ~ "vcp", cellline %in% c("tdpn1", "tdpn2") ~ "tardbp"), 
         genotype = case_when(cellline %in% c("cb1d", "cb1e") ~ "r155c", cellline %in% c("glia", "glib", "ncrme6", "ncrmc2") ~ "r191q", cellline %in% c("tdpn1", "tdpn2") ~ "g298s", TRUE ~ "ctrl"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), 
         isogenic = case_when(cellline %in% c("ncrme6", "ncrmc2") ~ "isogenic", mutation == "ctrl" ~ "ctrl", TRUE ~ "no"), parent_cellline = case_when(cellline %in% c("ncrmc2", "ncrme6") ~ "ctrl6", TRUE ~ cellline), treatment = factor(treatment, levels = c("untreated", "ml240"))) %>%
  select(sample, treatment, condition, mutation, cellline, genotype, isogenic, parent_cellline, database_dir)
whole.untreated.metadata = whole.metadata %>% filter(treatment == "untreated")
whole.ml240.metadata = whole.metadata %>% filter(treatment == "ml240")

# read_tsv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D', 'XIST')) %>% select(-gene_id) %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% 
#   mutate(male_counts = KDM5D + DDX3Y + RPS4Y1, female_counts = XIST, cellline = str_split_fixed(sample, "_", 3)[,1]) %>%
#   # filter(grepl("untreated",sample), grepl("whole",sample)) %>%
#   ggplot(aes(x = male_counts, y= female_counts)) + geom_point() + geom_text_repel(aes(label = cellline), max.overlaps = 15, size = 2.3) + theme_oz() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# 
# read_tsv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D', 'XIST')) %>% select(-gene_id) %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% 
#   mutate(male_counts = KDM5D + DDX3Y + RPS4Y1, female_counts = XIST, cellline = str_split_fixed(sample, "_", 3)[,1]) %>%
#   filter(grepl("tdpn1",sample)) %>%
#   ggplot(aes(x = male_counts, y= female_counts)) + geom_point() + geom_text_repel(aes(label = sample), max.overlaps = 15, size = 2.3) + theme_oz() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

# # Markmiller HNRNP mutant
# markmiller.metadata = read_csv(here(shared_path, "public-data/motor-neuron-stress-ipsc-markmiller-2021/sample-details/samplesheet.csv")) %>% 
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/motor-neuron-stress-ipsc-markmiller-2021",
#          condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), mutation = as.factor(mutation), fraction = as.factor(fraction)) %>%
#   select(sample, fraction, treatment, condition, mutation, cellline, database_dir)
# hnrnpa2b1.metadata = markmiller.metadata %>% filter(treatment == "untreated")

print("saving metadata")
save(all.metadata, metadata, ctrl.metadata, untreated.metadata, ctrl.untreated.metadata, ml240.metadata, whole.metadata, whole.untreated.metadata, whole.ml240.metadata,
     file = here(proj_path, "scripts/metadata_deseq_objects.RData"))
print("saved metadata")


# DESeq2 ----------------------------------------------------------------
# untreated.tx = DESeq.analysis(metadata = filter(all.metadata, treatment == "untreated"), design = ~ 1, transcript.level = TRUE, gene.level = FALSE, vsd_tx = TRUE)

# ### untreated  ----------------------------
untreated.nuc_vs_cyt = DESeq.analysis(metadata = untreated.metadata, design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)
untreated.nuc_vs_cyt.tx = DESeq.analysis(metadata = untreated.metadata, design = ~ 1, transcript.level = TRUE, gene.level = FALSE, vsd_tx = TRUE)

# nuc vs cyt
ctrl.untreated.nuc_vs_cyt = DESeq.analysis(metadata = filter(untreated.metadata, mutation == "ctrl"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
tardbp.untreated.nuc_vs_cyt = DESeq.analysis(metadata = filter(untreated.metadata, mutation == "tardbp"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
vcp_r155c.untreated.nuc_vs_cyt = DESeq.analysis(metadata = filter(untreated.metadata, genotype == "r155c", isogenic == "no"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
vcp_r191q.untreated.nuc_vs_cyt = DESeq.analysis(metadata = filter(untreated.metadata, genotype == "r191q", isogenic == "no"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
vcp_iso.untreated.nuc_vs_cyt = DESeq.analysis(metadata = filter(untreated.metadata, mutation == "vcp", isogenic == "isogenic"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")

# als vs ctrl
nuc.untreated.tardbp_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, mutation %in% c("ctrl", "tardbp"), fraction == "nucleus"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
cyt.untreated.tardbp_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, mutation %in% c("ctrl", "tardbp"), fraction == "cytoplasm"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
nuc.untreated.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl"), fraction == "nucleus"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
cyt.untreated.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl"), fraction == "cytoplasm"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
nuc.untreated.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl"), fraction == "nucleus"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
cyt.untreated.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl"), fraction == "cytoplasm"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
nuc.untreated.vcp_iso_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl"), fraction == "nucleus"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)
cyt.untreated.vcp_iso_vs_ctrl = DESeq.analysis(metadata = filter(untreated.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl"), fraction == "cytoplasm"),  design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE)

# Dual factor: ALS vs CTRL & nuc vs cyt
untreated.nuc_vs_cyt.tardbp_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)))
untreated.nuc_vs_cyt.tardbp_vs_ctrl.design_matrix = untreated.nuc_vs_cyt.tardbp_vs_ctrl.mm[, !apply(untreated.nuc_vs_cyt.tardbp_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
untreated.nuc_vs_cyt.tardbp_vs_ctrl = DESeq.analysis(metadata = mutate(filter(untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)), 
                                                     design = untreated.nuc_vs_cyt.tardbp_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
                                                     contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
                                                                  "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationtardbp.fractionnucleus.IRFinderIR")
untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(untreated.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)))
untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl.design_matrix = untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm[, !apply(untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = mutate(filter(untreated.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)), 
                                                     design = untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
                                                     contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
                                                                  "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_r155c.fractionnucleus.IRFinderIR")
untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)))
untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl.design_matrix = untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm[, !apply(untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = mutate(filter(untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)), 
                                                        design = untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
                                                        contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
                                                                     "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_r191q.fractionnucleus.IRFinderIR")
untreated.nuc_vs_cyt.vcp_iso_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(untreated.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl")), line_in_condition = recode_within(cellline, condition)))
untreated.nuc_vs_cyt.vcp_iso_vs_ctrl.design_matrix = untreated.nuc_vs_cyt.vcp_iso_vs_ctrl.mm[, !apply(untreated.nuc_vs_cyt.vcp_iso_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
untreated.nuc_vs_cyt.vcp_iso_vs_ctrl = DESeq.analysis(metadata = mutate(filter(untreated.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl")), line_in_condition = recode_within(cellline, condition)), 
                                                        design = untreated.nuc_vs_cyt.vcp_iso_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
                                                        contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
                                                                     "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_iso.fractionnucleus.IRFinderIR")

print("saving untreated")
save(untreated.nuc_vs_cyt, untreated.nuc_vs_cyt.tx, ctrl.untreated.nuc_vs_cyt, tardbp.untreated.nuc_vs_cyt, vcp_r155c.untreated.nuc_vs_cyt, vcp_r191q.untreated.nuc_vs_cyt, vcp_iso.untreated.nuc_vs_cyt, #hnrnpa2b1.untreated.nuc_vs_cyt, #als.untreated.nuc_vs_cyt,
     nuc.untreated.tardbp_vs_ctrl, cyt.untreated.tardbp_vs_ctrl, nuc.untreated.vcp_r155c_vs_ctrl, cyt.untreated.vcp_r155c_vs_ctrl, nuc.untreated.vcp_r191q_vs_ctrl, cyt.untreated.vcp_r191q_vs_ctrl, nuc.untreated.vcp_iso_vs_ctrl, cyt.untreated.vcp_iso_vs_ctrl,
     untreated.nuc_vs_cyt.vcp_r191q_vs_ctrl, untreated.nuc_vs_cyt.vcp_r155c_vs_ctrl, untreated.nuc_vs_cyt.vcp_iso_vs_ctrl, untreated.nuc_vs_cyt.tardbp_vs_ctrl, #untreated.nuc_vs_cyt.als_vs_ctrl, #untreated.nuc_vs_cyt.hnrnpa2b1_vs_ctrl,
     file = here(proj_path, "scripts/untreated_sample_deseq_objects.RData"))
print("complete untreated")

# combine mutants
untreated.nuc_vs_cyt.als_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(untreated.metadata, line_in_condition = recode_within(cellline, condition)))
untreated.nuc_vs_cyt.als_vs_ctrl.design_matrix = untreated.nuc_vs_cyt.als_vs_ctrl.mm[, !apply(untreated.nuc_vs_cyt.als_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
untreated.nuc_vs_cyt.als_vs_ctrl = DESeq.analysis(metadata = mutate(untreated.metadata, line_in_condition = recode_within(cellline, condition)), design = untreated.nuc_vs_cyt.als_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
                                                  contrast = c("Intercept"=0, "fractionnucleus"=0, "conditionals"=0, "conditionctrl.line_in_condition2"=0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0, "conditionals.line_in_condition3"=0, 
                                                               "conditionctrl.line_in_condition4"=0, "conditionals.line_in_condition4"=0, "conditionals.line_in_condition5"=0, "conditionals.line_in_condition6"=0, "conditionals.line_in_condition7"=0, 
                                                               "conditionals.line_in_condition8"=0, "fractionnucleus.conditionals"=1))
saveRDS(untreated.nuc_vs_cyt.als_vs_ctrl, file = here(proj_path, "scripts/untreated.nuc_vs_cyt.als_vs_ctrl.rds"))

# postmortem nygc transcript level
nygc.gene_tpm.male_counts = read_tsv(here(shared_path,"nygc-als-consortium/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
nygc_spinal_cord_metadata = read_csv(here(shared_path, "nygc-als-consortium/sample-details/nygc_spinal_cord_metadata.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"nygc-als-consortium")) %>% mutate(condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(nygc.gene_tpm.male_counts, sample, gender))
nygc_postmortem_spinal_cord.tx = DESeq.analysis(metadata = nygc_spinal_cord_metadata, design = ~ gender + library_prep + sample_source + condition, transcript.level = TRUE, vsd_tx = FALSE, gene.level = FALSE, contrast = "condition_als_vs_ctrl")
saveRDS(nygc_postmortem_spinal_cord.tx, here(camp_path, "projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.tx.rds"))


# ### ML240  ----------------------------
# untreated_ml240.nuc_vs_cyt = DESeq.analysis(metadata = metadata, design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm")
# ml240.nuc_vs_cyt = DESeq.analysis(metadata = ml240.metadata, design = ~ cellline + fraction, transcript.level = TRUE, vsd_tx = TRUE, contrast = "fraction_nucleus_vs_cytoplasm")
# ml240.nuc_vs_cyt.tx = DESeq.analysis(metadata = ml240.metadata, design = ~ 1, transcript.level = TRUE, gene.level = FALSE, vsd_tx = TRUE)

# # Dual factor: ALS vs CTRL & nuc vs cyt
# ml240.nuc_vs_cyt.tardbp_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(ml240.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)))
# ml240.nuc_vs_cyt.tardbp_vs_ctrl.design_matrix = ml240.nuc_vs_cyt.tardbp_vs_ctrl.mm[, !apply(ml240.nuc_vs_cyt.tardbp_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# ml240.nuc_vs_cyt.tardbp_vs_ctrl = DESeq.analysis(metadata = mutate(filter(ml240.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)), 
#                                                      design = ml240.nuc_vs_cyt.tardbp_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                      contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                   "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationtardbp.fractionnucleus.IRFinderIR")
# ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(ml240.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)))
# ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl.design_matrix = ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm[, !apply(ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = mutate(filter(ml240.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)), 
#                                                         design = ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                         contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                      "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_r155c.fractionnucleus.IRFinderIR")
# ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(ml240.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)))
# ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl.design_matrix = ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm[, !apply(ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = mutate(filter(ml240.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), line_in_condition = recode_within(cellline, condition)), 
#                                                         design = ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                         contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                      "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_r191q.fractionnucleus.IRFinderIR")
# ml240.nuc_vs_cyt.vcp_iso_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(ml240.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl")), line_in_condition = recode_within(cellline, condition)))
# ml240.nuc_vs_cyt.vcp_iso_vs_ctrl.design_matrix = ml240.nuc_vs_cyt.vcp_iso_vs_ctrl.mm[, !apply(ml240.nuc_vs_cyt.vcp_iso_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# ml240.nuc_vs_cyt.vcp_iso_vs_ctrl = DESeq.analysis(metadata = mutate(filter(ml240.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl")), line_in_condition = recode_within(cellline, condition)), 
#                                                       design = ml240.nuc_vs_cyt.vcp_iso_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                       contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                    "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationvcp_iso.fractionnucleus.IRFinderIR")
# # ml240.nuc_vs_cyt.als_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(ml240.metadata, mutation %in% c("ctrl","vcp"), isogenic %in% c("isogenic","ctrl")), line_in_condition = recode_within(cellline, condition)))
# # ml240.nuc_vs_cyt.als_vs_ctrl.design_matrix = ml240.nuc_vs_cyt.als_vs_ctrl.mm[, !apply(ml240.nuc_vs_cyt.als_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# # ml240.nuc_vs_cyt.als_vs_ctrl = DESeq.analysis(metadata = mutate(ml240.metadata, line_in_condition = recode_within(cellline, condition)), design = ml240.nuc_vs_cyt.als_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
# #                                                   contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
# #                                                                "fractionnucleus.conditionals"=1)) # run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * IRFinder, irfinder_contrast_variable = c("mutation","fraction"), irfinder_contrast_name = "mutationals.fractionnucleus.IRFinderIR")
# 
# # Single factor: nuc vs cyt
# ctrl.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, mutation == "ctrl"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# # vcp.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, mutation == "vcp"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# # vcp_noiso.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, mutation == "vcp", isogenic == "isogenic"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# # vcp_iso.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, mutation == "vcp", isogenic == "no"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# # tardbp.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, mutation == "tardbp"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# # als.ml240.nuc_vs_cyt = DESeq.analysis(metadata = filter(ml240.metadata, condition == "als"), design = ~ cellline + fraction, contrast = "fraction_nucleus_vs_cytoplasm", transcript.level = TRUE)#, run_irfinder = TRUE, irfinder_design = ~ fraction + fraction:IRFinder, irfinder_contrast_variable = "fraction")
# 
# # Single factor: ALS vs CTRL
# # nuc.ml240.vcp_vs_ctrl = DESeq.analysis(metadata = nucleus.ml240.metadata, design = ~ mutation, contrast = "mutation_vcp_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation")
# # cyt.ml240.vcp_vs_ctrl = DESeq.analysis(metadata = cytoplasm.ml240.metadata, design = ~ mutation, contrast = "mutation_vcp_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation")
# # nuc.ml240.tardbp_vs_ctrl = DESeq.analysis(metadata = nucleus.ml240.metadata, design = ~ mutation, contrast = "mutation_tardbp_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation")
# # cyt.ml240.tardbp_vs_ctrl = DESeq.analysis(metadata = cytoplasm.ml240.metadata, design = ~ mutation, contrast = "mutation_tardbp_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation")
# # nuc.ml240.als_vs_ctrl = DESeq.analysis(metadata = nucleus.ml240.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# # cyt.ml240.als_vs_ctrl = DESeq.analysis(metadata = cytoplasm.ml240.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# 
# # Three factor: ml240 vs untx & ALS vs CTRL & nuc vs cyt
# # ml240_vs_ml240.nuc_vs_cyt.vcp_vs_ctrl = DESeq.analysis(metadata = filter(fraction.metadata, mutation %in% c("ctrl","vcp")), design = ~ fraction * mutation * treatment, contrast = "fractionnucleus.mutationvcp.treatmentml240", transcript.level = TRUE,
# #                                                            run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * treatment * IRFinder, irfinder_contrast_variable = c("mutation","fraction","treatment"), irfinder_contrast_name = "mutationvcp.fractionnucleus.treatmentml240.IRFinderIR")
# # ml240_vs_ml240.nuc_vs_cyt.tardbp_vs_ctrl = DESeq.analysis(metadata = filter(fraction.metadata, mutation %in% c("ctrl","tardbp")), design = ~ fraction * mutation * treatment, contrast = "fractionnucleus.mutationtardbp.treatmentml240", transcript.level = TRUE,
# #                                                            run_irfinder = TRUE, irfinder_design = ~ mutation * fraction * treatment * IRFinder, irfinder_contrast_variable = c("mutation","fraction","treatment"), irfinder_contrast_name = "mutationtardbp.fractionnucleus.treatmentml240.IRFinderIR")
# # ml240_vs_ml240.nuc_vs_cyt.als_vs_ctrl = DESeq.analysis(metadata = filter(fraction.metadata, condition %in% c("ctrl","als")), design = ~ fraction * condition * treatment, contrast = "fractionnucleus.conditionals.treatmentml240", transcript.level = TRUE,
# #                                                               run_irfinder = TRUE, irfinder_design = ~ condition * fraction * treatment * IRFinder, irfinder_contrast_variable = c("condition","fraction","treatment"), irfinder_contrast_name = "conditionals.fractionnucleus.treatmentml240.IRFinderIR")
# 
# print("saving ml240 treated")
# save(ml240.nuc_vs_cyt, ml240.nuc_vs_cyt.tx, untreated_ml240.nuc_vs_cyt, ctrl.ml240.nuc_vs_cyt, #vcp_iso.ml240.nuc_vs_cyt, tardbp.ml240.nuc_vs_cyt, ml240.nuc_vs_cyt.vcp_vs_ctrl, ml240.nuc_vs_cyt.vcp_noiso_vs_ctrl,
#      # nuc.ml240.vcp_vs_ctrl, cyt.ml240.vcp_vs_ctrl, nuc.ml240.tardbp_vs_ctrl, cyt.ml240.tardbp_vs_ctrl, #nuc.ml240.als_vs_ctrl, cyt.ml240.als_vs_ctrl, #whole.ml240.als_vs_ctrl,
#      ml240.nuc_vs_cyt.vcp_r191q_vs_ctrl, ml240.nuc_vs_cyt.vcp_r155c_vs_ctrl, ml240.nuc_vs_cyt.vcp_iso_vs_ctrl, ml240.nuc_vs_cyt.tardbp_vs_ctrl, #ml240.nuc_vs_cyt.als_vs_ctrl, #als.ml240.nuc_vs_cyt,
#      # ml240_vs_ml240.nuc_vs_cyt.vcp_vs_ctrl, ml240_vs_ml240.nuc_vs_cyt.tardbp_vs_ctrl, ml240_vs_ml240.nuc_vs_cyt.als_vs_ctrl,
#      file = here(proj_path, "scripts/ml240_sample_deseq_objects.RData"))
# print("ml240 treated saved")
# 
# # Whole cell ALS vs CTRL ----------------------------
# whole.untreated.tardbp_vs_ctrl = DESeq.analysis(metadata = filter(whole.untreated.metadata, mutation %in% c("ctrl","tardbp")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.untreated.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = filter(whole.untreated.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.untreated.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = filter(whole.untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.untreated.vcp_iso_vs_ctrl = DESeq.analysis(metadata = filter(whole.untreated.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("isogenic","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# print("saving untreated whole")
# save(whole.untreated.tardbp_vs_ctrl, whole.untreated.vcp_r155c_vs_ctrl, whole.untreated.vcp_r191q_vs_ctrl, whole.untreated.vcp_iso_vs_ctrl,
#      file = here(proj_path, "scripts/whole_untreated_deseq_objects.RData"))
# print("untreated whole saved")
# 
# whole.ml240.tardbp_vs_ctrl = DESeq.analysis(metadata = filter(whole.ml240.metadata, mutation %in% c("ctrl","tardbp")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.ml240.vcp_r155c_vs_ctrl = DESeq.analysis(metadata = filter(whole.ml240.metadata, genotype %in% c("ctrl","r155c"), isogenic %in% c("no","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.ml240.vcp_r191q_vs_ctrl = DESeq.analysis(metadata = filter(whole.ml240.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("no","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# whole.ml240.vcp_iso_vs_ctrl = DESeq.analysis(metadata = filter(whole.ml240.metadata, genotype %in% c("ctrl","r191q"), isogenic %in% c("isogenic","ctrl")), design = ~ condition, contrast = "condition_als_vs_ctrl", transcript.level = TRUE, run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder, irfinder_contrast_variable = "condition")
# print("saving ml240 whole")
# save(whole.ml240.tardbp_vs_ctrl, whole.ml240.vcp_r155c_vs_ctrl, whole.ml240.vcp_r191q_vs_ctrl, whole.ml240.vcp_iso_vs_ctrl,
#      file = here(proj_path, "scripts/whole_ml240_sample_deseq_objects.RData"))
# print("ml240 whole saved")

# ML240 vs Untreated
whole.ctrl.ml240_vs_untreated = DESeq.analysis(metadata = filter(whole.metadata, mutation %in% c("ctrl")), design = ~ cellline + treatment, contrast = "treatment_ml240_vs_untreated", transcript.level = FALSE)
whole.tardbp.ml240_vs_untreated = DESeq.analysis(metadata = filter(whole.metadata, mutation %in% c("tardbp")), design = ~ cellline + treatment, contrast = "treatment_ml240_vs_untreated", transcript.level = FALSE)
whole.vcp_r155c.ml240_vs_untreated = DESeq.analysis(metadata = filter(whole.metadata, genotype %in% c("r155c"), isogenic %in% c("no")), design = ~ cellline + treatment, contrast = "treatment_ml240_vs_untreated", transcript.level = FALSE)
whole.vcp_r191q.ml240_vs_untreated = DESeq.analysis(metadata = filter(whole.metadata, genotype %in% c("r191q"), isogenic %in% c("no")), design = ~ cellline + treatment, contrast = "treatment_ml240_vs_untreated", transcript.level = FALSE)
whole.vcp_iso.ml240_vs_untreated = DESeq.analysis(metadata = filter(whole.metadata, genotype %in% c("r191q"), isogenic %in% c("isogenic")), design = ~ cellline + treatment, contrast = "treatment_ml240_vs_untreated", transcript.level = FALSE)
save(whole.ctrl.ml240_vs_untreated, whole.tardbp.ml240_vs_untreated, whole.vcp_r155c.ml240_vs_untreated, whole.vcp_r191q.ml240_vs_untreated, whole.vcp_iso.ml240_vs_untreated,
     file = here(proj_path, "scripts/whole_ml240_vs_untreated_deseq_objects.RData"))

# # ### Mass Spec ----------------------------
# 
# experimental design
sample_submission = read_excel(here(proj_path, "mass-spec/Nuc_Cyto_Sample_submission_5.4.22.xlsx"), skip = 1) %>% clean_names() %>% mutate(sample_name = gsub(" ","",tolower(sample_decoded)), cellline = str_split_fixed(sample_name,"_",3)[,1], treatment = str_split_fixed(sample_name,"_",3)[,2], fraction = str_split_fixed(sample_name,"_",3)[,3], mutation = case_when(mutation == "Wildtype" ~ "ctrl", TRUE ~ tolower(mutation))) %>%
  filter(!sample_name %in% c("ctrl6_ml240_nuclear", "glib_ml240_nuclear", "tdpn1_ml240_nuclear", "tdpn1_ml240_cytoplasm", "ctrl3_untreated_nuclear")) %>% # remove samples that failed qc
  select(sample_number, sample_name, cellline, treatment, fraction, mutation)
sample_replicates = bind_rows(mutate(sample_submission, technical_replicate = 1), mutate(sample_submission, technical_replicate = 2), mutate(sample_submission, technical_replicate = 3)) %>%
  mutate(lfq_id = paste0("lfq_intensity_", sample_number, "_", technical_replicate), lfq_name = gsub(" ", "", paste0("lfq_intensity_", sample_name, "_", technical_replicate)))

# MaxQuant LFQ intentisities
lfq_intensities = read_excel(here(proj_path, "mass-spec/IB2652_02092022.xlsx")) %>% clean_names() %>%
  select(!starts_with("lfq_intensity_2_")) %>% # remove untreated ctrl3 nuclear
  filter(reverse != "+", potential_contaminant != "+") %>%
  rename_at(vars(starts_with("lfq_intensity")), ~ str_replace_all(., setNames(sample_replicates$lfq_name, sample_replicates$lfq_id)))

lfq_intensities$protein_i_ds %>% duplicated() %>% any()
lfq_intensities %>% group_by(protein_i_ds) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
lfq_intensities.unique <- make_unique(lfq_intensities, "gene_names", "protein_i_ds", delim = ";")
lfq_intensities.unique$name %>% duplicated() %>% any()
lfq_intensities.unique %>% glimpse

sample_replicates$lfq_name[sample_replicates$lfq_name %in% colnames(lfq_intensities.unique)]
sample_replicates$lfq_name[!sample_replicates$lfq_name %in% colnames(lfq_intensities.unique)]

### Untreated ----------------------------

# all untreated samples for PCA
sample_replicates.untreated.nuc_vs_cyt = sample_replicates %>% filter(treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.untreated.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240"))
untreated.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.nuc_vs_cyt, experimental_design = sample_replicates.untreated.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                        design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.nuc_vs_cyt")

# nuc_vs_cyt
sample_replicates.untreated.ctrl_nuc_vs_cyt = sample_replicates %>% filter(mutation == "ctrl", treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.untreated.ctrl_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm|tdpn"))
untreated.ctrl_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.ctrl_nuc_vs_cyt, experimental_design = sample_replicates.untreated.ctrl_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                             design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.ctrl_nuc_vs_cyt")
# sample_replicates.untreated.tardbp_nuc_vs_cyt = sample_replicates %>% filter(mutation == "tardbp", treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.untreated.tardbp_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm|ctrl"))
# untreated.tardbp_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.tardbp_nuc_vs_cyt, experimental_design = sample_replicates.untreated.tardbp_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
#                                              design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.tardbp_nuc_vs_cyt")
# sample_replicates.untreated.vcp_r155c_nuc_vs_cyt = sample_replicates %>% filter(grepl("cb1",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.untreated.vcp_r155c_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|tdpn|ncrm|ctrl"))
# untreated.vcp_r155c_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r155c_nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_r155c_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
#                                                design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_r155c_nuc_vs_cyt")
# sample_replicates.untreated.vcp_r191q_nuc_vs_cyt = sample_replicates %>% filter(grepl("gli",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.untreated.vcp_r191q_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|cb1|tdpn|ncrm|ctrl"))
# untreated.vcp_r191q_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r191q_nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_r191q_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
#                                                   design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_r191q_nuc_vs_cyt")
# sample_replicates.untreated.vcp_iso_nuc_vs_cyt = sample_replicates %>% filter(grepl("ncrm",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.untreated.vcp_iso_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|tdpn|cb1|ctrl"))
# untreated.vcp_iso_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_iso_nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_iso_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
#                                                 design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_iso_nuc_vs_cyt")
#
# #als_vs_ctrl.nuc_vs_cyt
# sample_replicates.untreated.tardbp_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.tardbp","nuclear.tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.untreated.tardbp_vs_ctrl.nuc_vs_cyt)
# lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm"))
# untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.untreated.tardbp_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#                                                        model_matrix = untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.tardbp - cytoplasm.tardbp) - (nuclear.ctrl - cytoplasm.ctrl)",
#                                                        file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.tardbp_vs_ctrl.nuc_vs_cyt")

# sample_replicates.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|cb1", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt)
# lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|gli"))
# untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#                                                        model_matrix = untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
#                                                        file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt")

# sample_replicates.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|gli", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt)
# lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|cb1"))
# untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#                                                        model_matrix = untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
#                                                        file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt")

# sample_replicates.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|ncrm", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt)
# lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|gli|cb1"))
# untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#                                                        model_matrix = untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
#                                                        file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.vcp_iso_vs_ctrl.nuc_vs_cyt")
#
# print("saving untreated mass spec")
# save(untreated.nuc_vs_cyt.dep, untreated.ctrl_nuc_vs_cyt.dep, #untreated.tardbp_nuc_vs_cyt.dep, untreated.vcp_noiso_nuc_vs_cyt.dep, untreated.vcp_iso_nuc_vs_cyt.dep,
#      untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep, untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep, untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep, untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep,
#      file = here(proj_path, "scripts/dep_mass_spec_objects.RData"))
# print("untreated mass spec saved")

# sample_replicates.untreated.als_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(als = case_when(mutation == "ctrl"~"ctrl",TRUE~"als"), condition = factor(paste(fraction,als,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.als","nuclear.als"))) %>% #, line_in_condition = recode_within(cellline, condition)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)#, line_in_condition)
# untreated.als_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + cellline, data = sample_replicates.untreated.als_vs_ctrl.nuc_vs_cyt)
# lfq_intensities.unique.untreated.als_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("ml240"))
# untreated.als_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.als_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.untreated.als_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#                                                        model_matrix = untreated.als_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.als - cytoplasm.als) - (nuclear.ctrl - cytoplasm.ctrl)",
#                                                        file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.als_vs_ctrl.nuc_vs_cyt")
# saveRDS(untreated.als_vs_ctrl.nuc_vs_cyt.dep, file = here(proj_path, "scripts/untreated.als_vs_ctrl.nuc_vs_cyt.dep.rds"))


# #als_vs_ctrl.nuc_cyt
# sample_replicates.untreated.tardbp_vs_ctrl.nuc_cyt = sample_replicates %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(mutation, levels = c("ctrl","tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_cyt = lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm"))
# untreated.tardbp_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_cyt, experimental_design = sample_replicates.untreated.tardbp_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "tardbp",
#                                         design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.nuc_cyt.tardbp_vs_ctrl")
# 
# sample_replicates.untreated.vcp_r155c_vs_ctrl.nuc_cyt = sample_replicates %>% filter(grepl("ctrl|cb1", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|gli"))
# untreated.vcp_r155c_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_cyt, experimental_design = sample_replicates.untreated.vcp_r155c_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
#                                                     design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.nuc_cyt.tardbp_vs_ctrl")
# 
# sample_replicates.untreated.vcp_r191q_vs_ctrl.nuc_cyt = sample_replicates %>% filter(grepl("ctrl|gli", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|cb1"))
# untreated.vcp_r191q_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_cyt, experimental_design = sample_replicates.untreated.vcp_r191q_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
#                                                     design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.nuc_cyt.vcp_r191q_vs_ctrl")
# 
# sample_replicates.untreated.vcp_iso_vs_ctrl.nuc_cyt = sample_replicates %>% filter(grepl("ctrl|ncrm", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_cyt = lfq_intensities.unique %>% select(!matches("ml240|tdpn|gli|cb1"))
# untreated.vcp_iso_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_cyt, experimental_design = sample_replicates.untreated.vcp_iso_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
#                                                     design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated.nuc_cyt.vcp_iso_vs_ctrl")
# 
# save(untreated.tardbp_vs_ctrl.nuc_cyt.dep, untreated.vcp_r155c_vs_ctrl.nuc_cyt.dep, untreated.vcp_r191q_vs_ctrl.nuc_cyt.dep, untreated.vcp_iso_vs_ctrl.nuc_cyt.dep,
#      file = here(proj_path, "scripts/dep_mass_spec.untreated.als_vs_ctrl.nuc_cyt.RData"))
#
#
# ### ML240 ----------------------------
# 
# # all untreated & ml240 samples for PCA
# sample_replicates.untreated_ml240.nuc_vs_cyt = sample_replicates %>% filter(treatment %in% c("untreated", "ml240")) %>% group_by(treatment, cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# untreated_ml240.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique, experimental_design = sample_replicates.untreated_ml240.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
#                                     design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "untreated_ml240.nuc_vs_cyt")
# 
# # all ml240 samples for PCA
# sample_replicates.ml240.nuc_vs_cyt = sample_replicates %>% filter(treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.ml240.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated"))
# ml240.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.nuc_vs_cyt, experimental_design = sample_replicates.ml240.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
#                                     design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.nuc_vs_cyt")
# 
# 
# # nuc_vs_cyt
# sample_replicates.ml240.ctrl_nuc_vs_cyt = sample_replicates %>% filter(mutation == "ctrl", treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
# lfq_intensities.unique.ml240.ctrl_nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated|gli|cb1|ncrm|tdpn"))
# ml240.ctrl_nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.ctrl_nuc_vs_cyt, experimental_design = sample_replicates.ml240.ctrl_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
#                                          design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.ctrl_nuc_vs_cyt")
# 
# 
# #als_vs_ctrl.nuc_vs_cyt
# sample_replicates.ml240.tardbp_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.tardbp","nuclear.tardbp")), line_in_condition = recode_within(cellline, mutation)) %>% 
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.ml240.tardbp_vs_ctrl.nuc_vs_cyt) 
# lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated|gli|cb1|ncrm"))
# ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.ml240.tardbp_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#            model_matrix = ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.tardbp - cytoplasm.tardbp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.tardbp_vs_ctrl.nuc_vs_cyt")
# sample_replicates.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|cb1", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>% 
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt) 
# lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|gli"))
# ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#             model_matrix = ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt")
# sample_replicates.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|gli", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>% 
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt) 
# lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|cb1"))
# ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#             model_matrix = ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt")
# sample_replicates.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt = sample_replicates %>% filter(grepl("ctrl|ncrm", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
#   mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>% 
#   select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
# ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = sample_replicates.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt) 
# lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt = lfq_intensities.unique %>% select(!matches("untreated|tdpn|gli|cb1"))
# ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt, experimental_design = sample_replicates.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
#               model_matrix = ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "ml240.vcp_iso_vs_ctrl.nuc_vs_cyt")
# 
# print("saving ml240 mass spec")
# save(untreated_ml240.nuc_vs_cyt.dep, ml240.nuc_vs_cyt.dep, ml240.ctrl_nuc_vs_cyt.dep, #ml240.tardbp_nuc_vs_cyt.dep, ml240.vcp_noiso_nuc_vs_cyt.dep, ml240.vcp_iso_nuc_vs_cyt.dep,
#      ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep, ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep, ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep, ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep,
#      file = here(proj_path, "scripts/dep_mass_spec_ml240_objects.RData"))
# print("ml240 mass spec saved")

# nuc_cyt.ml240_vs_untreated
sample_replicates.ctrl.nuc_cyt.ml240_vs_untreated = sample_replicates %>% filter(mutation %in% c("ctrl")) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>% mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.ctrl.nuc_cyt.ml240_vs_untreated = lfq_intensities.unique %>% select(!matches("tdpn|gli|cb1|ncrm"))
ctrl.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = lfq_intensities.unique.ctrl.nuc_cyt.ml240_vs_untreated, experimental_design = sample_replicates.ctrl.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "ctrl.nuc_cyt.ml240_vs_untreated")

sample_replicates.tardbp.nuc_cyt.ml240_vs_untreated = sample_replicates %>% filter(mutation %in% c("tardbp")) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>% mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.tardbp.nuc_cyt.ml240_vs_untreated = lfq_intensities.unique %>% select(!matches("ctrl|gli|cb1|ncrm"))
tardbp.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = lfq_intensities.unique.tardbp.nuc_cyt.ml240_vs_untreated, experimental_design = sample_replicates.tardbp.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                     design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "tardbp.nuc_cyt.ml240_vs_untreated")

sample_replicates.vcp_r155c.nuc_cyt.ml240_vs_untreated = sample_replicates %>% filter(grepl("cb1", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>% mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.vcp_r155c.nuc_cyt.ml240_vs_untreated = lfq_intensities.unique %>% select(!matches("ctrl|tdpn|ncrm|gli"))
vcp_r155c.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = lfq_intensities.unique.vcp_r155c.nuc_cyt.ml240_vs_untreated, experimental_design = sample_replicates.vcp_r155c.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "vcp_r155c.nuc_cyt.ml240_vs_untreated")

sample_replicates.vcp_r191q.nuc_cyt.ml240_vs_untreated = sample_replicates %>% filter(grepl("gli", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>% mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.vcp_r191q.nuc_cyt.ml240_vs_untreated = lfq_intensities.unique %>% select(!matches("ctrl|tdpn|ncrm|cb1"))
vcp_r191q.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = lfq_intensities.unique.vcp_r191q.nuc_cyt.ml240_vs_untreated, experimental_design = sample_replicates.vcp_r191q.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "vcp_r191q.nuc_cyt.ml240_vs_untreated")

sample_replicates.vcp_iso.nuc_cyt.ml240_vs_untreated = sample_replicates %>% filter(grepl("ncrm", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>% mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
lfq_intensities.unique.vcp_iso.nuc_cyt.ml240_vs_untreated = lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|cb1"))
vcp_iso.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = lfq_intensities.unique.vcp_iso.nuc_cyt.ml240_vs_untreated, experimental_design = sample_replicates.vcp_iso.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "vcp_iso.nuc_cyt.ml240_vs_untreated")

save(ctrl.nuc_cyt.ml240_vs_untreated.dep, tardbp.nuc_cyt.ml240_vs_untreated.dep, vcp_r155c.nuc_cyt.ml240_vs_untreated.dep, vcp_r191q.nuc_cyt.ml240_vs_untreated.dep, vcp_iso.nuc_cyt.ml240_vs_untreated.dep,
     file = here(proj_path, "scripts/nuc_cyt.ml240_vs_untreated.dep.RData"))

# # Rebuttal Reviewer 1, Comment 4
# nuc_who.untreated.metadata = read_csv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   filter(fraction != "cytoplasm", sample != "ctrl3_untreated_nucleus") %>% #failed nucleus QC - clusters with whole so remove
#   mutate(database_dir = here(collab_path,"motor-neuron-vcp-inhibitor-harley-2022"), mutation = case_when(cellline %in% c("ctrl1", "ctrl3", "ctrl4", "ctrl6") ~ "ctrl", cellline %in% c("glia", "glib", "cb1d", "cb1e", "ncrme6", "ncrmc2") ~ "vcp", cellline %in% c("tdpn1", "tdpn2") ~ "tardbp"), 
#          genotype = case_when(cellline %in% c("cb1d", "cb1e") ~ "r155c", cellline %in% c("glia", "glib", "ncrme6", "ncrmc2") ~ "r191q", cellline %in% c("tdpn1", "tdpn2") ~ "g298s", TRUE ~ "ctrl"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), 
#          isogenic = case_when(cellline %in% c("ncrme6", "ncrmc2") ~ "isogenic", mutation == "ctrl" ~ "ctrl", TRUE ~ "no"), parent_cellline = case_when(cellline %in% c("ncrmc2", "ncrme6") ~ "ctrl6", TRUE ~ cellline),
#          treatment = factor(treatment, levels = c("untreated", "ml240")), fraction = gsub("nuclear","nucleus",fraction), fraction = factor(fraction, levels = c("whole", "nucleus"))) %>%
#   select(sample, fraction, treatment, condition, mutation, cellline, genotype, isogenic, parent_cellline, database_dir) %>% filter(treatment == "untreated")
# untreated.nuc_vs_who.tardbp_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(nuc_who.untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)))
# untreated.nuc_vs_who.tardbp_vs_ctrl.design_matrix = untreated.nuc_vs_who.tardbp_vs_ctrl.mm[, !apply(untreated.nuc_vs_who.tardbp_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# untreated.nuc_vs_who.tardbp_vs_ctrl = DESeq.analysis(metadata = mutate(filter(nuc_who.untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)), 
#                                                      design = untreated.nuc_vs_who.tardbp_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                      contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                   "fractionnucleus.conditionals"=1)) 
# untreated.nuc_vs_who.tardbp_vs_ctrl$res.tx %>% filter(transcript_id %in% "ENST00000619417")
# 
# cyt_who.untreated.metadata = read_csv(here(collab_path, "motor-neuron-vcp-inhibitor-harley-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   filter(fraction != "nucleus", sample != "ctrl3_untreated_nucleus") %>% #failed nucleus QC - clusters with whole so remove
#   mutate(database_dir = here(collab_path,"motor-neuron-vcp-inhibitor-harley-2022"), mutation = case_when(cellline %in% c("ctrl1", "ctrl3", "ctrl4", "ctrl6") ~ "ctrl", cellline %in% c("glia", "glib", "cb1d", "cb1e", "ncrme6", "ncrmc2") ~ "vcp", cellline %in% c("tdpn1", "tdpn2") ~ "tardbp"), 
#          genotype = case_when(cellline %in% c("cb1d", "cb1e") ~ "r155c", cellline %in% c("glia", "glib", "ncrme6", "ncrmc2") ~ "r191q", cellline %in% c("tdpn1", "tdpn2") ~ "g298s", TRUE ~ "ctrl"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"), condition = factor(condition, levels = c("ctrl", "als")), 
#          isogenic = case_when(cellline %in% c("ncrme6", "ncrmc2") ~ "isogenic", mutation == "ctrl" ~ "ctrl", TRUE ~ "no"), parent_cellline = case_when(cellline %in% c("ncrmc2", "ncrme6") ~ "ctrl6", TRUE ~ cellline),
#          treatment = factor(treatment, levels = c("untreated", "ml240")), fraction = gsub("nuclear","nucleus",fraction), fraction = factor(fraction, levels = c("whole", "cytoplasm"))) %>%
#   select(sample, fraction, treatment, condition, mutation, cellline, genotype, isogenic, parent_cellline, database_dir) %>% filter(treatment == "untreated")
# untreated.cyt_vs_who.tardbp_vs_ctrl.mm <- model.matrix( ~ fraction + condition + condition:line_in_condition + fraction:condition, data=mutate(filter(cyt_who.untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)))
# untreated.cyt_vs_who.tardbp_vs_ctrl.design_matrix = untreated.cyt_vs_who.tardbp_vs_ctrl.mm[, !apply(untreated.cyt_vs_who.tardbp_vs_ctrl.mm==0, 2, all)]  # remove missing covariates from design matrix
# untreated.cyt_vs_who.tardbp_vs_ctrl = DESeq.analysis(metadata = mutate(filter(cyt_who.untreated.metadata, mutation %in% c("ctrl","tardbp")), line_in_condition = recode_within(cellline, condition)), 
#                                                      design = untreated.cyt_vs_who.tardbp_vs_ctrl.design_matrix, transcript.level = TRUE, combine.contrasts = TRUE,
#                                                      contrast = c("Intercept"= 0,"fractionnucleus"= 0,"conditionals"= 0, "conditionctrl.line_in_condition2"= 0, "conditionals.line_in_condition2"=0, "conditionctrl.line_in_condition3"=0 ,"conditionctrl.line_in_condition4"=0, 
#                                                                   "fractionnucleus.conditionals"=1)) 
# untreated.cyt_vs_who.tardbp_vs_ctrl$res.tx %>% filter(transcript_id %in% "ENST00000619417") 


# Co-immunoprecipitation VCP --------------------------------------------------------
# experimental design
coip_samples = read_excel(here(proj_path, "mass-spec/vcp_co_ip/fairouz_lfq_ib2608.xlsx"), sheet = "annotation") %>% clean_names() %>% 
  mutate(cellline = tolower(str_split_fixed(sample_name," ",3)[,1]), treatment = str_split_fixed(sample_name," ",3)[,2], treatment = case_when(treatment == "UT"~"untreated",treatment == "ML240"~"ml240"), fraction_replicate = str_split_fixed(sample_name," ",3)[,3], fraction = str_split_fixed(fraction_replicate,"_",2)[,1], 
         fraction = case_when(fraction == "N"~"nuclear", fraction=="C"~"cytoplasm"), technical_replicate = str_split_fixed(fraction_replicate,"_",2)[,2], 
         mutation = case_when(grepl("ctrl",cellline) ~ "ctrl", grepl("tdpn",cellline) ~ "tardbp", TRUE ~ "vcp"), sample_number = str_split_fixed(lfq_id," ",3)[,3], sample_number = str_split_fixed(sample_number,"_",2)[,1], sample_name = paste(cellline, treatment, fraction, sep="_"), 
         lfq_id = tolower(gsub(" ","_",lfq_id)), lfq_name = paste0("lfq_intensity_", sample_name, "_", technical_replicate)) %>%
  # filter(!sample_name %in% c("ctrl6_ml240_nuclear", "glib_ml240_nuclear", "tdpn1_ml240_nuclear", "tdpn1_ml240_cytoplasm", "ctrl3_untreated_nuclear")) %>% # remove samples that failed qc
  select(sample_number, sample_name, cellline, treatment, fraction, mutation, technical_replicate, lfq_id, lfq_name)

# MaxQuant LFQ intentisities
coip_lfq_intensities = read_excel(here(proj_path, "mass-spec/vcp_co_ip/fairouz_lfq_ib2608.xlsx"), sheet = "LFQ") %>% clean_names() %>%
  # select(!starts_with("lfq_intensity_2_")) %>% # remove untreated ctrl3 nuclear
  filter(is.na(reverse), is.na(potential_contaminant)) %>%
  rename_at(vars(starts_with("lfq_intensity")), ~ str_replace_all(., setNames(coip_samples$lfq_name, coip_samples$lfq_id)))

coip_lfq_intensities$protein_i_ds %>% duplicated() %>% any()
coip_lfq_intensities %>% group_by(protein_i_ds) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
coip_lfq_intensities.unique <- make_unique(coip_lfq_intensities, "gene_names", "protein_i_ds", delim = ";")
coip_lfq_intensities.unique$name %>% duplicated() %>% any()
coip_lfq_intensities.unique %>% glimpse

coip_samples$lfq_name[coip_samples$lfq_name %in% colnames(coip_lfq_intensities.unique)]
coip_samples$lfq_name[!coip_samples$lfq_name %in% colnames(coip_lfq_intensities.unique)]

### Untreated Co-IP ----------------------------

# all untreated samples for PCA
coip_samples.untreated.nuc_vs_cyt = coip_samples %>% filter(treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240"))
coip_untreated.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.nuc_vs_cyt, experimental_design = coip_samples.untreated.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                        design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.nuc_vs_cyt")

# nuc_vs_cyt
coip_samples.untreated.ctrl_nuc_vs_cyt = coip_samples %>% filter(mutation == "ctrl", treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.ctrl_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm|tdpn"))
coip_untreated.ctrl_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.ctrl_nuc_vs_cyt, experimental_design = coip_samples.untreated.ctrl_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                             design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.ctrl_nuc_vs_cyt")
coip_samples.untreated.tardbp_nuc_vs_cyt = coip_samples %>% filter(mutation == "tardbp", treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.tardbp_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm|ctrl"))
coip_untreated.tardbp_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.tardbp_nuc_vs_cyt, experimental_design = coip_samples.untreated.tardbp_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
                                             design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.tardbp_nuc_vs_cyt")
coip_samples.untreated.vcp_r155c_nuc_vs_cyt = coip_samples %>% filter(grepl("cb1",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.vcp_r155c_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|tdpn|ncrm|ctrl"))
coip_untreated.vcp_r155c_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r155c_nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_r155c_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
                                               design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_r155c_nuc_vs_cyt")
coip_samples.untreated.vcp_r191q_nuc_vs_cyt = coip_samples %>% filter(grepl("gli",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.vcp_r191q_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|cb1|tdpn|ncrm|ctrl"))
coip_untreated.vcp_r191q_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r191q_nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_r191q_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
                                                  design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_r191q_nuc_vs_cyt")
coip_samples.untreated.vcp_iso_nuc_vs_cyt = coip_samples %>% filter(grepl("ncrm",cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.untreated.vcp_iso_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|tdpn|cb1|ctrl"))
coip_untreated.vcp_iso_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_iso_nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_iso_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "cytoplasm", numerator = "nuclear",
                                                design_formula = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_iso_nuc_vs_cyt")

#als_vs_ctrl.nuc_vs_cyt
coip_samples.untreated.tardbp_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.tardbp","nuclear.tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.untreated.tardbp_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm"))
coip_untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.untreated.tardbp_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
                                                       model_matrix = untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.tardbp - cytoplasm.tardbp) - (nuclear.ctrl - cytoplasm.ctrl)",
                                                       file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.tardbp_vs_ctrl.nuc_vs_cyt")

coip_samples.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|cb1", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|gli"))
coip_untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
                                                       model_matrix = untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
                                                       file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt")

coip_samples.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|gli", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|cb1"))
coip_untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
                                                       model_matrix = untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
                                                       file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt")

coip_samples.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|ncrm", cellline), treatment == "untreated") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|gli|cb1"))
coip_untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.untreated.vcp_iso_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
                                                       model_matrix = untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)",
                                                       file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.vcp_iso_vs_ctrl.nuc_vs_cyt")

print("saving untreated mass spec")
save(coip_untreated.nuc_vs_cyt.dep, coip_untreated.ctrl_nuc_vs_cyt.dep, #coip_untreated.tardbp_nuc_vs_cyt.dep, coip_untreated.vcp_noiso_nuc_vs_cyt.dep, coip_untreated.vcp_iso_nuc_vs_cyt.dep,
     coip_untreated.tardbp_vs_ctrl.nuc_vs_cyt.dep, coip_untreated.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep, coip_untreated.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep, coip_untreated.vcp_iso_vs_ctrl.nuc_vs_cyt.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec_objects.RData"))
print("untreated mass spec saved")


#als_vs_ctrl.nuc_cyt
coip_samples.untreated.tardbp_vs_ctrl.nuc_cyt = coip_samples %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|gli|cb1|ncrm"))
coip_untreated.tardbp_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.tardbp_vs_ctrl.nuc_cyt, experimental_design = coip_samples.untreated.tardbp_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "tardbp",
                                        design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.nuc_cyt.tardbp_vs_ctrl")

coip_samples.untreated.vcp_r155c_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|cb1", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|gli"))
coip_untreated.vcp_r155c_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r155c_vs_ctrl.nuc_cyt, experimental_design = coip_samples.untreated.vcp_r155c_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                    design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.nuc_cyt.tardbp_vs_ctrl")

coip_samples.untreated.vcp_r191q_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|gli", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|ncrm|cb1"))
coip_untreated.vcp_r191q_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_r191q_vs_ctrl.nuc_cyt, experimental_design = coip_samples.untreated.vcp_r191q_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                    design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.nuc_cyt.vcp_r191q_vs_ctrl")

coip_samples.untreated.vcp_iso_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|ncrm", cellline), treatment == "untreated") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("ml240|tdpn|gli|cb1"))
coip_untreated.vcp_iso_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.untreated.vcp_iso_vs_ctrl.nuc_cyt, experimental_design = coip_samples.untreated.vcp_iso_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                    design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated.nuc_cyt.vcp_iso_vs_ctrl")

save(coip_untreated.tardbp_vs_ctrl.nuc_cyt.dep, coip_untreated.vcp_r155c_vs_ctrl.nuc_cyt.dep, coip_untreated.vcp_r191q_vs_ctrl.nuc_cyt.dep, coip_untreated.vcp_iso_vs_ctrl.nuc_cyt.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec.untreated.als_vs_ctrl.nuc_cyt.RData"))


### ML240 Co-IP ----------------------------

# all untreated & ml240 samples for PCA
coip_samples.untreated_ml240.nuc_vs_cyt = coip_samples %>% filter(treatment %in% c("untreated", "ml240")) %>% group_by(treatment, cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_untreated_ml240.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique, experimental_design = coip_samples.untreated_ml240.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                    design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_untreated_ml240.nuc_vs_cyt")

# all ml240 samples for PCA
coip_samples.ml240.nuc_vs_cyt = coip_samples %>% filter(treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.ml240.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated"))
coip_ml240.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.nuc_vs_cyt, experimental_design = coip_samples.ml240.nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                    design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.nuc_vs_cyt")


# nuc_vs_cyt
coip_samples.ml240.ctrl_nuc_vs_cyt = coip_samples %>% filter(mutation == "ctrl", treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  select(label = lfq_name, condition = fraction, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.ml240.ctrl_nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|gli|cb1|ncrm|tdpn"))
coip_ml240.ctrl_nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.ctrl_nuc_vs_cyt, experimental_design = coip_samples.ml240.ctrl_nuc_vs_cyt, show_proteins_per_threshold = TRUE, threshold = 10, denomenator = "cytoplasm", numerator = "nuclear",
                                         design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.ctrl_nuc_vs_cyt")


#als_vs_ctrl.nuc_vs_cyt
coip_samples.ml240.tardbp_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.tardbp","nuclear.tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.ml240.tardbp_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|gli|cb1|ncrm"))
coip_ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.ml240.tardbp_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
           model_matrix = ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.tardbp - cytoplasm.tardbp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.tardbp_vs_ctrl.nuc_vs_cyt")
coip_samples.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|cb1", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|gli"))
coip_ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
            model_matrix = ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt")
coip_samples.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|gli", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|cb1"))
coip_ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
            model_matrix = ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt")
coip_samples.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt = coip_samples %>% filter(grepl("ctrl|ncrm", cellline), treatment == "ml240") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(paste(fraction,mutation,sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.vcp","nuclear.vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm = model.matrix(~0 + condition + technical_replicate + line_in_condition, data = coip_samples.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt)
coip_lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|gli|cb1"))
coip_ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt, experimental_design = coip_samples.ml240.vcp_iso_vs_ctrl.nuc_vs_cyt, show_proteins_per_threshold = FALSE, threshold = 5,
              model_matrix = ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep.mm, contrast = "(nuclear.vcp - cytoplasm.vcp) - (nuclear.ctrl - cytoplasm.ctrl)", file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.vcp_iso_vs_ctrl.nuc_vs_cyt")

print("saving ml240 mass spec")
save(coip_untreated_ml240.nuc_vs_cyt.dep, coip_ml240.nuc_vs_cyt.dep, coip_ml240.ctrl_nuc_vs_cyt.dep, #coip_ml240.tardbp_nuc_vs_cyt.dep, coip_ml240.vcp_noiso_nuc_vs_cyt.dep, coip_ml240.vcp_iso_nuc_vs_cyt.dep,
     coip_ml240.tardbp_vs_ctrl.nuc_vs_cyt.dep, coip_ml240.vcp_r191q_vs_ctrl.nuc_vs_cyt.dep, coip_ml240.vcp_r155c_vs_ctrl.nuc_vs_cyt.dep, coip_ml240.vcp_iso_vs_ctrl.nuc_vs_cyt.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec_ml240_objects.RData"))
print("ml240 mass spec saved")

#als_vs_ctrl.nuc_cyt
coip_samples.ml240.tardbp_vs_ctrl.nuc_cyt = coip_samples %>% filter(mutation %in% c("tardbp","ctrl"), treatment == "ml240") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","tardbp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|gli|cb1|ncrm"))
coip_ml240.tardbp_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.tardbp_vs_ctrl.nuc_cyt, experimental_design = coip_samples.ml240.tardbp_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "tardbp",
                                                     design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.nuc_cyt.tardbp_vs_ctrl")

coip_samples.ml240.vcp_r155c_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|cb1", cellline), treatment == "ml240") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|gli"))
coip_ml240.vcp_r155c_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_r155c_vs_ctrl.nuc_cyt, experimental_design = coip_samples.ml240.vcp_r155c_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                        design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.nuc_cyt.tardbp_vs_ctrl")

coip_samples.ml240.vcp_r191q_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|gli", cellline), treatment == "ml240") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|ncrm|cb1"))
coip_ml240.vcp_r191q_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_r191q_vs_ctrl.nuc_cyt, experimental_design = coip_samples.ml240.vcp_r191q_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                        design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.nuc_cyt.vcp_r191q_vs_ctrl")

coip_samples.ml240.vcp_iso_vs_ctrl.nuc_cyt = coip_samples %>% filter(grepl("ctrl|ncrm", cellline), treatment == "ml240") %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(mutation, levels = c("ctrl","vcp")), line_in_condition = recode_within(cellline, mutation)) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation, line_in_condition)
coip_lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_cyt = coip_lfq_intensities.unique %>% select(!matches("untreated|tdpn|gli|cb1"))
coip_ml240.vcp_iso_vs_ctrl.nuc_cyt.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ml240.vcp_iso_vs_ctrl.nuc_cyt, experimental_design = coip_samples.ml240.vcp_iso_vs_ctrl.nuc_cyt, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "ctrl", numerator = "vcp",
                                                      design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip_ml240.nuc_cyt.vcp_iso_vs_ctrl")

save(coip_ml240.tardbp_vs_ctrl.nuc_cyt.dep, coip_ml240.vcp_r155c_vs_ctrl.nuc_cyt.dep, coip_ml240.vcp_r191q_vs_ctrl.nuc_cyt.dep, coip_ml240.vcp_iso_vs_ctrl.nuc_cyt.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec.ml240.als_vs_ctrl.nuc_cyt.RData"))


# Co-IP ml240_vs_untreated --------

# Combine nuclear & cytoplasm samples

# CTRL
coip_samples.ctrl.nuc_cyt.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("ctrl")) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.ctrl.nuc_cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("tdpn|gli|cb1|ncrm"))
coip.ctrl.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ctrl.nuc_cyt.ml240_vs_untreated, experimental_design = coip_samples.ctrl.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                        design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.ctrl.nuc_cyt.ml240_vs_untreated", plot_dep = FALSE)
# TARDBP
coip_samples.tardbp.nuc_cyt.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("tardbp")) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.tardbp.nuc_cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|gli|cb1|ncrm"))
coip.tardbp.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.tardbp.nuc_cyt.ml240_vs_untreated, experimental_design = coip_samples.tardbp.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                        design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.tardbp.nuc_cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP R155C
coip_samples.vcp_r155c.nuc_cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("cb1", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.vcp_r155c.nuc_cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|ncrm"))
coip.vcp_r155c.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r155c.nuc_cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_r155c.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                          design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r155c.nuc_cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP R191Q
coip_samples.vcp_r191q.nuc_cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("gli", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.vcp_r191q.nuc_cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|cb1|ncrm"))
coip.vcp_r191q.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r191q.nuc_cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_r191q.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                             design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r191q.nuc_cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP knock-in
coip_samples.vcp_iso.nuc_cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("ncrm", cellline)) %>% group_by(cellline, fraction, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, fraction, mutation)
coip_lfq_intensities.unique.vcp_iso.nuc_cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|cb1"))
coip.vcp_iso.nuc_cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_iso.nuc_cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_iso.nuc_cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                             design_formula = ~0 + condition + fraction + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_iso.nuc_cyt.ml240_vs_untreated", plot_dep = FALSE)

save(coip.ctrl.nuc_cyt.ml240_vs_untreated.dep, coip.tardbp.nuc_cyt.ml240_vs_untreated.dep, coip.vcp_r155c.nuc_cyt.ml240_vs_untreated.dep, coip.vcp_r191q.nuc_cyt.ml240_vs_untreated.dep, coip.vcp_iso.nuc_cyt.ml240_vs_untreated.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec.nuc_cyt.ml240_vs_untreated.RData"))

# Nuclear samples separately
# CTRL
coip_samples.ctrl.nuc.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("ctrl"), fraction == "nuclear") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.ctrl.nuc.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("tdpn|gli|cb1|ncrm|cytoplasm"))
coip.ctrl.nuc.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ctrl.nuc.ml240_vs_untreated, experimental_design = coip_samples.ctrl.nuc.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                        design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.ctrl.nuc.ml240_vs_untreated", plot_dep = FALSE)
# TARDBP
coip_samples.tardbp.nuc.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("tardbp"), fraction == "nuclear") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.tardbp.nuc.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|gli|cb1|ncrm|cytoplasm"))
coip.tardbp.nuc.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.tardbp.nuc.ml240_vs_untreated, experimental_design = coip_samples.tardbp.nuc.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                          design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.tardbp.nuc.ml240_vs_untreated", plot_dep = FALSE)
# VCP R155C
coip_samples.vcp_r155c.nuc.ml240_vs_untreated = coip_samples %>% filter(grepl("cb1", cellline), fraction == "nuclear") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_r155c.nuc.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|ncrm|cytoplasm"))
coip.vcp_r155c.nuc.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r155c.nuc.ml240_vs_untreated, experimental_design = coip_samples.vcp_r155c.nuc.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                             design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r155c.nuc.ml240_vs_untreated", plot_dep = FALSE)
# VCP R191Q
coip_samples.vcp_r191q.nuc.ml240_vs_untreated = coip_samples %>% filter(grepl("gli", cellline), fraction == "nuclear") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_r191q.nuc.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|cb1|ncrm|cytoplasm"))
coip.vcp_r191q.nuc.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r191q.nuc.ml240_vs_untreated, experimental_design = coip_samples.vcp_r191q.nuc.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                             design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r191q.nuc.ml240_vs_untreated", plot_dep = FALSE)
# VCP knock-in
coip_samples.vcp_iso.nuc.ml240_vs_untreated = coip_samples %>% filter(grepl("ncrm", cellline), fraction == "nuclear") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_iso.nuc.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|cb1|cytoplasm"))
coip.vcp_iso.nuc.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_iso.nuc.ml240_vs_untreated, experimental_design = coip_samples.vcp_iso.nuc.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                           design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_iso.nuc.ml240_vs_untreated", plot_dep = FALSE)

# Cytoplasm samples separately
# CTRL
coip_samples.ctrl.cyt.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("ctrl"), fraction == "cytoplasm") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.ctrl.cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("tdpn|gli|cb1|ncrm|nuclear"))
coip.ctrl.cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.ctrl.cyt.ml240_vs_untreated, experimental_design = coip_samples.ctrl.cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                    design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.ctrl.cyt.ml240_vs_untreated", plot_dep = FALSE)
# TARDBP
coip_samples.tardbp.cyt.ml240_vs_untreated = coip_samples %>% filter(mutation %in% c("tardbp"), fraction == "cytoplasm") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.tardbp.cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|gli|cb1|ncrm|nuclear"))
coip.tardbp.cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.tardbp.cyt.ml240_vs_untreated, experimental_design = coip_samples.tardbp.cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                      design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.tardbp.cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP R155C
coip_samples.vcp_r155c.cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("cb1", cellline), fraction == "cytoplasm") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_r155c.cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|ncrm|nuclear"))
coip.vcp_r155c.cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r155c.cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_r155c.cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                         design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r155c.cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP R191Q
coip_samples.vcp_r191q.cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("gli", cellline), fraction == "cytoplasm") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_r191q.cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|cb1|ncrm|nuclear"))
coip.vcp_r191q.cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_r191q.cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_r191q.cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 5, denomenator = "untreated", numerator = "ml240",
                                                         design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_r191q.cyt.ml240_vs_untreated", plot_dep = FALSE)
# VCP knock-in
coip_samples.vcp_iso.cyt.ml240_vs_untreated = coip_samples %>% filter(grepl("ncrm", cellline), fraction == "cytoplasm") %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup %>%
  mutate(condition = factor(treatment, levels = c("untreated","ml240"))) %>%
  select(label = lfq_name, condition, replicate, technical_replicate, cellline, treatment, mutation)
coip_lfq_intensities.unique.vcp_iso.cyt.ml240_vs_untreated = coip_lfq_intensities.unique %>% select(!matches("ctrl|tdpn|gli|cb1|nuclear"))
coip.vcp_iso.cyt.ml240_vs_untreated.dep = DEP.analysis(protein_groups = coip_lfq_intensities.unique.vcp_iso.cyt.ml240_vs_untreated, experimental_design = coip_samples.vcp_iso.cyt.ml240_vs_untreated, show_proteins_per_threshold = TRUE, threshold = 4, denomenator = "untreated", numerator = "ml240",
                                                       design_formula = ~0 + condition + technical_replicate + cellline, file_path = here(proj_path, "mass-spec"), file_prefix = "coip.vcp_iso.cyt.ml240_vs_untreated", plot_dep = FALSE)

save(coip.ctrl.nuc.ml240_vs_untreated.dep, coip.tardbp.nuc.ml240_vs_untreated.dep, coip.vcp_r155c.nuc.ml240_vs_untreated.dep, coip.vcp_r191q.nuc.ml240_vs_untreated.dep, coip.vcp_iso.nuc.ml240_vs_untreated.dep,
     coip.ctrl.cyt.ml240_vs_untreated.dep, coip.tardbp.cyt.ml240_vs_untreated.dep, coip.vcp_r155c.cyt.ml240_vs_untreated.dep, coip.vcp_r191q.cyt.ml240_vs_untreated.dep, coip.vcp_iso.cyt.ml240_vs_untreated.dep,
     file = here(proj_path, "scripts/dep_coip_mass_spec.nuc_cyt_separate.ml240_vs_untreated.RData"))

