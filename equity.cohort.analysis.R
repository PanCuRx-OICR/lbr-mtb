#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(tidyr)
library(cowplot)
library(ggalluvial)
library(rlist)
library(PNWColors)
library(ggpubr)
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(DescTools)
library(survival)
library(ggfortify)
library(survminer)
library(Cairo)

'%ni%' <- function(x,y)!('%in%'(x,y))

fontSize = 28

#### 

basedir <- "~/Documents/scripts/bitbucket/"
setwd(basedir)
source(paste0(basedir, "/BTC.functions.R"))


BinomCI(x = 28, n = 55, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

#### FILE INTAKE ####

cohort <- fread("treatment-equity/source-data/cohort.txt")

cohort <- cohort[cohort$BTC_ID != "",]
cohort$BTCID <- gsub("_","",cohort$BTC_ID)

all_comer_btcids <- cohort$BTC_ID[cohort$COHORT %in% c('All-Comers')]

summary <- fread("temp_data/BTC.summary.csv")
summary$mmr_bin <- "MSS"
summary$mmr_bin[summary$mmr_score >= 4] <- "MSI"

all_somatic_variants <- fread("temp_data/BTC.somatic.variants.csv")

MTB_variants <- read_json("legresley/source-data/MTB_summary.json", simplifyVector = FALSE)
MTB_by_sample_raw <- pull_actions_from_MTB_json(MTB_variants)

oncokb_genes <- fread("temp_data/oncokb_biomarker_drug_associations.tsv")
INCLUDES_BILIARY = c("All Solid Tumors", "All Solid Tumors (excluding Colorectal Cancer)","All Solid Tumors (excluding Bladder Cancer)","All Solid Tumors (excluding Pancreatic Adenocarcinoma, Non-Small Cell Lung Cancer)","All Solid Tumors (excluding Thyroid Cancer, Non-Small Cell Lung Cancer)", "Biliary Tract Cancer, NOS", "Cholangiocarcinoma", "Intrahepatic Cholangiocarcinoma, Cholangiocarcinoma", "Hepatobiliary Cancer, Tubular Adenoma of the Colon, Esophagogastric Cancer, Anal Cancer, Small Bowel Cancer, Gastrointestinal Neuroendocrine Tumor")

THESE_LEVELS = c("1","2","3")
TRUNCATING_MUTATIONS <- c("deletion breakpoint","stopgain")

oncogenic_variants <- fread("temp_data/oncokb_oncogenic.list.txt",header=F)
names(oncogenic_variants) <- c("variant","gene")

fusions_raw <- fread('temp_data/BTC.star.fusions.txt', header=F)
names(fusions_raw)[1] <- 'donor'
fusions <- separate(fusions_raw, V3, c("gene","gene2"), sep = "--")

other_biomarkers <- fread('treatment-equity/source-data/other_biomarkers.txt', header=T)

#### profile ####

full_cohort_n <- nrow(cohort)
cat(full_cohort_n)
passing_samples_n <- nrow(cohort %>% filter(BIOPSY == "QC Pass"))
cat(passing_samples_n)

time_to_clinical_report <- as.numeric( as.Date(cohort$Clinical_Report_Issued) - as.Date(cohort$Specimen_Collection_Date) )
#as.numeric( as.Date("2023-02-16") - as.Date("2023-03-11") )
#as.numeric( as.Date("2023-08-09") - as.Date("2023-06-01") )

median(time_to_clinical_report, na.rm = T)
median(time_to_clinical_report, na.rm = T) - IQR(time_to_clinical_report, na.rm = T)
median(time_to_clinical_report, na.rm = T) + IQR(time_to_clinical_report, na.rm = T)

#### ONCOKB analysis ####

oncokb_genes <- oncokb_genes[oncokb_genes$Level %in% THESE_LEVELS,]

somatic_variants_oncoKB_filtered <- all_somatic_variants %>% filter(gene %in% oncokb_genes$Gene & 
                                donor %in% cohort$BTCID 
                              )

## add fusions from RNA
fusion_genes <- unique(oncokb_genes$Gene[oncokb_genes$Alterations == "Fusions"])
actionable_fusions <- unique(fusions[fusions$gene %in% fusion_genes ,-c(2,4)]) %>% filter(donor %in% cohort$BTCID)

actionable_variants_df <- get_oncokb_actionable_variants( somatic_variants_oncoKB_filtered,  
                                                          oncokb_genes, oncogenic_variants, 
                                                          TRUNCATING_MUTATIONS, INCLUDES_BILIARY, 
                                                          cohort=cohort, summary=summary, 
                                                          actionable_fusions=actionable_fusions, other_biomarkers=other_biomarkers)

#BTC0023 exonic TMB is 6
actionable_variants_df<- actionable_variants_df[actionable_variants_df$donor != "BTC0023" ,]

write.table(
  actionable_variants_df,
  file = "treatment-equity/results/actionable_variants.txt",
  row.names = FALSE, quote = FALSE, sep = "\t"
)

tableS1_oncokb <- actionable_variants_df %>% select(donor, gene, alteration_type, treatment) %>% unique()

tableS1_oncokb <- tableS1_oncokb %>%
  group_by(donor ,  gene  ,    alteration_type ) %>%
  summarise(treatmentd = paste(treatment, collapse = "; "), .groups = "drop")

tableS1_oncokb$alteration_type[tableS1_oncokb$gene == 'FGFR2'] <- 'fusion'
tableS1_oncokb$alteration_type[tableS1_oncokb$gene == 'TMB'] <- 'TMB-H'
tableS1_oncokb$alteration_type[tableS1_oncokb$gene == 'MSI'] <- 'MSI-H'

tableS1_oncokb$alteration_type[tableS1_oncokb$alteration_type == 'strong amplification'] <- 'amplification'
tableS1_oncokb$alteration_type[tableS1_oncokb$alteration_type == 'homozygous deletion'] <- 'deletion'

write.table(
  tableS1_oncokb,
  file = "treatment-equity/results/tableS1_oncokb.txt",
  row.names = FALSE, quote = FALSE, sep = "\t"
)


# seperately annotate other biomarkers
actionable_variants_df$alteration_class[actionable_variants_df$alteration_type == 'FISH'] <- "FISH"

actionable_variants <- unique(actionable_variants_df[,c("donor","gene","alteration_class","oncokb_level")])

# don't track specific mutations, lump as 'hotspot'
unique(actionable_variants$alteration_class)
actionable_variants$alteration_class[actionable_variants$alteration_class %in% c("V600", "V600E", "R132", "R132C", "R132H", "E545K", "H1047R", "G12D", "G12R", "E542K" )] <- "Hotspot"
actionable_variants$alteration_class[actionable_variants$alteration_class == 'Truncating Mutations'] <- "Truncating"


# only keep the best OncoKb hit for any mutation
actionable_variants <- actionable_variants %>% group_by(donor, gene, alteration_class) %>% summarize(top_oncokb_level = min(oncokb_level))
actionable_variants$top_oncokb_level[actionable_variants$top_oncokb_level == 3.5] <- '3B'

# incl. donors without actions in the onco plot
no_variants <- cohort[cohort$BTCID %ni% unique(actionable_variants$donor) & cohort$BIOPSY == "QC Pass",]
names(no_variants)[names(no_variants) == "BTCID"] <- "donor"
no_variants$gene <-NA
no_variants$alteration_class <- NA
no_variants$top_oncokb_level <- NA
no_variants <- no_variants[,c("donor","gene","alteration_class","top_oncokb_level")]

all_variants <- rbind.data.frame(actionable_variants, no_variants)

all_variants <- make_oncoplot_factor_order(all_variants)


actionable_frequency <- get_gene_frequencies_oncokb(all_variants)




outcomes <- cohort %>% filter(BIOPSY == "QC Pass") %>% select(BTCID, BTC_ID, OUTCOME, REASON, Funding, ACTION_AT_MTB)
outcomes$REASON[outcomes$OUTCOME %in% c("No Recurrence/Progression", "Deceased Without 2L", "Best Supportive Care")] <- NA
outcomes$REASON[outcomes$REASON == ''] <- NA
outcomes$Funding[outcomes$Funding == ''] <- NA

outcomes$donor_factor    <- factor(outcomes$BTCID,  levels = levels(unique(all_variants$donor_factor)))

this_theme <- theme(        
                  panel.background = element_rect(fill='transparent'), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color='transparent'), #transparent plot bg
                  panel.grid.major.x = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(color = NA,fill=NA), #transparent legend bg
                  legend.box.background = element_rect(color = NA,fill=NA), #transparent legend panel
                  
                  plot.title = element_blank(),
                  legend.spacing = unit(0, 'cm'),)

onco_plot <- 
  ggplot(all_variants , aes(x=donor_factor, y=gene_factor, fill=alteration_class)) + 
    geom_tile( na.rm = TRUE, color="white") +
    
    labs(x="Sample",y="",fill="") + 
    theme_bw(base_size = fontSize) + 
    this_theme + 
    theme(
      legend.position = "top" ,
          plot.margin = unit(c(0,1,0,1), "lines"),
          axis.text.y = element_text(face = "italic"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
         
         
    )+
  scale_fill_discrete( na.value = "white", na.translate = F) +
  scale_y_discrete(limits=rev, breaks=na.omit(unique(all_variants$gene_factor)))



fz_plot <- 
  ggplot(actionable_frequency , aes(x=frequency*100, y=gene_factor, fill=top_oncokb_level)) + 
    geom_bar(stat = "identity")+
    scale_y_discrete(limits=rev) +
    theme_bw(base_size = fontSize)+ 
    labs(x="Cohort Frequency (%)",fill="OncoKb Level") + 
    guides(fill = guide_legend(title.position="top")) +
    #scale_fill_grey()+
    this_theme + 
    scale_fill_manual(labels=c('1', '2', '3B', 'None'), values=c("indianred","#FF9900","#FFFF66","lightgrey")) + 
    theme(
      legend.position = "top" ,
          plot.margin = unit(c(0,1,0,-1), "lines"),
          axis.text.y =element_blank(),
          axis.title.y =element_blank(),
          axis.ticks.y =element_blank(),
    )

## comment if you'd like to grey out sample's w/out L2
treatment <- outcomes %>% select(donor_factor, REASON, OUTCOME)
treatment$REASON[treatment$REASON %in% c('FOLFOX', "Capecitabine", "GemCis Rechallenge")] <- NA
treatment$REASON[treatment$OUTCOME %in% c("No Recurrence/Progression", "Deceased Without 2L", "Best Supportive Care")] <- 'No 2L Therapy'
treatment <- treatment %>% select(donor_factor, REASON)
treatment$REASON2 <- treatment$REASON
treatment <- melt(treatment, id.vars= c("donor_factor"))
treatment$value[treatment$donor_factor == 'BTC0038' & treatment$variable == 'REASON2'] <- 'TDX1'

outcomes_plot <- 
  ggplot(treatment , aes(x=donor_factor, y=variable, fill=value)) + 
    geom_tile( na.rm = TRUE, color="white") +
    geom_text(label=ifelse(treatment$value == 'No 2L Therapy' & treatment$variable == 'REASON2',"X",""), alpha=0.5, size=8, vjust = 1.1) +
    labs(x="Sample",y='Targeted\nTherapy',fill="") + 
    theme_bw(base_size = fontSize) + 
    this_theme + 
  scale_fill_manual(values = c("forestgreen", "#64298e", "#94da40", "white","#1b511d", "#73c3e6", "#104b6d", "#f996f1", "orange"), na.value = 'white', na.translate = F) +  
    theme(
      plot.margin = unit(c(0,1,0,1), "lines"),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      legend.position = "bottom" ,
      panel.grid.major.y = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1),
      
    )
outcomes_legend =  get_plot_component(outcomes_plot, pattern = "guide-box-bottom" )
outcomes_plot = outcomes_plot + guides(fill="none")

funding <- outcomes %>% select(donor_factor, Funding, OUTCOME)
funding$Funding[funding$OUTCOME %in% c("No Recurrence/Progression", "Deceased Without 2L", "Best Supportive Care")] <- 'No 2L Therapy'
funding <- funding %>% select(donor_factor, Funding)
funding$funding2 <- funding$Funding
funding <- melt(funding, id.vars= c("donor_factor"))
funding$value[funding$donor_factor == 'BTC0038' & funding$variable == 'funding2'] <- 'Compassionate'

Funding_plot <- 
  ggplot(funding , aes(x=donor_factor, y=variable, fill=value)) + 
  geom_tile( na.rm = TRUE, color="white") +
  geom_text(label=ifelse(funding$value == 'No 2L Therapy' & funding$variable == 'funding2',"X",""), alpha=0.5, size=8, vjust = 1.1) +
  labs(x="",y="Funding",fill="") + 
  theme_bw(base_size = fontSize) + 
  this_theme + 
  scale_fill_manual(values = c("#00496f", "#0f85a0", "white", "#ed8b00", "#dd4124"), na.value = 'white', na.translate = F) +  
  theme(
    plot.margin = unit(c(0,1,0,1), "lines"),
    panel.grid.major.y = element_blank(),
    
    legend.position = "bottom" ,
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)

    axis.text.x =element_blank(),
    axis.title.x =element_blank(),
    axis.ticks.x =element_blank()
  )
funding_legend =  get_plot_component(Funding_plot, pattern = "guide-box-bottom" )
Funding_plot = Funding_plot + guides(fill="none")



width_gap = -0.11
gap = -0.24


svg(filename = paste0("treatment-equity/results/BTC.oncoplot.svg"), width = 18, height = 11)
print(
plot_grid(onco_plot,NULL, fz_plot,
          NULL, NULL, NULL,
          outcomes_plot,NULL, NULL, 
          NULL, NULL, NULL,
          Funding_plot,NULL, NULL, 
          NULL, NULL, NULL,
          outcomes_legend, NULL, NULL,
          funding_legend, NULL, NULL,
          rel_heights = c(0.8, gap,0.3, gap,0.3,0, 0.1, 0.05 ), rel_widths = c(1,width_gap,0.4 ), ncol = 3,align = 'hv', axis = "rlbt")
)
dev.off()

length(unique(actionable_variants_df$gene))
hit_samples_n <- length(unique(actionable_variants_df$donor))
target_tally <- unique(actionable_variants_df[, c('donor', 'gene')] ) %>% group_by(donor) %>% tally()
multihit_samples_n <- nrow(target_tally %>% filter(n > 1))

cat(multihit_samples_n, " of ", hit_samples_n)
BinomCI(x = multihit_samples_n, n = hit_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

tier_frequency <- get_tier_frequencies_oncokb(variants_by_sample=all_variants)

tier_one_n <- unlist(tier_frequency %>% filter(top_oncokb_level == 1) %>% select(tally))
cat(tier_one_n, " of ", full_cohort_n)
BinomCI(x = tier_one_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

tier_one_or_two_n <- unlist((tier_frequency %>% filter(top_oncokb_level == 2) %>% select(tally)) + tier_one_n)
cat(tier_one_or_two_n, " of ", full_cohort_n)
BinomCI(x = tier_one_or_two_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

tier_one_to_three_n <- unlist((tier_frequency %>% filter(top_oncokb_level == "3B") %>% select(tally)) + tier_one_or_two_n)
cat(tier_one_to_three_n, " of ", full_cohort_n)
BinomCI(x = tier_one_to_three_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))


hits_by_gene <-all_variants %>% filter(!is.na(gene)) %>% group_by(gene) %>% tally()
hits_by_gene <- hits_by_gene[order(-hits_by_gene$n),]
hits_by_gene

BinomCI(x = 5, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
BinomCI(x = 4, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
BinomCI(x = 3, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))


#### oncokb mtb comparison ####

  MTB_tally <- MTB_by_sample_raw %>% group_by(BTC_ID) %>% tally()
  MTB_tally$mtb_action <- "MTB match"
  
  actions_by_donor <- actionable_variants %>% group_by(donor) %>% tally()
  actions_by_donor$one_action <- "OncoKB Matched"

  cohort_outcomes <- left_join(outcomes, actions_by_donor, by=c("BTCID"="donor"))
  cohort_outcomes <- left_join(cohort_outcomes, MTB_tally, by=c("BTC_ID"="BTC_ID"))
  
  cohort_outcomes$one_action[is.na(cohort_outcomes$one_action)] <- "No OncoKB Match"
  cohort_outcomes$mtb_action[is.na(cohort_outcomes$mtb_action)] <- "No MTB Match"
  
  cohort_outcomes_tally <- cohort_outcomes %>% group_by(one_action, mtb_action) %>% tally()
  
  
  png(paste0("treatment-equity/results/BTC.balloon.png") ,width = 400, height = 300) #, bg = "transparent"
  
  ggplot(cohort_outcomes_tally, aes(x = one_action, y = mtb_action)) +
    geom_point(aes(size = n, color=mtb_action)) +
    scale_size(range = c(10, 60)) +
    scale_color_manual(values = pnw_palette(name="Bay",n=11)[c(3,8)]) +  
    guides(color='none', size='none')+
    geom_text(aes(
      label = n,
      y = as.numeric(as.factor(mtb_action)) , 
      x = as.numeric(as.factor(one_action)))
    ) + theme_bw(base_size = 20) +
    
    
    theme(    plot.margin = unit(c(0,1,0,1), "lines"),
              panel.grid.major = element_blank(),
              axis.title =element_blank(),
              legend.position = "bottom" ,
              axis.text.x = element_text(angle = 10, vjust = 1, hjust=1)
    )
  
  dev.off()
  
  write.table(
    cohort_outcomes,
    file = "treatment-equity/results/cohort_outcomes.txt",
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  outcome_tally <- cohort_outcomes %>% group_by(ACTION_AT_MTB, mtb_action) %>% tally()
  recommendations_n <- sum(outcome_tally$n[outcome_tally$ACTION_AT_MTB == "Recommendation Made"])
  alive_at_mtb_n <- sum(outcome_tally$n[outcome_tally$ACTION_AT_MTB != "Deceased before Recommendation"])
  
  cat(recommendations_n, " of ", full_cohort_n)
  BinomCI(x = recommendations_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

  #sub for all-comers (post 2023)
  outcome_tally_sub <- cohort_outcomes %>% filter(BTC_ID %in% all_comer_btcids) %>% group_by(ACTION_AT_MTB, mtb_action) %>% tally()
  recommendations_sub_n <- sum(outcome_tally_sub$n[outcome_tally_sub$ACTION_AT_MTB == "Recommendation Made"])
  alive_at_mtb_sub_n <- sum(outcome_tally_sub$n[outcome_tally_sub$ACTION_AT_MTB != "Deceased before Recommendation"])
  
  cat(recommendations_sub_n, " of ", alive_at_mtb_sub_n)
  BinomCI(x = recommendations_sub_n, n = alive_at_mtb_sub_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  mtb_onco_agree_n <- sum(cohort_outcomes_tally$n[(cohort_outcomes_tally$one_action == "OncoKB Matched" & 
                          cohort_outcomes_tally$mtb_action == "MTB match") | 
                            (cohort_outcomes_tally$one_action == "No OncoKB Match" & 
                               cohort_outcomes_tally$mtb_action == "No MTB Match")])
  
  cat(mtb_onco_agree_n, " of ", full_cohort_n)
  BinomCI(x = mtb_onco_agree_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  mtb_rescue_n <- sum(cohort_outcomes_tally$n[(cohort_outcomes_tally$one_action == "No OncoKB Match" & 
                                                       cohort_outcomes_tally$mtb_action == "MTB match")])
  
  cat(mtb_rescue_n, " of ", full_cohort_n)
  BinomCI(x = mtb_rescue_n, n = full_cohort_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  #### MTB oncoplot ####
  
  MTB_new_markers <- MTB_by_sample_raw %>% filter(gene %ni% c(unique(oncokb_genes$Gene), "TMB", "MMRD"))
  MTB_new_markers <- MTB_new_markers[,c("BTC_ID",	"gene",	"variant_type",	"variant", "therapy")]
  MTB_new_markers$OncoKb <- "new biomarker"
  
  MTB_in_oncoKB <- MTB_by_sample_raw %>% filter(gene %in% unique(oncokb_genes$Gene))
  MTB_in_oncoKB$BTCID <- gsub("_","",MTB_in_oncoKB$BTC_ID)

  actionable_variants_full <- fread('treatment-equity/results/actionable_variants.txt')
  actionable_variants_full <- unique(actionable_variants_full[,c('donor','gene', 'alteration_type')])
  
  MTB_known_actions <- left_join(MTB_in_oncoKB, actionable_variants_full, by=c("BTCID"="donor","gene"="gene"))
  MTB_unknown_actions <- MTB_known_actions %>% filter((variant != alteration_type & variant_type != alteration_type & gene != alteration_type) | is.na(alteration_type))
  MTB_unknown_actions <- MTB_unknown_actions[,c("BTC_ID",	"gene",	"variant_type",	"variant", "therapy")]
  MTB_unknown_actions$OncoKb <- "unknown variant effect"
  
  MTB_variants <- rbind.data.frame(MTB_new_markers, MTB_unknown_actions)
  
  write.table(
    MTB_variants[,-4],
    file = "treatment-equity/results/tableS1_MTB.txt",
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  MTB_variants <- cohort %>% filter(BIOPSY == "QC Pass") %>% select(BTCID, BTC_ID ) %>% left_join(MTB_variants[,-5], by=c("BTC_ID"="BTC_ID"))
  
  names(MTB_variants)[1] <- c("donor")
  MTB_variants <- make_oncoplot_factor_order(MTB_variants)
  
  MTB_variants$btc_factor    <- factor(MTB_variants$donor,  levels = levels(unique(all_variants$donor_factor)))
  
  names(MTB_variants)[6] <- c("top_oncokb_level")
  MTB_frequency <- get_gene_frequencies_oncokb(MTB_variants)
  
  
  black_onco_plot <- 
    ggplot(all_variants , aes(x=donor_factor, y=gene_factor, fill=alteration_class)) + 
    geom_tile( na.rm = TRUE, color="white") +
    
    labs(x="Sample",y="",fill="") + 
    theme_bw(base_size = fontSize) + 
    this_theme + 
    theme(
      legend.position = "top" ,
      plot.margin = unit(c(0,1,0,1), "lines"),
      axis.text.y = element_text(face = "italic"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      
      
    )+
    guides(fill='none')+
    scale_fill_manual(values = c(rep("black", 7)), na.value = "white", na.translate = F) +
    scale_y_discrete(limits=rev, breaks=na.omit(unique(all_variants$gene_factor)))
  
  
 MTB_plot <- 
    ggplot(MTB_variants , aes(x=btc_factor, y=gene_factor, fill=variant_type)) + 
    geom_tile( na.rm = TRUE, color="white") +
    
    labs(x="Sample",y="",fill="") + 
    theme_bw(base_size = fontSize) + 
    this_theme + 
    theme(
      legend.position = "bottom" ,
      plot.margin = unit(c(0,1,0,1), "lines"),
      axis.text.y = element_text(face = "italic"),
      axis.title.x=element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      
      
    )+
    scale_fill_discrete( na.value = "white", na.translate = F) +
    scale_y_discrete(limits=rev, breaks=na.omit(unique(MTB_variants$gene_factor)))
  
  
  
  MTB_fz_plot <- 
    ggplot(MTB_frequency , aes(x=frequency*100, y=gene_factor, fill=top_oncokb_level)) + 
    geom_bar(stat = "identity")+
    scale_y_discrete(limits=rev) +
    theme_bw(base_size = fontSize)+ 
    labs(x="Cohort Frequency (%)", fill="") + 
    guides(fill = guide_legend(title.position="top")) +
    scale_fill_grey()+
    this_theme + 
    scale_fill_discrete( na.value = "white", na.translate = F) +
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    theme(
      legend.position = "top" ,
      plot.margin = unit(c(0,1,0,-1), "lines"),
      axis.text.y =element_blank(),
      axis.title.y =element_blank(),
      axis.ticks.y =element_blank(),
    )
  
  
  
  width_gap = -0.11
  gap = -0.22
  
  png(paste0("treatment-equity/results/BTC.MTB.png") ,width = 1600, height = 1200) #, bg = "transparent"
  
  plot_grid(black_onco_plot,NULL, NULL,
            NULL, NULL, NULL,
            MTB_plot,NULL, MTB_fz_plot, 
            
            rel_heights = c(0.5, gap,0.5 ), rel_widths = c(1,width_gap,0.4 ), ncol = 3,align = 'hv', axis = "rlbt")
  
  dev.off()
  
  
  #MTAP counted in figure
  BinomCI(x = 6, n = mtb_rescue_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  MTB_frequency[order(-MTB_frequency$tally),]
  BinomCI(x = 11, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  #PTEN
  BinomCI(x = 3, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  

  #### SANKEY ####
  
  alluvial <- cohort %>% select(BTC_ID, BIOPSY, ACTION_AT_MTB, OUTCOME, REASON)
  alluvial[alluvial == ""] <- NA
  alluvial[alluvial == "Deceased before Recommendation"] <- "Deceased Before Recommendation"
  
  biopsy_tally <- alluvial %>% group_by(BIOPSY) %>% tally()
  biopsy_tally$frequency <- round((biopsy_tally$n / sum(biopsy_tally$n)) * 100)
  biopsy_tally$biopsy_string <- paste0(biopsy_tally$BIOPSY, "\n", biopsy_tally$frequency, "% (n=", biopsy_tally$n, ")")
  this_tally <- biopsy_tally %>% select(BIOPSY, biopsy_string) 
  alluvial <- left_join(alluvial, this_tally, by=c("BIOPSY"="BIOPSY"))
  
  action_tally <- alluvial %>% group_by(ACTION_AT_MTB) %>% tally()
  action_tally <- action_tally[!is.na(action_tally$ACTION_AT_MTB),]
  action_tally$frequency <- round((action_tally$n / sum(action_tally$n)) * 100)
  action_tally$action_string <- paste0(action_tally$ACTION_AT_MTB, "\n", action_tally$frequency, "% (n=", action_tally$n, ")")
  this_tally <- action_tally %>% select(ACTION_AT_MTB, action_string) 
  alluvial <- left_join(alluvial, this_tally, by=c("ACTION_AT_MTB"="ACTION_AT_MTB"))
  
  outcome_tally <- alluvial %>% group_by(OUTCOME) %>% tally()
  outcome_tally <- outcome_tally[!is.na(outcome_tally$OUTCOME),]
  outcome_tally$frequency <- round((outcome_tally$n / sum(outcome_tally$n)) * 100)
  outcome_tally$outcome_string <- paste0(outcome_tally$OUTCOME, "\n", outcome_tally$frequency, "% (n=", outcome_tally$n, ")")
  this_tally <- outcome_tally %>% select(OUTCOME, outcome_string) 
  alluvial <- left_join(alluvial, this_tally, by=c("OUTCOME"="OUTCOME"))
  
  reason_tally <- alluvial %>% group_by(REASON) %>% tally()
  reason_tally <- reason_tally[!is.na(reason_tally$REASON),]
  reason_tally$frequency <- round((reason_tally$n / sum(reason_tally$n)) * 100)
  reason_tally$reason_string <- paste0(reason_tally$REASON, "\n", reason_tally$frequency, "% (n=", reason_tally$n, ")")
  this_tally <- reason_tally %>% select(REASON, reason_string) 
  alluvial <- left_join(alluvial, this_tally, by=c("REASON"="REASON"))
  

  sankey_data <- alluvial %>%
    make_long(biopsy_string,  action_string, outcome_string, reason_string)
  
  unique(c(sankey_data$node, sankey_data$next_node))

 level_order <- c(
   
   "QC Fail\n4% (n=2)",
   "QC Pass\n96% (n=53)",
   
   "Deceased Before Recommendation\n8% (n=4)",
   "No Recommendation\n19% (n=10)" ,
   "Recommendation Made\n74% (n=39)",
   
   #"No 2L Therapy\n57% (n=28)",
   "Deceased or Best Supportive Care Without 2L\n10% (n=5)",
   "No Recurrence/Progression\n49% (n=24)"  ,
   "2L Therapy\n10% (n=5)",
   "Non-targeted 2L Therapy\n14% (n=7)",
   "Targeted Therapy\n16% (n=8)",
   
   
   "FOLFOX\n45% (n=9)",
   "GemCis Rechallenge\n10% (n=2)",
   "Capecitabine\n5% (n=1)" ,
   
   
   "Nivolumab\n5% (n=1)",
   "Everolimus\n5% (n=1)",
   "Pemigatinib\n10% (n=2)",
   "Durvalumab\n5% (n=1)",
   "Trastuzumab\n5% (n=1)",
   "Selpercatinib\n5% (n=1)",
   "Study Drug\n5% (n=1)" 
  )
 
  sankey_data$node    <- factor(sankey_data$node, levels = level_order)
  sankey_data$next_node    <- factor(sankey_data$next_node, levels = level_order)
  
  png(paste0("treatment-equity/results/BTC.targeted.sankey.png"), width = 1400, height = 1000, bg = "transparent")
  
  ggplot(sankey_data, aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
    geom_sankey() +
    geom_sankey_label(size=6, fill="white")+
    theme_sankey(base_size = 25)+
    
    scale_fill_manual(values = pnw_palette(name="Bay", n=22 )[c(1,6,11,18,22,4,14,2:3,5,7:10,12:13,15:17,19:21)], na.value = 'white', na.translate = F) +  
    scale_x_discrete(labels = c("Biopsy", "Action at MTB", "Treatment Type", "Treatment")) +
    theme(legend.position = "none", 
          axis.title.x = element_blank())

  
  dev.off()
  

  recommendations_n <- sum(outcome_tally$n[outcome_tally$ACTION_AT_MTB %in% c("Recommendation Made", "Deceased Before Recommendation") ])
  
  outcome_tally <- alluvial %>% group_by(OUTCOME, ACTION_AT_MTB) %>% tally() #
  no_recurrence_n <- na.omit(outcome_tally$n[outcome_tally$ACTION_AT_MTB == "Recommendation Made" & outcome_tally$OUTCOME == "No Recurrence/Progression"  ])
  cat(no_recurrence_n, " of ", recommendations_n)
  BinomCI(x = no_recurrence_n, n = recommendations_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))

  ## patients with targetable alterations either received targeted therapy in the first-line setting or received second-line treatments
  yes_2L_n <- sum(outcome_tally$n[outcome_tally$ACTION_AT_MTB == "Recommendation Made" & outcome_tally$OUTCOME %in% c("Targeted Therapy","Non-targeted 2L Therapy")]  , na.rm = T)
  yes_targeted_n <- sum(outcome_tally$n[outcome_tally$ACTION_AT_MTB == "Recommendation Made" & outcome_tally$OUTCOME %in% c("Targeted Therapy") ], na.rm = T)
  cat(yes_targeted_n, " of ", yes_2L_n)
  BinomCI(x = yes_targeted_n, n = yes_2L_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  
  ## received targeted first-line treatment
  BinomCI(x = 2, n = 55, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  ## received targeted therapies in the second- or third-line
  BinomCI(x = 6, n = yes_2L_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  
  #funded for CCA by the public health care system
  #BinomCI(x = 0, n = 8, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  #funded compassionately through programs from the pharmaceutical companies
  BinomCI(x = 6, n = 9, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  
  ## had targetable alterations, progressed on first-line therapy, but did not receive targeted therapies
  
  no_targeted_n <- sum(outcome_tally$n[outcome_tally$OUTCOME %in% c("No Targeted Therapy") & outcome_tally$ACTION_AT_MTB == "Recommendation Made"], na.rm = T)
  
  BinomCI(x = no_targeted_n, n = 55, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  ## no recommendation
  no_reco_n <- sum(outcome_tally$n[ outcome_tally$ACTION_AT_MTB == "No Recommendation"], na.rm = T)
  no_targeted_no_reco_n <- sum(outcome_tally$n[outcome_tally$OUTCOME %in% c("No Targeted Therapy") & outcome_tally$ACTION_AT_MTB == "No Recommendation"], na.rm = T)
  cat(no_targeted_no_reco_n, " of ", no_reco_n)
  BinomCI(x = no_targeted_no_reco_n, n = no_reco_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
  
  
  
  #### SWIMMERS ####

  swimmers <- cohort 
  swimmers <- swimmers[swimmers$OUTCOME %in% c("Targeted Therapy" ,  "Non-targeted 2L Therapy"), ]
  
  #swimmers$OUTCOME[swimmers$OUTCOME %in% c("2L Therapy")] <- "Non-targeted 2L Therapy"
  
  swimmers$full_length <- as.numeric(as.Date(swimmers$Last_follow_up) - as.Date(swimmers$Initial_Dx)) / 365
  swimmers <- swimmers[order(swimmers$full_length),]
  swimmers$BTCID_factor    <- factor(swimmers$BTCID, levels = swimmers$BTCID)
  

  swimmers$L1[swimmers$L1 == ''] <- NA
  swimmers$L2[swimmers$L2 == ''] <- NA
  swimmers$L3[swimmers$L3 == ''] <- NA
  
  for(this_swimmer_index in c(1:nrow(swimmers))){
    cat(this_swimmer_index, " ", swimmers$BTC_ID[this_swimmer_index], "\n")

    if(!is.na(swimmers$L3[this_swimmer_index]) & is.na(swimmers$L3_end[this_swimmer_index])){
      swimmers$L3_end[this_swimmer_index] <- swimmers$Last_follow_up[this_swimmer_index]
    } else if(!is.na(swimmers$L2[this_swimmer_index]) & is.na(swimmers$L2_End[this_swimmer_index])){
      swimmers$L2_End[this_swimmer_index] <- swimmers$Last_follow_up[this_swimmer_index]
      
    }else if(!is.na(swimmers$L1[this_swimmer_index]) & is.na(swimmers$L1_end[this_swimmer_index])){
      swimmers$L1_end[this_swimmer_index] <- swimmers$Last_follow_up[this_swimmer_index]
      
    }
    
    
  }
  
  swimmers$L1_start_adj <- as.numeric(as.Date(swimmers$L1_start) - as.Date(swimmers$Initial_Dx)) / 365
  swimmers$L1_end_adj <- as.numeric(as.Date(swimmers$L1_end) - as.Date(swimmers$Initial_Dx))/ 365
  
  swimmers$L2_start_adj <- as.numeric(as.Date(swimmers$L2_Start) - as.Date(swimmers$Initial_Dx)) / 365
  swimmers$L2_end_adj <- as.numeric(as.Date(swimmers$L2_End) - as.Date(swimmers$Initial_Dx)) / 365
  
  swimmers$L3_start_adj <- as.numeric(as.Date(swimmers$L3_Start) - as.Date(swimmers$Initial_Dx))/ 365
  
  
  swimmers$Deceased[swimmers$Deceased == ''] <- 'Alive'

  swimmers$L2[swimmers$L2 == 'Best Supportive Care'] <- NA
  
  swimmers$L1[swimmers$L1 %in% c("FOLFOX","FOLFIRINOX")] <- "FOLFOX or FOLFIRINOX"
  swimmers$L2[swimmers$L2 %in% c("FOLFOX","FOLFIRINOX")] <- "FOLFOX or FOLFIRINOX"
  swimmers$L3[swimmers$L3 %in% c("FOLFOX","FOLFIRINOX")] <- "FOLFOX or FOLFIRINOX"
  
  swimmers$L3[swimmers$L3 %in% c("Durva","GemCis","GemCisDurva")] <- "GemCis and/or Durva"
  swimmers$L2[swimmers$L2 %in% c("Durva","GemCis","GemCisDurva")] <- "GemCis and/or Durva"
  swimmers$L1[swimmers$L1 %in% c("Durva","GemCis","GemCisDurva")] <- "GemCis and/or Durva"
  
  svg(paste0("treatment-equity/results/BTC.targeted.swimmers.svg"), width = 10, height = 5)
  
  ggplot(swimmers, aes(y=BTCID_factor, x=full_length) ) + 
    geom_bar(stat="identity", aes(fill=OUTCOME)) +
    
    geom_segment(aes(x=full_length, xend=full_length+(80/365), alpha=ifelse(Deceased == "Deceased", "Alive", "Deceased"), ), arrow = arrow(length = unit(0.2, "cm"))) + 
    geom_point(aes(shape=Deceased), size=3) + 
    
    geom_segment(aes(x=L1_start_adj, xend=L1_end_adj, y=BTCID_factor, yend=BTCID_factor, color=L1), linewidth=3, na.rm = TRUE)+
    geom_segment(aes(x=L2_start_adj, xend=L2_end_adj, y=BTCID_factor, yend=BTCID_factor, color=L2), linewidth=3, na.rm = TRUE)+
    geom_segment(aes(x=L3_start_adj, xend=full_length, y=BTCID_factor, yend=BTCID_factor, color=L3), linewidth=3, na.rm = TRUE)+

    theme_bw(base_size = 18) +
    guides(alpha="none",  shape= "none", fill='none')+
    guides(color=guide_legend(nrow=3,byrow=TRUE)) +
    
    scale_shape_manual( values=c(NA, 20)) +
    scale_fill_manual(values = c("grey", "black")) +  
    scale_alpha_discrete(range = c(0, 1))+
    scale_color_manual(values = c("#48bf8e", "#782857", "#a3c9fe", "#3a508a", "#69e04c", "#8301bd", "#20d8fd", "#f90da0", "#207a3f", "#df7acb", "#1c35b7", "#a3c541"), na.value = 'white', na.translate = F) +  
    
    labs(x="Time Since Diagnosis (Years)", y='', shape='Status', fill="", color='Treatment') +
    theme(legend.position = "top",
         # legend.background = element_rect(size=0.5,   linetype="solid",  colour ="black")
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          panel.border = element_blank()
          )

    dev.off()
  
    swimmers_n <- nrow(swimmers)
    cat(swimmers_n, " of ", passing_samples_n)
    BinomCI(x = swimmers_n, n = passing_samples_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
    
    reco_tar_swimmers_n <- length(swimmers$BTC_ID[swimmers$ACTION_AT_MTB == 'Recommendation Made' & swimmers$OUTCOME == 'Targeted Therapy' ])
    reco_tar_swimmers_alive_n <- length(swimmers$BTC_ID[swimmers$ACTION_AT_MTB == 'Recommendation Made' & swimmers$OUTCOME == 'Targeted Therapy' & swimmers$Deceased == 'Deceased'])
    cat(reco_tar_swimmers_alive_n, " of ", reco_tar_swimmers_n)
    BinomCI(x = reco_tar_swimmers_alive_n, n = reco_tar_swimmers_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
    
    no_reco_swimmers_n <- length(swimmers$BTC_ID[swimmers$ACTION_AT_MTB %in% c('Recommendation Made','No Recommendation') & swimmers$OUTCOME == 'No Targeted Therapy' ])
    no_reco_swimmers_alive_n <- length(swimmers$BTC_ID[swimmers$ACTION_AT_MTB %in% c('Recommendation Made','No Recommendation') & swimmers$OUTCOME == 'No Targeted Therapy' & swimmers$Deceased == 'Deceased'])
    cat(no_reco_swimmers_alive_n, " of ", no_reco_swimmers_n)
    BinomCI(x = no_reco_swimmers_alive_n, n = no_reco_swimmers_n, conf.level = 0.95, sides = c("two.sided"), method = c("wilson"))
    
    
    
    ##  
    swimmers$time <- as.numeric(as.Date(swimmers$Last_follow_up) - as.Date(swimmers$Initial_Dx)) 
    swimmers$trt[swimmers$OUTCOME == "Targeted Therapy"] <- 1
    swimmers$trt[swimmers$OUTCOME == "No Targeted Therapy"] <- 0
    
    swimmers$status[swimmers$Deceased == "Deceased"] <- 1
    swimmers$status[swimmers$Deceased == "Alive"] <- 0
    
    
    km <- with(swimmers, Surv(time, status))
    km_trt_fit <- survfit(Surv(time, status) ~ OUTCOME, data=swimmers)
    summary(km_trt_fit)
  
    surv_pvalue(km_trt_fit)
    
    autoplot(km_trt_fit) + theme_bw()
    
    
    swimmers$L2_months <-    (swimmers$L2_end_adj - swimmers$L2_start_adj) * 12
   min( swimmers$L2_months[swimmers$OUTCOME == "Non-targeted 2L Therapy"] , na.rm =T)
   max( swimmers$L2_months[swimmers$OUTCOME == "Non-targeted 2L Therapy"] , na.rm =T)
   
    
    