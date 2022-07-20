library(dplyr)
library(stringr)
library(readxl)
library(MSnbase)
library(dplyr)
library(readr)
library(DEP)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(VIM)
library(PEMM)
library(tidyverse)
library(missForest)
library(scales)
library(matrixStats)

#Defining experiment labeling
wt_ed <-
  tibble(
    label = c(
      "NoBiotinPD_3",
      "NoBiotinPD_4",
      "NoBiotinPD_1",
      "NoBiotinPD_2",
      "NoBiotinPD_5",
      "NoBiotinPD_6",
      "NoHpg_3",
      "NoHpg_4",
      "NoHpg_1",
      "NoHpg_2",
      "NoHpg_5",
      "NoHpg_6",
      "0Hr_3",
      "0Hr_4",
      "0Hr_1",
      "0Hr_2",
      "0Hr_5",
      "0Hr_6",
      "05Hr_3",
      "05Hr_4",
      "05Hr_1",
      "05Hr_2",
      "05Hr_5",
      "05Hr_6",
      "10Hr_3",
      "10Hr_4",
      "10Hr_1",
      "10Hr_2",
      "10Hr_5",
      "10Hr_6",
      "15Hr_3",
      "15Hr_4",
      "15Hr_1",
      "15Hr_2",
      "15Hr_5",
      "15Hr_6",
      "20Hr_3",
      "20Hr_4",
      "20Hr_1",
      "20Hr_2",
      "20Hr_5",
      "20Hr_6",
      "30Hr_3",
      "30Hr_4",
      "30Hr_1",
      "30Hr_2",
      "30Hr_5",
      "30Hr_6"
    )
  ) %>%
  separate(label, c("condition", "replicate"), remove = FALSE)

a2234d_ed <-
  tibble(
    label = c(
      "NoBiotinPD_1",
      "NoBiotinPD_2",
      "NoBiotinPD_3",
      "NoBiotinPD_4",
      "NoBiotinPD_5",
      "NoBiotinPD_6",
      "NoHpg_1",
      "NoHpg_2",
      "NoHpg_3",
      "NoHpg_4",
      "NoHpg_5",
      "NoHpg_6",
      "0Hr_1",
      "0Hr_2",
      "0Hr_3",
      "0Hr_4",
      "0Hr_5",
      "0Hr_6",
      "05Hr_1",
      "05Hr_2",
      "05Hr_3",
      "05Hr_4",
      "05Hr_5",
      "05Hr_6",
      "10Hr_1",
      "10Hr_2",
      "10Hr_3",
      "10Hr_4",
      "10Hr_5",
      "10Hr_6",
      "15Hr_1",
      "15Hr_2",
      "15Hr_3",
      "15Hr_4",
      "15Hr_5",
      "15Hr_6",
      "20Hr_1",
      "20Hr_2",
      "20Hr_3",
      "20Hr_4",
      "20Hr_5",
      "20Hr_6",
      "30Hr_1",
      "30Hr_2",
      "30Hr_3",
      "30Hr_4",
      "30Hr_5",
      "30Hr_6"
    )
  ) %>%
  separate(label, c("condition", "replicate"), remove = FALSE)

c1264r_ed <-
  tibble(
    label = c(
      "NoBiotinPD_1",
      "NoBiotinPD_2",
      "NoBiotinPD_3",
      "NoBiotinPD_4",
      "NoBiotinPD_5",
      "NoBiotinPD_6",
      "NoHpg_1",
      "NoHpg_2",
      "NoHpg_3",
      "NoHpg_4",
      "NoHpg_5",
      "NoHpg_6",
      "0Hr_1",
      "0Hr_2",
      "0Hr_3",
      "0Hr_4",
      "0Hr_5",
      "0Hr_6",
      "05Hr_1",
      "05Hr_2",
      "05Hr_3",
      "05Hr_4",
      "05Hr_5",
      "05Hr_6",
      "10Hr_1",
      "10Hr_2",
      "10Hr_3",
      "10Hr_4",
      "10Hr_5",
      "10Hr_6",
      "15Hr_1",
      "15Hr_2",
      "15Hr_3",
      "15Hr_4",
      "15Hr_5",
      "15Hr_6",
      "20Hr_1",
      "20Hr_2",
      "20Hr_3",
      "20Hr_4",
      "20Hr_5",
      "20Hr_6",
      "30Hr_1",
      "30Hr_2",
      "30Hr_3",
      "30Hr_4",
      "30Hr_5",
      "30Hr_6"
    )
  ) %>%
  separate(label, c("condition", "replicate"), remove = FALSE)

ml240_ed <-
  tibble(
    label = c(
      "NoBiotinPD_1",
      "NoBiotinPD_2",
      "NoHpg_1",
      "NoHpg_2",
      "0Hr_1",
      "0Hr_2",
      "05Hr_1",
      "05Hr_2",
      "10Hr_1",
      "10Hr_2",
      "15Hr_1",
      "15Hr_2",
      "20Hr_1",
      "20Hr_2",
      "30Hr_1",
      "30Hr_2"
    )
  ) %>%
  separate(label, c("condition", "replicate"), remove = FALSE)



#import and tidy data

wt_input <- read_csv(file.choose())
a2234d_input <- read_csv(file.choose())
c1264r_input <- read_csv(file.choose())
ml240_input <- read_csv(file.choose())

wt_unique <- make_unique(wt_input, "Gene Symbol", "Accession", delim = ";")
a2234d_unique <- make_unique(a2234d_input, "Gene Symbol", "Accession", delim = ";")
c1264r_unique <- make_unique(c1264r_input, "Gene Symbol", "Accession", delim = ";")
ml240_unique <- make_unique(ml240_input, "Gene Symbol", "Accession", delim = ";")

wt_columns <- grep("Abundance:", colnames(wt_unique))
a2234d_columns <- grep("Abundance:", colnames(a2234d_unique))
c1264r_columns <- grep("Abundance:", colnames(c1264r_unique))
ml240_columns <- grep("Abundance:", colnames(ml240_unique))

colnames(wt_unique)[wt_columns] <-
  paste("TMT", wt_ed$label, sep = ".")
colnames(a2234d_unique)[a2234d_columns] <-
  paste("TMT", a2234d_ed$label, sep = ".")
colnames(c1264r_unique)[c1264r_columns] <-
  paste("TMT", c1264r_ed$label, sep = ".")
colnames(ml240_unique)[ml240_columns] <-
  paste("TMT", ml240_ed$label, sep = ".")

wt_unique <- wt_unique[rowSums(is.na(wt_unique[,57:104])) != ncol(wt_unique[,57:104]),]
a2234d_unique <- a2234d_unique[rowSums(is.na(a2234d_unique[,57:104])) != ncol(a2234d_unique[,57:104]),]
c1264r_unique <- c1264r_unique[rowSums(is.na(c1264r_unique[,64:115])) != ncol(c1264r_unique[,64:115]),]
ml240_unique <- ml240_unique[rowSums(is.na(ml240_unique[,58:73])) != ncol(ml240_unique[,58:73]),]

# Define summerized experiment
wt_se <- make_se(wt_unique, wt_columns, wt_ed)
a2234d_se <- make_se(a2234d_unique, a2234d_columns, a2234d_ed)
c1264r_se <- make_se(c1264r_unique, c1264r_columns, c1264r_ed)
ml240_se <- make_se(ml240_unique, ml240_columns, ml240_ed)

# Plot and filter out proteins that contain missing values
wt_filt <- filter_proteins(wt_se, "fraction", min = 0.5)
a2234d_filt <- filter_proteins(a2234d_se, "fraction", min = 0.5)
c1264r_filt <- filter_proteins(c1264r_se, "fraction", min = 0.5)

wt_filt_diff <- test_diff(wt_filt, type ="manual", test = "NoBiotinPD_vs_NoHpg")
a2234d_filt_diff <- test_diff(a2234d_filt, type ="manual", test = "NoBiotinPD_vs_NoHpg")
c1264r_filt_diff <- test_diff(c1264r_filt, type ="manual", test = "NoBiotinPD_vs_NoHpg")

wt_filt_dep <- add_rejections(wt_filt_diff, alpha = 0.05, lfc = 1)
a2234d_filt_dep <- add_rejections(a2234d_filt_diff, alpha = 0.05, lfc = 1)
c1264r_filt_dep <- add_rejections(c1264r_filt_diff, alpha = 0.05, lfc = 1)

wt_filt_table <- get_results(wt_filt_dep)
a2234d_filt_table <- get_results(a2234d_filt_dep)
c1264r_filt_table <- get_results(c1264r_filt_dep)

wt_filt_table <- wt_filt_table %>%
  mutate(NoBiotinPD_fdr = p.adjust(wt_filt_table$NoBiotinPD_vs_NoHpg_p.val, method ="fdr"))
a2234d_filt_table <- a2234d_filt_table %>%
  mutate(NoBiotinPD_fdr = p.adjust(a2234d_filt_table$NoBiotinPD_vs_NoHpg_p.val, method ="fdr"))
c1264r_filt_table <- c1264r_filt_table %>%
  mutate(NoBiotinPD_fdr = p.adjust(c1264r_filt_table$NoBiotinPD_vs_NoHpg_p.val, method ="fdr"))

#Complex volcano plot annotation

#Upload Table of MCP Sorted Interactors

MCP_interactors <- read_csv(file.choose()) %>%
  as_tibble()

#Filt Datasets
#Annotate proteins above cutoff
wt_filt_table <- wt_filt_table %>%  
  mutate(NoBPD_accent =
           if_else(NoBiotinPD_fdr <= 0.05 &
                     (NoBiotinPD_vs_NoHpg_ratio >= 2*sd(wt_filt_ratios$NoBiotinPD_vs_NoHpg_ratio))
                   , 
                   TRUE, FALSE))

#Annotate proteins from MCP paper
wt_filt_table <- wt_filt_table %>%  
  mutate(MCP_accent =
           if_else(name %in% MCP_interactors$Gene
                   , 
                   TRUE, FALSE))

#Annotate TG
wt_filt_table <- wt_filt_table %>%  
  mutate(TG_accent =
           if_else(name == "TG"
                   , 
                   TRUE, FALSE))

ggplot(data = wt_filt_table, mapping = aes(x = NoBiotinPD_vs_NoHpg_ratio, y = -log10(NoBiotinPD_fdr))) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(2*sd(wt_filt_ratios$NoBiotinPD_vs_NoHpg_ratio)), linetype = 2, color = "grey40") +   
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  scale_color_manual(values = c("grey60", "dodgerblue"), guide = F) + 
  geom_point(data = wt_filt_table %>%
               dplyr::filter(NoBPD_accent),
             color = "dodgerblue") +
  geom_point(data = wt_filt_table %>%
               dplyr::filter(MCP_accent),
             color = "forestgreen", size = 2.5) +
  geom_point(data = wt_filt_table %>%
               dplyr::filter(TG_accent),
             color = "red4", size = 3) +
  geom_text_repel(data = wt_filt_table %>%
                    dplyr::filter(MCP_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black", min.segment.length = 0) +
  geom_text_repel(data = wt_filt_table %>%
                    dplyr::filter(TG_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black") +
  annotate("text", x = -3, y = -log10(0.01), label = "(-) Biotin PD vs (-) Hpg", 
           vjust = -35, hjust = 0.12) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17)) +
  labs(x = expression("Log"[2]*"Fold Change"), 
       y = expression("-Log"[10]*"FDR")) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

a2234d_filt_table <- a2234d_filt_table %>%  
  mutate(NoBPD_accent =
           if_else(NoBiotinPD_fdr <= 0.05 &
                     (NoBiotinPD_vs_NoHpg_ratio >= 2*sd(a2234d_filt_ratios$NoBiotinPD_vs_NoHpg_ratio))
                   , 
                   TRUE, FALSE))

#Annotate proteins from MCP paper
a2234d_filt_table <- a2234d_filt_table %>%  
  mutate(MCP_accent =
           if_else(name %in% MCP_interactors$Gene
                   , 
                   TRUE, FALSE))

#Annotate TG
a2234d_filt_table <- a2234d_filt_table %>%  
  mutate(TG_accent =
           if_else(name == "TG"
                   , 
                   TRUE, FALSE))


ggplot(data = a2234d_filt_table, mapping = aes(x = NoBiotinPD_vs_NoHpg_ratio, y = -log10(NoBiotinPD_fdr))) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(2*sd(a2234d_filt_ratios$NoBiotinPD_vs_NoHpg_ratio)), linetype = 2, color = "grey40") +   
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  scale_color_manual(values = c("grey60", "dodgerblue"), guide = F) + 
  geom_point(data = a2234d_filt_table %>%
               dplyr::filter(NoBPD_accent),
             color = "dodgerblue") +
  geom_point(data = a2234d_filt_table %>%
               dplyr::filter(MCP_accent),
             color = "forestgreen", size = 2.5) +
  geom_point(data = a2234d_filt_table %>%
               dplyr::filter(TG_accent),
             color = "red4", size = 3) +
  geom_text_repel(data = a2234d_filt_table %>%
                    dplyr::filter(MCP_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black") +
  geom_text_repel(data = a2234d_filt_table %>%
                    dplyr::filter(TG_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black") +
  annotate("text", x = -3, y = -log10(0.01), label = "(-) Biotin PD vs (-) Hpg", 
           vjust = -35, hjust = 0.12) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17)) +
  labs(x = expression("Log"[2]*"Fold Change"), 
       y = expression("-Log"[10]*"FDR")) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

c1264r_filt_table <- c1264r_filt_table %>%  
  mutate(NoBPD_accent =
           if_else(NoBiotinPD_fdr <= 0.05 &
                     (NoBiotinPD_vs_NoHpg_ratio >= 2*sd(c1264r_filt_ratios$NoBiotinPD_vs_NoHpg_ratio))
                   , 
                   TRUE, FALSE))

#Annotate proteins from MCP paper
c1264r_filt_table <- c1264r_filt_table %>%  
  mutate(MCP_accent =
           if_else(name %in% MCP_interactors$Gene
                   , 
                   TRUE, FALSE))

#Annotate TG
c1264r_filt_table <- c1264r_filt_table %>%  
  mutate(TG_accent =
           if_else(name == "TG"
                   , 
                   TRUE, FALSE))

ggplot(data = c1264r_filt_table, mapping = aes(x = NoBiotinPD_vs_NoHpg_ratio, y = -log10(NoBiotinPD_fdr))) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(2*sd(c1264r_filt_ratios$NoBiotinPD_vs_NoHpg_ratio)), linetype = 2, color = "grey40") +   
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  scale_color_manual(values = c("grey60", "dodgerblue"), guide = F) + 
  geom_point(data = c1264r_filt_table %>%
               dplyr::filter(NoBPD_accent),
             color = "dodgerblue") + 
  geom_point(data = c1264r_filt_table %>%
               dplyr::filter(MCP_accent),
             color = "forestgreen", size = 2.5) +
  geom_point(data = c1264r_filt_table %>%
               dplyr::filter(TG_accent),
             color = "red4", size = 3) +
  geom_text_repel(data = c1264r_filt_table %>%
                    dplyr::filter(MCP_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black") +
  geom_text_repel(data = c1264r_filt_table %>%
                    dplyr::filter(TG_accent),
                  aes(label = name),
                  hjust = 1,
                  segment.color = NA,
                  color = "black") +
  annotate("text", x = -3, y = -log10(0.05), label = "(-) Biotin PD vs (-) Hpg", 
           vjust = -35, hjust = 0.12) +
  scale_x_continuous(breaks = c(-3, -2, -1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)) +
  labs(x = expression("Log"[2]*"Fold Change"), 
       y = expression("-Log"[10]*"FDR")) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

#Upload Table of Pathway Sorted Interactors

interactors <- read_csv(file.choose()) %>%
  as_tibble()

interactors_ribosome <- interactors$Ribosome %>%
  na.omit()
interactors_import <- interactors$Import %>%
  na.omit()
interactors_NGlyco <- interactors$Glycosylation %>%
  na.omit()
interactors_Hsp <- interactors$`Hsp40/70/90 & Other Chaperoning Activity` %>%
  na.omit()
interactors_redox <- interactors$`Disulfide Bond Formation` %>%
  na.omit()
interactors_collagen <- interactors$`Proline Isomerization, Hydroxylation, & Collagen Processing` %>%
  na.omit()
interactors_traffic <- interactors$Trafficking %>% 
  na.omit()
interactors_ERAD <- interactors$`Protein Ubiquitination & Proteasomal Degradation` %>%
  na.omit()
interactors_autophagy <- interactors$`Autophagy & Lysosomal Processing` %>% 
  na.omit()
interactors_misc <- interactors$Miscellaneous %>%
  na.omit()

ordered_interactors <- c("TG", interactors_ribosome, interactors_import, interactors_NGlyco, interactors_Hsp, interactors_redox, interactors_collagen, interactors_traffic,
                         interactors_ERAD, interactors_autophagy, interactors_misc)

#Pull TMT chase data and format For heatmaps

wt_unique_tmt <- select(wt_unique, starts_with("TMT.")) %>%
  mutate(name = wt_unique$name, ID = wt_unique$ID) %>%
  filter(grepl(paste(ordered_interactors, collapse = "|"), name)) %>%
  as_tibble() %>%
  relocate(name, ID)

a2234d_unique_tmt <- select(a2234d_unique, starts_with("TMT.")) %>%
  mutate(name = a2234d_unique$name, ID = a2234d_unique$ID) %>%
  filter(grepl(paste(ordered_interactors, collapse = "|"), name)) %>%
  as_tibble() %>%
  relocate(name, ID)

c1264r_unique_tmt <- select(c1264r_unique, starts_with("TMT.")) %>%
  mutate(name = c1264r_unique$name, ID = c1264r_unique$ID) %>%
  filter(grepl(paste(ordered_interactors, collapse = "|"), name)) %>%
  as_tibble() %>%
  relocate(name, ID)

ml240_unique_tmt <- select(ml240_unique, starts_with("TMT.")) %>%
  mutate(name = ml240_unique$name, ID = ml240_unique$ID) %>%
  filter(grepl(paste(ordered_interactors, collapse = "|"), name)) %>%
  as_tibble() %>%
  relocate(name, ID)

wt_unique_tmt[,3:50] <- log(wt_unique_tmt[,3:50], base = 2)
a2234d_unique_tmt[,3:50] <- log(a2234d_unique_tmt[,3:50], base = 2)
c1264r_unique_tmt[,3:50] <- log(c1264r_unique_tmt[,3:50], base = 2)
ml240_unique_tmt[,3:18] <- log(ml240_unique_tmt[,3:18], base = 2)

#replace negative values with zero
wt_unique_tmt[wt_unique_tmt < 0] <- 0
a2234d_unique_tmt[a2234d_unique_tmt < 0] <- 0
c1264r_unique_tmt[c1264r_unique_tmt < 0] <- 0
ml240_unique_tmt[ml240_unique_tmt < 0] <- 0

#replace NAs
wt_unique_tmt <- wt_unique_tmt %>%
  mutate(TMT.NoHpg_1 = replace_na(wt_unique_tmt$TMT.NoHpg_1, median(wt_unique_tmt$TMT.NoHpg_1, na.rm = TRUE)),
         TMT.NoHpg_2 = replace_na(wt_unique_tmt$TMT.NoHpg_2, median(wt_unique_tmt$TMT.NoHpg_2, na.rm = TRUE)),
         TMT.NoHpg_3 = replace_na(wt_unique_tmt$TMT.NoHpg_3, median(wt_unique_tmt$TMT.NoHpg_3, na.rm = TRUE)),
         TMT.NoHpg_4 = replace_na(wt_unique_tmt$TMT.NoHpg_4, median(wt_unique_tmt$TMT.NoHpg_4, na.rm = TRUE)),
         TMT.NoHpg_5 = replace_na(wt_unique_tmt$TMT.NoHpg_5, median(wt_unique_tmt$TMT.NoHpg_5, na.rm = TRUE)),
         TMT.NoHpg_6 = replace_na(wt_unique_tmt$TMT.NoHpg_6, median(wt_unique_tmt$TMT.NoHpg_6, na.rm = TRUE)))

a2234d_unique_tmt <- a2234d_unique_tmt %>%
  mutate(TMT.NoHpg_1 = replace_na(a2234d_unique_tmt$TMT.NoHpg_1, median(a2234d_unique_tmt$TMT.NoHpg_1, na.rm = TRUE)),
         TMT.NoHpg_2 = replace_na(a2234d_unique_tmt$TMT.NoHpg_2, median(a2234d_unique_tmt$TMT.NoHpg_2, na.rm = TRUE)),
         TMT.NoHpg_3 = replace_na(a2234d_unique_tmt$TMT.NoHpg_3, median(a2234d_unique_tmt$TMT.NoHpg_3, na.rm = TRUE)),
         TMT.NoHpg_4 = replace_na(a2234d_unique_tmt$TMT.NoHpg_4, median(a2234d_unique_tmt$TMT.NoHpg_4, na.rm = TRUE)),
         TMT.NoHpg_5 = replace_na(a2234d_unique_tmt$TMT.NoHpg_5, median(a2234d_unique_tmt$TMT.NoHpg_5, na.rm = TRUE)),
         TMT.NoHpg_6 = replace_na(a2234d_unique_tmt$TMT.NoHpg_6, median(a2234d_unique_tmt$TMT.NoHpg_6, na.rm = TRUE)))

c1264r_unique_tmt <- c1264r_unique_tmt %>%
  mutate(TMT.NoHpg_1 = replace_na(c1264r_unique_tmt$TMT.NoHpg_1, median(c1264r_unique_tmt$TMT.NoHpg_1, na.rm = TRUE)),
         TMT.NoHpg_2 = replace_na(c1264r_unique_tmt$TMT.NoHpg_2, median(c1264r_unique_tmt$TMT.NoHpg_2, na.rm = TRUE)),
         TMT.NoHpg_3 = replace_na(c1264r_unique_tmt$TMT.NoHpg_3, median(c1264r_unique_tmt$TMT.NoHpg_3, na.rm = TRUE)),
         TMT.NoHpg_4 = replace_na(c1264r_unique_tmt$TMT.NoHpg_4, median(c1264r_unique_tmt$TMT.NoHpg_4, na.rm = TRUE)),
         TMT.NoHpg_5 = replace_na(c1264r_unique_tmt$TMT.NoHpg_5, median(c1264r_unique_tmt$TMT.NoHpg_5, na.rm = TRUE)),
         TMT.NoHpg_6 = replace_na(c1264r_unique_tmt$TMT.NoHpg_6, median(c1264r_unique_tmt$TMT.NoHpg_6, na.rm = TRUE)))

ml240_unique_tmt <- ml240_unique_tmt %>%
  mutate(TMT.NoHpg_1 = replace_na(ml240_unique_tmt$TMT.NoHpg_1, median(ml240_unique_tmt$TMT.NoHpg_1, na.rm = TRUE)),
         TMT.NoHpg_2 = replace_na(ml240_unique_tmt$TMT.NoHpg_2, median(ml240_unique_tmt$TMT.NoHpg_2, na.rm = TRUE)))

#Pull values for TG Normalization

wt_tg_vec1 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_1")) %>%
  as.numeric()

wt_tg_vec2 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_2")) %>%
  as.numeric()

wt_tg_vec3 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_3")) %>%
  as.numeric()

wt_tg_vec4 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_4")) %>%
  as.numeric()

wt_tg_vec5 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_5")) %>%
  as.numeric()

wt_tg_vec6 <- wt_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_6")) %>%
  as.numeric()  


a2234d_tg_vec1 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_1")) %>%
  as.numeric()

a2234d_tg_vec2 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_2")) %>%
  as.numeric()

a2234d_tg_vec3 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_3")) %>%
  as.numeric()

a2234d_tg_vec4 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_4")) %>%
  as.numeric()

a2234d_tg_vec5 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_5")) %>%
  as.numeric()

a2234d_tg_vec6 <- a2234d_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_6")) %>%
  as.numeric()

c1264r_tg_vec1 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_1")) %>%
  as.numeric()

c1264r_tg_vec2 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_2")) %>%
  as.numeric()

c1264r_tg_vec3 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_3")) %>%
  as.numeric()

c1264r_tg_vec4 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_4")) %>%
  as.numeric()

c1264r_tg_vec5 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_5")) %>%
  as.numeric()

c1264r_tg_vec6 <- c1264r_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_6")) %>%
  as.numeric()

ml240_tg_vec1 <- ml240_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_1")) %>%
  as.numeric()

ml240_tg_vec2 <- ml240_unique_tmt %>%
  filter(name == "TG") %>%
  select(ends_with("Hr_2")) %>%
  as.numeric()

#Normalize values to Tg abundance replicate by replicate

wt_unique_tmt <- mutate(wt_unique_tmt, 
                        X0Hr_Tg_Norm_1 = wt_unique_tmt$TMT.0Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[1]), 
                        X05Hr_Tg_Norm_1 = wt_unique_tmt$TMT.05Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[2]), 
                        X10Hr_Tg_Norm_1 = wt_unique_tmt$TMT.10Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[3]), 
                        X15Hr_Tg_Norm_1 = wt_unique_tmt$TMT.15Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[4]), 
                        X20Hr_Tg_Norm_1 = wt_unique_tmt$TMT.20Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[5]), 
                        X30Hr_Tg_Norm_1 = wt_unique_tmt$TMT.30Hr_1 * (mean(wt_tg_vec1)/wt_tg_vec1[6]),
                        X0Hr_Tg_Norm_2 = wt_unique_tmt$TMT.0Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[1]), 
                        X05Hr_Tg_Norm_2 = wt_unique_tmt$TMT.05Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[2]), 
                        X10Hr_Tg_Norm_2 = wt_unique_tmt$TMT.10Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[3]), 
                        X15Hr_Tg_Norm_2 = wt_unique_tmt$TMT.15Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[4]), 
                        X20Hr_Tg_Norm_2 = wt_unique_tmt$TMT.20Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[5]),
                        X30Hr_Tg_Norm_2 = wt_unique_tmt$TMT.30Hr_2 * (mean(wt_tg_vec2)/wt_tg_vec2[6]),
                        X0Hr_Tg_Norm_3 = wt_unique_tmt$TMT.0Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[1]), 
                        X05Hr_Tg_Norm_3 = wt_unique_tmt$TMT.05Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[2]), 
                        X10Hr_Tg_Norm_3 = wt_unique_tmt$TMT.10Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[3]), 
                        X15Hr_Tg_Norm_3 = wt_unique_tmt$TMT.15Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[4]), 
                        X20Hr_Tg_Norm_3 = wt_unique_tmt$TMT.20Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[5]),
                        X30Hr_Tg_Norm_3 = wt_unique_tmt$TMT.30Hr_3 * (mean(wt_tg_vec3)/wt_tg_vec3[6]),
                        X0Hr_Tg_Norm_4 = wt_unique_tmt$TMT.0Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[1]), 
                        X05Hr_Tg_Norm_4 = wt_unique_tmt$TMT.05Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[2]), 
                        X10Hr_Tg_Norm_4 = wt_unique_tmt$TMT.10Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[3]), 
                        X15Hr_Tg_Norm_4 = wt_unique_tmt$TMT.15Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[4]), 
                        X20Hr_Tg_Norm_4 = wt_unique_tmt$TMT.20Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[5]),
                        X30Hr_Tg_Norm_4 = wt_unique_tmt$TMT.30Hr_4 * (mean(wt_tg_vec4)/wt_tg_vec4[6]),
                        X0Hr_Tg_Norm_5 = wt_unique_tmt$TMT.0Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[1]), 
                        X05Hr_Tg_Norm_5 = wt_unique_tmt$TMT.05Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[2]), 
                        X10Hr_Tg_Norm_5 = wt_unique_tmt$TMT.10Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[3]), 
                        X15Hr_Tg_Norm_5 = wt_unique_tmt$TMT.15Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[4]), 
                        X20Hr_Tg_Norm_5 = wt_unique_tmt$TMT.20Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[5]),
                        X30Hr_Tg_Norm_5 = wt_unique_tmt$TMT.30Hr_5 * (mean(wt_tg_vec5)/wt_tg_vec5[6]),
                        X0Hr_Tg_Norm_6 = wt_unique_tmt$TMT.0Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[1]), 
                        X05Hr_Tg_Norm_6 = wt_unique_tmt$TMT.05Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[2]), 
                        X10Hr_Tg_Norm_6 = wt_unique_tmt$TMT.10Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[3]), 
                        X15Hr_Tg_Norm_6 = wt_unique_tmt$TMT.15Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[4]), 
                        X20Hr_Tg_Norm_6 = wt_unique_tmt$TMT.20Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[5]),
                        X30Hr_Tg_Norm_6 = wt_unique_tmt$TMT.30Hr_6 * (mean(wt_tg_vec6)/wt_tg_vec6[6]))

a2234d_unique_tmt <- mutate(a2234d_unique_tmt, 
                            X0Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.0Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[1]), 
                            X05Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.05Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[2]), 
                            X10Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.10Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[3]), 
                            X15Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.15Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[4]), 
                            X20Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.20Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[5]), 
                            X30Hr_Tg_Norm_1 = a2234d_unique_tmt$TMT.30Hr_1 * (mean(a2234d_tg_vec1)/a2234d_tg_vec1[6]),
                            X0Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.0Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[1]), 
                            X05Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.05Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[2]), 
                            X10Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.10Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[3]), 
                            X15Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.15Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[4]), 
                            X20Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.20Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[5]),
                            X30Hr_Tg_Norm_2 = a2234d_unique_tmt$TMT.30Hr_2 * (mean(a2234d_tg_vec2)/a2234d_tg_vec2[6]),
                            X0Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.0Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[1]), 
                            X05Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.05Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[2]), 
                            X10Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.10Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[3]), 
                            X15Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.15Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[4]), 
                            X20Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.20Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[5]),
                            X30Hr_Tg_Norm_3 = a2234d_unique_tmt$TMT.30Hr_3 * (mean(a2234d_tg_vec3)/a2234d_tg_vec3[6]),
                            X0Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.0Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[1]), 
                            X05Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.05Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[2]), 
                            X10Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.10Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[3]), 
                            X15Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.15Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[4]), 
                            X20Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.20Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[5]),
                            X30Hr_Tg_Norm_4 = a2234d_unique_tmt$TMT.30Hr_4 * (mean(a2234d_tg_vec4)/a2234d_tg_vec4[6]),
                            X0Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.0Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[1]), 
                            X05Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.05Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[2]), 
                            X10Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.10Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[3]), 
                            X15Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.15Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[4]), 
                            X20Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.20Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[5]),
                            X30Hr_Tg_Norm_5 = a2234d_unique_tmt$TMT.30Hr_5 * (mean(a2234d_tg_vec5)/a2234d_tg_vec5[6]),
                            X0Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.0Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[1]), 
                            X05Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.05Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[2]), 
                            X10Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.10Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[3]), 
                            X15Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.15Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[4]), 
                            X20Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.20Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[5]),
                            X30Hr_Tg_Norm_6 = a2234d_unique_tmt$TMT.30Hr_6 * (mean(a2234d_tg_vec6)/a2234d_tg_vec6[6]))

c1264r_unique_tmt <- mutate(c1264r_unique_tmt, 
                            X0Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.0Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[1]), 
                            X05Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.05Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[2]), 
                            X10Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.10Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[3]), 
                            X15Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.15Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[4]), 
                            X20Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.20Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[5]), 
                            X30Hr_Tg_Norm_1 = c1264r_unique_tmt$TMT.30Hr_1 * (mean(c1264r_tg_vec1)/c1264r_tg_vec1[6]),
                            X0Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.0Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[1]), 
                            X05Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.05Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[2]), 
                            X10Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.10Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[3]), 
                            X15Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.15Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[4]), 
                            X20Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.20Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[5]),
                            X30Hr_Tg_Norm_2 = c1264r_unique_tmt$TMT.30Hr_2 * (mean(c1264r_tg_vec2)/c1264r_tg_vec2[6]),
                            X0Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.0Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[1]), 
                            X05Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.05Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[2]), 
                            X10Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.10Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[3]), 
                            X15Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.15Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[4]), 
                            X20Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.20Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[5]),
                            X30Hr_Tg_Norm_3 = c1264r_unique_tmt$TMT.30Hr_3 * (mean(c1264r_tg_vec3)/c1264r_tg_vec3[6]),
                            X0Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.0Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[1]), 
                            X05Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.05Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[2]), 
                            X10Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.10Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[3]), 
                            X15Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.15Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[4]), 
                            X20Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.20Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[5]),
                            X30Hr_Tg_Norm_4 = c1264r_unique_tmt$TMT.30Hr_4 * (mean(c1264r_tg_vec4)/c1264r_tg_vec4[6]),
                            X0Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.0Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[1]), 
                            X05Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.05Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[2]), 
                            X10Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.10Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[3]), 
                            X15Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.15Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[4]), 
                            X20Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.20Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[5]),
                            X30Hr_Tg_Norm_5 = c1264r_unique_tmt$TMT.30Hr_5 * (mean(c1264r_tg_vec5)/c1264r_tg_vec5[6]),
                            X0Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.0Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[1]), 
                            X05Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.05Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[2]), 
                            X10Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.10Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[3]), 
                            X15Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.15Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[4]), 
                            X20Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.20Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[5]),
                            X30Hr_Tg_Norm_6 = c1264r_unique_tmt$TMT.30Hr_6 * (mean(c1264r_tg_vec6)/c1264r_tg_vec6[6]))

ml240_unique_tmt <- mutate(ml240_unique_tmt, 
                           X0Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.0Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[1]), 
                           X05Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.05Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[2]), 
                           X10Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.10Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[3]), 
                           X15Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.15Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[4]), 
                           X20Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.20Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[5]), 
                           X30Hr_Tg_Norm_1 = ml240_unique_tmt$TMT.30Hr_1 * (mean(ml240_tg_vec1)/ml240_tg_vec1[6]),
                           X0Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.0Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[1]), 
                           X05Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.05Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[2]), 
                           X10Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.10Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[3]), 
                           X15Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.15Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[4]), 
                           X20Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.20Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[5]),
                           X30Hr_Tg_Norm_2 = ml240_unique_tmt$TMT.30Hr_2 * (mean(ml240_tg_vec2)/ml240_tg_vec2[6]))

# Calculate ratio enrichment of proteins vs NoHpg and then average
wt_unique_tmt <- mutate(wt_unique_tmt, 
                        X0Hr_ratio_1 = X0Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                        X05Hr_ratio_1 = X05Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                        X10Hr_ratio_1 = X10Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                        X15Hr_ratio_1 = X15Hr_Tg_Norm_1 - TMT.NoHpg_1,
                        X20Hr_ratio_1 = X20Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                        X30Hr_ratio_1 = X30Hr_Tg_Norm_1 - TMT.NoHpg_1,
                        X0Hr_ratio_2 = X0Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                        X05Hr_ratio_2 = X05Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                        X10Hr_ratio_2 = X10Hr_Tg_Norm_2 - TMT.NoHpg_2,
                        X15Hr_ratio_2 = X15Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                        X20Hr_ratio_2 = X20Hr_Tg_Norm_2 - TMT.NoHpg_2,
                        X30Hr_ratio_2 = X30Hr_Tg_Norm_2 - TMT.NoHpg_2,
                        X0Hr_ratio_3 = X0Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                        X05Hr_ratio_3 = X05Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                        X10Hr_ratio_3 = X10Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                        X15Hr_ratio_3 = X15Hr_Tg_Norm_3 - TMT.NoHpg_3,
                        X20Hr_ratio_3 = X20Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                        X30Hr_ratio_3 = X30Hr_Tg_Norm_3 - TMT.NoHpg_3,
                        X0Hr_ratio_4 = X0Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                        X05Hr_ratio_4 = X05Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                        X10Hr_ratio_4 = X10Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                        X15Hr_ratio_4 = X15Hr_Tg_Norm_4 - TMT.NoHpg_4,
                        X20Hr_ratio_4 = X20Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                        X30Hr_ratio_4 = X30Hr_Tg_Norm_4 - TMT.NoHpg_4,
                        X0Hr_ratio_5 = X0Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                        X05Hr_ratio_5 = X05Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                        X10Hr_ratio_5 = X10Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                        X15Hr_ratio_5 = X15Hr_Tg_Norm_5 - TMT.NoHpg_5,
                        X20Hr_ratio_5 = X20Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                        X30Hr_ratio_5 = X30Hr_Tg_Norm_5 - TMT.NoHpg_5,
                        X0Hr_ratio_6 = X0Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                        X05Hr_ratio_6 = X05Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                        X10Hr_ratio_6 = X10Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                        X15Hr_ratio_6 = X15Hr_Tg_Norm_6 - TMT.NoHpg_6,
                        X20Hr_ratio_6 = X20Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                        X30Hr_ratio_6 = X30Hr_Tg_Norm_6 - TMT.NoHpg_6)


a2234d_unique_tmt <- mutate(a2234d_unique_tmt, 
                            X0Hr_ratio_1 = X0Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X05Hr_ratio_1 = X05Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X10Hr_ratio_1 = X10Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X15Hr_ratio_1 = X15Hr_Tg_Norm_1 - TMT.NoHpg_1,
                            X20Hr_ratio_1 = X20Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X30Hr_ratio_1 = X30Hr_Tg_Norm_1 - TMT.NoHpg_1,
                            X0Hr_ratio_2 = X0Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X05Hr_ratio_2 = X05Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X10Hr_ratio_2 = X10Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X15Hr_ratio_2 = X15Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X20Hr_ratio_2 = X20Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X30Hr_ratio_2 = X30Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X0Hr_ratio_3 = X0Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X05Hr_ratio_3 = X05Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X10Hr_ratio_3 = X10Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X15Hr_ratio_3 = X15Hr_Tg_Norm_3 - TMT.NoHpg_3,
                            X20Hr_ratio_3 = X20Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X30Hr_ratio_3 = X30Hr_Tg_Norm_3 - TMT.NoHpg_3,
                            X0Hr_ratio_4 = X0Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X05Hr_ratio_4 = X05Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X10Hr_ratio_4 = X10Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X15Hr_ratio_4 = X15Hr_Tg_Norm_4 - TMT.NoHpg_4,
                            X20Hr_ratio_4 = X20Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X30Hr_ratio_4 = X30Hr_Tg_Norm_4 - TMT.NoHpg_4,
                            X0Hr_ratio_5 = X0Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X05Hr_ratio_5 = X05Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X10Hr_ratio_5 = X10Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X15Hr_ratio_5 = X15Hr_Tg_Norm_5 - TMT.NoHpg_5,
                            X20Hr_ratio_5 = X20Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X30Hr_ratio_5 = X30Hr_Tg_Norm_5 - TMT.NoHpg_5,
                            X0Hr_ratio_6 = X0Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X05Hr_ratio_6 = X05Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X10Hr_ratio_6 = X10Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X15Hr_ratio_6 = X15Hr_Tg_Norm_6 - TMT.NoHpg_6,
                            X20Hr_ratio_6 = X20Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X30Hr_ratio_6 = X30Hr_Tg_Norm_6 - TMT.NoHpg_6)

c1264r_unique_tmt <- mutate(c1264r_unique_tmt, 
                            X0Hr_ratio_1 = X0Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X05Hr_ratio_1 = X05Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X10Hr_ratio_1 = X10Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X15Hr_ratio_1 = X15Hr_Tg_Norm_1 - TMT.NoHpg_1,
                            X20Hr_ratio_1 = X20Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                            X30Hr_ratio_1 = X30Hr_Tg_Norm_1 - TMT.NoHpg_1,
                            X0Hr_ratio_2 = X0Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X05Hr_ratio_2 = X05Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X10Hr_ratio_2 = X10Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X15Hr_ratio_2 = X15Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                            X20Hr_ratio_2 = X20Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X30Hr_ratio_2 = X30Hr_Tg_Norm_2 - TMT.NoHpg_2,
                            X0Hr_ratio_3 = X0Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X05Hr_ratio_3 = X05Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X10Hr_ratio_3 = X10Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X15Hr_ratio_3 = X15Hr_Tg_Norm_3 - TMT.NoHpg_3,
                            X20Hr_ratio_3 = X20Hr_Tg_Norm_3 - TMT.NoHpg_3, 
                            X30Hr_ratio_3 = X30Hr_Tg_Norm_3 - TMT.NoHpg_3,
                            X0Hr_ratio_4 = X0Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X05Hr_ratio_4 = X05Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X10Hr_ratio_4 = X10Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X15Hr_ratio_4 = X15Hr_Tg_Norm_4 - TMT.NoHpg_4,
                            X20Hr_ratio_4 = X20Hr_Tg_Norm_4 - TMT.NoHpg_4, 
                            X30Hr_ratio_4 = X30Hr_Tg_Norm_4 - TMT.NoHpg_4,
                            X0Hr_ratio_5 = X0Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X05Hr_ratio_5 = X05Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X10Hr_ratio_5 = X10Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X15Hr_ratio_5 = X15Hr_Tg_Norm_5 - TMT.NoHpg_5,
                            X20Hr_ratio_5 = X20Hr_Tg_Norm_5 - TMT.NoHpg_5, 
                            X30Hr_ratio_5 = X30Hr_Tg_Norm_5 - TMT.NoHpg_5,
                            X0Hr_ratio_6 = X0Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X05Hr_ratio_6 = X05Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X10Hr_ratio_6 = X10Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X15Hr_ratio_6 = X15Hr_Tg_Norm_6 - TMT.NoHpg_6,
                            X20Hr_ratio_6 = X20Hr_Tg_Norm_6 - TMT.NoHpg_6, 
                            X30Hr_ratio_6 = X30Hr_Tg_Norm_6 - TMT.NoHpg_6)

ml240_unique_tmt <- mutate(ml240_unique_tmt, 
                           X0Hr_ratio_1 = X0Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                           X05Hr_ratio_1 = X05Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                           X10Hr_ratio_1 = X10Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                           X15Hr_ratio_1 = X15Hr_Tg_Norm_1 - TMT.NoHpg_1,
                           X20Hr_ratio_1 = X20Hr_Tg_Norm_1 - TMT.NoHpg_1, 
                           X30Hr_ratio_1 = X30Hr_Tg_Norm_1 - TMT.NoHpg_1,
                           X0Hr_ratio_2 = X0Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                           X05Hr_ratio_2 = X05Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                           X10Hr_ratio_2 = X10Hr_Tg_Norm_2 - TMT.NoHpg_2,
                           X15Hr_ratio_2 = X15Hr_Tg_Norm_2 - TMT.NoHpg_2, 
                           X20Hr_ratio_2 = X20Hr_Tg_Norm_2 - TMT.NoHpg_2,
                           X30Hr_ratio_2 = X30Hr_Tg_Norm_2 - TMT.NoHpg_2)

#Process Averages, SD, and SEMs of Raw FCs

wt_unique_tmt <- mutate(wt_unique_tmt, 
                        X0Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_6")], na.rm=TRUE),
                        X05Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_6")], na.rm=TRUE),
                        X10Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_6")], na.rm=TRUE),
                        X15Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_6")], na.rm=TRUE),
                        X20Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_6")], na.rm=TRUE),
                        X30Hr_ratio_Avg = rowMeans(wt_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_6")], na.rm=TRUE))

a2234d_unique_tmt <- mutate(a2234d_unique_tmt, 
                            X0Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_5","X0Hr_ratio_6")], na.rm=TRUE),
                            X05Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_5","X05Hr_ratio_6")], na.rm=TRUE),
                            X10Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_5","X10Hr_ratio_6")], na.rm=TRUE),
                            X15Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_5","X15Hr_ratio_6")], na.rm=TRUE),
                            X20Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_5","X20Hr_ratio_6")], na.rm=TRUE),
                            X30Hr_ratio_Avg = rowMeans(a2234d_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_5","X30Hr_ratio_6")], na.rm=TRUE))

c1264r_unique_tmt <- mutate(c1264r_unique_tmt, 
                            X0Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_5","X0Hr_ratio_6")], na.rm=TRUE),
                            X05Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_5","X05Hr_ratio_6")], na.rm=TRUE),
                            X10Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_5","X10Hr_ratio_6")], na.rm=TRUE),
                            X15Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_5","X15Hr_ratio_6")], na.rm=TRUE),
                            X20Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_5","X20Hr_ratio_6")], na.rm=TRUE),
                            X30Hr_ratio_Avg = rowMeans(c1264r_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_5","X30Hr_ratio_6")], na.rm=TRUE))

ml240_unique_tmt <- mutate(ml240_unique_tmt, 
                           X0Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2")], na.rm=TRUE),
                           X05Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2")], na.rm=TRUE),
                           X10Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2")], na.rm=TRUE),
                           X15Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2")], na.rm=TRUE),
                           X20Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2")], na.rm=TRUE),
                           X30Hr_ratio_Avg = rowMeans(ml240_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2")], na.rm=TRUE))

wt_unique_tmt <- wt_unique_tmt %>%
  mutate(X0Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_6")]), na.rm=TRUE),
         X05Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_6")]), na.rm=TRUE),
         X10Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_6")]), na.rm=TRUE),
         X15Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_6")]), na.rm=TRUE),
         X20Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_6")]), na.rm=TRUE),
         X30Hr_ratio_sd = rowSds(as.matrix(wt_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_6")]), na.rm=TRUE))


a2234d_unique_tmt <- a2234d_unique_tmt %>%
  mutate(X0Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_5","X0Hr_ratio_6")]), na.rm=TRUE),
         X05Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_5","X05Hr_ratio_6")]), na.rm=TRUE),
         X10Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_5","X10Hr_ratio_6")]), na.rm=TRUE),
         X15Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_5","X15Hr_ratio_6")]), na.rm=TRUE),
         X20Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_5","X20Hr_ratio_6")]), na.rm=TRUE),
         X30Hr_ratio_sd = rowSds(as.matrix(a2234d_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_5","X30Hr_ratio_6")]), na.rm=TRUE))

c1264r_unique_tmt <- c1264r_unique_tmt %>%
  mutate(X0Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2","X0Hr_ratio_3","X0Hr_ratio_4","X0Hr_ratio_5","X0Hr_ratio_6")]), na.rm=TRUE),
         X05Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2","X05Hr_ratio_3","X05Hr_ratio_4","X05Hr_ratio_5","X05Hr_ratio_6")]), na.rm=TRUE),
         X10Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2","X10Hr_ratio_3","X10Hr_ratio_4","X10Hr_ratio_5","X10Hr_ratio_6")]), na.rm=TRUE),
         X15Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2","X15Hr_ratio_3","X15Hr_ratio_4","X15Hr_ratio_5","X15Hr_ratio_6")]), na.rm=TRUE),
         X20Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2","X20Hr_ratio_3","X20Hr_ratio_4","X20Hr_ratio_5","X20Hr_ratio_6")]), na.rm=TRUE),
         X30Hr_ratio_sd = rowSds(as.matrix(c1264r_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2","X30Hr_ratio_3","X30Hr_ratio_4","X30Hr_ratio_5","X30Hr_ratio_6")]), na.rm=TRUE))

ml240_unique_tmt <- ml240_unique_tmt %>%
  mutate(X0Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X0Hr_ratio_1","X0Hr_ratio_2")]), na.rm=TRUE),
         X05Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X05Hr_ratio_1","X05Hr_ratio_2")]), na.rm=TRUE),
         X10Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X10Hr_ratio_1","X10Hr_ratio_2")]), na.rm=TRUE),
         X15Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X15Hr_ratio_1","X15Hr_ratio_2")]), na.rm=TRUE),
         X20Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X20Hr_ratio_1","X20Hr_ratio_2")]), na.rm=TRUE),
         X30Hr_ratio_sd = rowSds(as.matrix(ml240_unique_tmt[,c("X30Hr_ratio_1","X30Hr_ratio_2")]), na.rm=TRUE))

wt_unique_tmt <- wt_unique_tmt %>%
  mutate(X0Hr_ratio_sem = X0Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))),
         X05Hr_ratio_sem = X05Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))),
         X10Hr_ratio_sem = X10Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))),
         X15Hr_ratio_sem = X15Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))),
         X20Hr_ratio_sem = X20Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))),
         X30Hr_ratio_sem = X30Hr_ratio_sd / sqrt(rowSums(!is.na(select(wt_unique_tmt, ends_with("Avg"))))))

a2234d_unique_tmt <- a2234d_unique_tmt %>%
  mutate(X0Hr_ratio_sem = X0Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))),
         X05Hr_ratio_sem = X05Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))),
         X10Hr_ratio_sem = X10Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))),
         X15Hr_ratio_sem = X15Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))),
         X20Hr_ratio_sem = X20Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))),
         X30Hr_ratio_sem = X30Hr_ratio_sd / sqrt(rowSums(!is.na(select(a2234d_unique_tmt, ends_with("Avg"))))))

c1264r_unique_tmt <- c1264r_unique_tmt %>%
  mutate(X0Hr_ratio_sem = X0Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))),
         X05Hr_ratio_sem = X05Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))),
         X10Hr_ratio_sem = X10Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))),
         X15Hr_ratio_sem = X15Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))),
         X20Hr_ratio_sem = X20Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))),
         X30Hr_ratio_sem = X30Hr_ratio_sd / sqrt(rowSums(!is.na(select(c1264r_unique_tmt, ends_with("Avg"))))))

ml240_unique_tmt <- ml240_unique_tmt %>%
  mutate(X0Hr_ratio_sem = X0Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))),
         X05Hr_ratio_sem = X05Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))),
         X10Hr_ratio_sem = X10Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))),
         X15Hr_ratio_sem = X15Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))),
         X20Hr_ratio_sem = X20Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))),
         X30Hr_ratio_sem = X30Hr_ratio_sd / sqrt(rowSums(!is.na(select(ml240_unique_tmt, ends_with("Avg"))))))

#Format SEM values for heatmap
wt_sem <- select(wt_unique_tmt, ends_with("ratio_sem")) 
a2234d_sem <- select(a2234d_unique_tmt, ends_with("ratio_sem")) 
c1264r_sem <- select(c1264r_unique_tmt, contains("ratio_sem")) 

wt_sem_t <- wt_sem %>%
  t() %>%
  as.data.frame()

a2234d_sem_t <- a2234d_sem %>%
  t() %>%
  as.data.frame()

c1264r_sem_t <- c1264r_sem %>%
  t() %>%
  as.data.frame()

colnames(wt_sem_t) <- wt_unique_tmt$name
wt_sem_t <- mutate(wt_sem_t, Sample = colnames(wt_sem))

colnames(a2234d_sem_t) <- a2234d_unique_tmt$name
a2234d_sem_t <- mutate(a2234d_sem_t, Sample = colnames(a2234d_sem))

colnames(c1264r_sem_t) <- c1264r_unique_tmt$name
c1264r_sem_t <- mutate(c1264r_sem_t, Sample = colnames(c1264r_sem))

wt_sem_long <- wt_sem_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "SEM")

a2234d_sem_long <- a2234d_sem_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "SEM")

c1264r_sem_long <- c1264r_sem_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "SEM")

#Filter and plot SEM values based on pathway
wt_tg_sem <- filter(wt_sem_long, Gene == "TG")
wt_ribosome_sem <- filter(wt_sem_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
wt_import_sem <- filter(wt_sem_long, grepl(paste(interactors_import, collapse="|"), Gene))
wt_nglyco_sem <- filter(wt_sem_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
wt_hsp_sem <- filter(wt_sem_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
wt_redox_sem <- filter(wt_sem_long, grepl(paste(interactors_redox, collapse="|"), Gene))
wt_collagen_sem <- filter(wt_sem_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
wt_traffic_sem <- filter(wt_sem_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
wt_ERAD_sem <- filter(wt_sem_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
wt_autophagy_sem <- filter(wt_sem_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
wt_misc_sem <- filter(wt_sem_long, grepl(paste(interactors_misc, collapse="|"), Gene))
wt_total_sem <- rbind(wt_tg_sem, wt_ribosome_sem, wt_import_sem, wt_nglyco_sem, wt_hsp_sem, wt_redox_sem, wt_collagen_sem, wt_traffic_sem, wt_ERAD_sem, wt_autophagy_sem, wt_misc_sem)

a2234d_tg_sem <- filter(a2234d_sem_long, Gene == "TG")
a2234d_ribosome_sem <- filter(a2234d_sem_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
a2234d_import_sem <- filter(a2234d_sem_long, grepl(paste(interactors_import, collapse="|"), Gene))
a2234d_nglyco_sem <- filter(a2234d_sem_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
a2234d_hsp_sem <- filter(a2234d_sem_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
a2234d_redox_sem <- filter(a2234d_sem_long, grepl(paste(interactors_redox, collapse="|"), Gene))
a2234d_collagen_sem <- filter(a2234d_sem_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
a2234d_traffic_sem <- filter(a2234d_sem_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
a2234d_ERAD_sem <- filter(a2234d_sem_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
a2234d_autophagy_sem <- filter(a2234d_sem_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
a2234d_misc_sem <- filter(a2234d_sem_long, grepl(paste(interactors_misc, collapse="|"), Gene))
a2234d_total_sem <- rbind(a2234d_tg_sem, a2234d_ribosome_sem, a2234d_import_sem, a2234d_nglyco_sem, a2234d_hsp_sem, a2234d_redox_sem, a2234d_collagen_sem, a2234d_traffic_sem, a2234d_ERAD_sem, a2234d_autophagy_sem, a2234d_misc_sem)

c1264r_tg_sem <- filter(c1264r_sem_long, Gene == "TG")
c1264r_ribosome_sem <- filter(c1264r_sem_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
c1264r_import_sem <- filter(c1264r_sem_long, grepl(paste(interactors_import, collapse="|"), Gene))
c1264r_nglyco_sem <- filter(c1264r_sem_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
c1264r_hsp_sem <- filter(c1264r_sem_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
c1264r_redox_sem <- filter(c1264r_sem_long, grepl(paste(interactors_redox, collapse="|"), Gene))
c1264r_collagen_sem <- filter(c1264r_sem_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
c1264r_traffic_sem <- filter(c1264r_sem_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
c1264r_ERAD_sem <- filter(c1264r_sem_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
c1264r_autophagy_sem <- filter(c1264r_sem_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
c1264r_misc_sem <- filter(c1264r_sem_long, grepl(paste(interactors_misc, collapse="|"), Gene))
c1264r_total_sem <- rbind(c1264r_tg_sem, c1264r_ribosome_sem, c1264r_import_sem, c1264r_nglyco_sem, c1264r_hsp_sem, c1264r_redox_sem, c1264r_collagen_sem, c1264r_traffic_sem, c1264r_ERAD_sem, c1264r_autophagy_sem, c1264r_misc_sem)

#Plot SEM of ratios by pathways
ggplot(data = wt_total_sem, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_sem", "X05Hr_ratio_sem", "X10Hr_ratio_sem", "X15Hr_ratio_sem", "X20Hr_ratio_sem", "X30Hr_ratio_sem")),
                                          y = factor(Gene, level = rev(ordered_interactors)),  fill = SEM))+
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="gray28", high = "navyblue", midpoint = 0, na.value = "black", breaks = c(0.25, 0.5, 0.75, 1, 1.25)) +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs (−) Hpg - SEM")

ggplot(data = a2234d_total_sem, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_sem", "X05Hr_ratio_sem", "X10Hr_ratio_sem", "X15Hr_ratio_sem", "X20Hr_ratio_sem", "X30Hr_ratio_sem")),
                                              y = factor(Gene, level = rev(ordered_interactors)),  fill = SEM))+
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="gray28", high = "navyblue", midpoint = 0, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs (−) Hpg - SEM")

ggplot(data = c1264r_total_sem, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_sem","X05Hr_ratio_sem", "X10Hr_ratio_sem", "X15Hr_ratio_sem", "X20Hr_ratio_sem", "X30Hr_ratio_sem")),
                                              y = factor(Gene, level = rev(ordered_interactors)),  fill = SEM))+
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="gray28", high = "navyblue", midpoint = 0, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs (−) Hpg - SEM")

#Format ratio averages w/o scaling
wt_unique_raw <- select(wt_unique_tmt, ends_with("ratio_Avg"))
a2234d_unique_raw <- select(a2234d_unique_tmt, ends_with("ratio_Avg"))
c1264r_unique_raw <- select(c1264r_unique_tmt, contains("ratio_Avg"))

c1264r_ind_raw <- select(c1264r_unique_tmt, contains("ratio_"))

wt_unique_raw_t <- wt_unique_raw %>%
  t() %>%
  as.data.frame()

a2234d_unique_raw_t <- a2234d_unique_raw %>%
  t() %>%
  as.data.frame()

c1264r_unique_raw_t <- c1264r_unique_raw %>%
  t() %>%
  as.data.frame()

c1264r_ind_raw_t <- c1264r_ind_raw %>%
  t() %>%
  as.data.frame()

colnames(wt_unique_raw_t) <- wt_unique_tmt$name
wt_unique_raw_t <- mutate(wt_unique_raw_t, Sample = colnames(wt_unique_raw))

colnames(a2234d_unique_raw_t) <- a2234d_unique_tmt$name
a2234d_unique_raw_t <- mutate(a2234d_unique_raw_t, Sample = colnames(a2234d_unique_raw))

colnames(c1264r_unique_raw_t) <- c1264r_unique_tmt$name
c1264r_unique_raw_t <- mutate(c1264r_unique_raw_t, Sample = colnames(c1264r_unique_raw))

wt_unique_raw_long <- wt_unique_raw_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

a2234d_unique_raw_long <- a2234d_unique_raw_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

c1264r_unique_raw_long <- c1264r_unique_raw_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

#Filter by pathway and plot unscaled data
wt_tg_raw <- filter(wt_unique_raw_long, Gene == "TG")
wt_ribosome_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
wt_import_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_import, collapse="|"), Gene))
wt_nglyco_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
wt_hsp_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
wt_redox_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_redox, collapse="|"), Gene))
wt_collagen_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
wt_traffic_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
wt_ERAD_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
wt_autophagy_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
wt_misc_raw <- filter(wt_unique_raw_long, grepl(paste(interactors_misc, collapse="|"), Gene))
wt_total_raw <- rbind(wt_tg_raw, wt_ribosome_raw, wt_import_raw, wt_nglyco_raw, wt_hsp_raw, wt_redox_raw, wt_collagen_raw, wt_traffic_raw, wt_ERAD_raw, wt_autophagy_raw, wt_misc_raw)

a2234d_tg_raw <- filter(a2234d_unique_raw_long, Gene == "TG")
a2234d_ribosome_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
a2234d_import_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_import, collapse="|"), Gene))
a2234d_nglyco_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
a2234d_hsp_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
a2234d_redox_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_redox, collapse="|"), Gene))
a2234d_collagen_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
a2234d_traffic_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
a2234d_ERAD_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
a2234d_autophagy_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
a2234d_misc_raw <- filter(a2234d_unique_raw_long, grepl(paste(interactors_misc, collapse="|"), Gene))
a2234d_total_raw <- rbind(a2234d_tg_raw, a2234d_ribosome_raw, a2234d_import_raw, a2234d_nglyco_raw, a2234d_hsp_raw, a2234d_redox_raw, a2234d_collagen_raw, a2234d_traffic_raw, a2234d_ERAD_raw, a2234d_autophagy_raw, a2234d_misc_raw)

c1264r_tg_raw <- filter(c1264r_unique_raw_long, Gene == "TG")
c1264r_ribosome_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
c1264r_import_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_import, collapse="|"), Gene))
c1264r_nglyco_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
c1264r_hsp_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
c1264r_redox_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_redox, collapse="|"), Gene))
c1264r_collagen_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
c1264r_traffic_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
c1264r_ERAD_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
c1264r_autophagy_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
c1264r_misc_raw <- filter(c1264r_unique_raw_long, grepl(paste(interactors_misc, collapse="|"), Gene))
c1264r_total_raw <- rbind(c1264r_tg_raw, c1264r_ribosome_raw, c1264r_import_raw, c1264r_nglyco_raw, c1264r_hsp_raw, c1264r_redox_raw, c1264r_collagen_raw, c1264r_traffic_raw, c1264r_ERAD_raw, c1264r_autophagy_raw, c1264r_misc_raw)

ggplot(data = wt_total_raw, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                          y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="navyblue", mid="gray95", high = "red4", midpoint = 0, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs No Hpg")

ggplot(data = a2234d_total_raw, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                              y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="navyblue", mid="gray95", high = "red4", midpoint = 0, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs No Hpg")

ggplot(data = c1264r_total_raw, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg","X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                              y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradient2(low ="navyblue", mid="gray95", high = "red4", midpoint = 0, na.value = "black", breaks = c(-4,-3,-2,-1,0,1,2,3)) +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Log2 FC vs No Hpg")

#Pull ratio averages and scale -1 to 0 to 1 on aggregate & format for heatmaps
wt_unique_tmt_fullscale_avg <- wt_unique_tmt %>%
  select(ends_with("ratio_Avg")) %>%
  apply(1, function(x) {ifelse(x < 0, -x/min(x, na.rm = TRUE), x/max(x, na.rm = TRUE))}) %>%
  t()  %>%
  as.data.frame()

a2234d_unique_tmt_fullscale_avg <- a2234d_unique_tmt %>%
  select(ends_with("ratio_Avg")) %>%
  apply(1, function(x) {ifelse(x < 0, -x/min(x, na.rm = TRUE), x/max(x, na.rm = TRUE))}) %>%
  t()  %>%
  as.data.frame()

c1264r_unique_tmt_fullscale_avg <- c1264r_unique_tmt %>%
  select(ends_with("ratio_Avg")) %>%
  apply(1, function(x) {ifelse(x < 0, -x/min(x, na.rm = TRUE), x/max(x, na.rm = TRUE))}) %>%
  t()  %>%
  as.data.frame() 

ml240_unique_tmt_fullscale_avg <- ml240_unique_tmt %>%
  select(contains("ratio_Avg")) %>%
  apply(1, function(x) {ifelse(x < 0, -x/min(x, na.rm = TRUE), x/max(x, na.rm = TRUE))}) %>%
  t()  %>%
  as.data.frame()

wt_unique_tmt_fullscale_avg_t <- wt_unique_tmt_fullscale_avg %>%
  t() %>%
  as.data.frame()

a2234d_unique_tmt_fullscale_avg_t <- a2234d_unique_tmt_fullscale_avg %>%
  t() %>%
  as.data.frame()

c1264r_unique_tmt_fullscale_avg_t <- c1264r_unique_tmt_fullscale_avg %>%
  t() %>%
  as.data.frame()

ml240_unique_tmt_fullscale_avg_t <- ml240_unique_tmt_fullscale_avg %>%
  t() %>%
  as.data.frame()

colnames(wt_unique_tmt_fullscale_avg_t) <- wt_unique_tmt$name
wt_unique_tmt_fullscale_avg_t <- mutate(wt_unique_tmt_fullscale_avg_t, Sample = rownames(wt_unique_tmt_fullscale_avg_t))

colnames(a2234d_unique_tmt_fullscale_avg_t) <- a2234d_unique_tmt$name
a2234d_unique_tmt_fullscale_avg_t <- mutate(a2234d_unique_tmt_fullscale_avg_t, Sample = rownames(a2234d_unique_tmt_fullscale_avg_t))

colnames(c1264r_unique_tmt_fullscale_avg_t) <- c1264r_unique_tmt$name
c1264r_unique_tmt_fullscale_avg_t <- mutate(c1264r_unique_tmt_fullscale_avg_t, Sample = colnames(c1264r_unique_tmt_fullscale_avg))

colnames(ml240_unique_tmt_fullscale_avg_t) <- ml240_unique_tmt$name
ml240_unique_tmt_fullscale_avg_t <- mutate(ml240_unique_tmt_fullscale_avg_t, Sample = colnames(ml240_unique_tmt_fullscale_avg))

wt_unique_tmt_fullscale_avg_long <- wt_unique_tmt_fullscale_avg_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

a2234d_unique_tmt_fullscale_avg_long <- a2234d_unique_tmt_fullscale_avg_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

c1264r_unique_tmt_fullscale_avg_long <- c1264r_unique_tmt_fullscale_avg_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

ml240_unique_tmt_fullscale_avg_long <- ml240_unique_tmt_fullscale_avg_t %>%
  pivot_longer(cols = -matches("Sample"), names_to = "Gene", values_to = "Enrichment")

wt_tg_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, Gene == "TG")
wt_ribosome_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
wt_import_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_import, collapse="|"), Gene))
wt_nglyco_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
wt_hsp_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
wt_redox_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_redox, collapse="|"), Gene))
wt_collagen_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
wt_traffic_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
wt_ERAD_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
wt_autophagy_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
wt_misc_fullscale_avg <- filter(wt_unique_tmt_fullscale_avg_long, grepl(paste(interactors_misc, collapse="|"), Gene))
wt_total_fullscale_avg <- rbind(wt_tg_fullscale_avg, wt_ribosome_fullscale_avg, wt_import_fullscale_avg,  wt_nglyco_fullscale_avg, wt_hsp_fullscale_avg, 
                                wt_redox_fullscale_avg, wt_collagen_fullscale_avg, wt_traffic_fullscale_avg, wt_ERAD_fullscale_avg, 
                                wt_autophagy_fullscale_avg, wt_misc_fullscale_avg)

a2234d_tg_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, Gene == "TG")
a2234d_ribosome_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
a2234d_import_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_import, collapse="|"), Gene))
a2234d_nglyco_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
a2234d_hsp_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
a2234d_redox_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_redox, collapse="|"), Gene))
a2234d_collagen_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
a2234d_traffic_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
a2234d_ERAD_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
a2234d_autophagy_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
a2234d_misc_fullscale_avg <- filter(a2234d_unique_tmt_fullscale_avg_long, grepl(paste(interactors_misc, collapse="|"), Gene))
a2234d_total_fullscale_avg <- rbind(a2234d_tg_fullscale_avg, a2234d_ribosome_fullscale_avg, a2234d_import_fullscale_avg,  a2234d_nglyco_fullscale_avg, a2234d_hsp_fullscale_avg, 
                                    a2234d_redox_fullscale_avg, a2234d_collagen_fullscale_avg, a2234d_traffic_fullscale_avg, a2234d_ERAD_fullscale_avg, 
                                    a2234d_autophagy_fullscale_avg, a2234d_misc_fullscale_avg)

c1264r_tg_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, Gene == "TG")
c1264r_ribosome_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
c1264r_import_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_import, collapse="|"), Gene))
c1264r_nglyco_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
c1264r_hsp_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
c1264r_redox_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_redox, collapse="|"), Gene))
c1264r_collagen_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
c1264r_traffic_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
c1264r_ERAD_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
c1264r_autophagy_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
c1264r_misc_fullscale_avg <- filter(c1264r_unique_tmt_fullscale_avg_long, grepl(paste(interactors_misc, collapse="|"), Gene))
c1264r_total_fullscale_avg <- rbind(c1264r_tg_fullscale_avg, c1264r_ribosome_fullscale_avg, c1264r_import_fullscale_avg,  c1264r_nglyco_fullscale_avg, c1264r_hsp_fullscale_avg, 
                                    c1264r_redox_fullscale_avg, c1264r_collagen_fullscale_avg, c1264r_traffic_fullscale_avg, c1264r_ERAD_fullscale_avg, 
                                    c1264r_autophagy_fullscale_avg, c1264r_misc_fullscale_avg)

ml240_tg_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, Gene == "TG")
ml240_ribosome_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ribosome, collapse="|"), Gene))
ml240_import_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_import, collapse="|"), Gene))
ml240_nglyco_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_NGlyco, collapse="|"), Gene))
ml240_hsp_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_Hsp, collapse="|"), Gene))
ml240_redox_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_redox, collapse="|"), Gene))
ml240_collagen_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_collagen, collapse="|"), Gene))
ml240_traffic_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_traffic, collapse="|"), Gene))
ml240_ERAD_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_ERAD, collapse="|"), Gene))
ml240_autophagy_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_autophagy, collapse="|"), Gene))
ml240_misc_fullscale_avg <- filter(ml240_unique_tmt_fullscale_avg_long, grepl(paste(interactors_misc, collapse="|"), Gene))
ml240_total_fullscale_avg <- rbind(ml240_tg_fullscale_avg, ml240_ribosome_fullscale_avg, ml240_import_fullscale_avg,  ml240_nglyco_fullscale_avg, ml240_hsp_fullscale_avg, 
                                   ml240_redox_fullscale_avg, ml240_collagen_fullscale_avg, ml240_traffic_fullscale_avg, ml240_ERAD_fullscale_avg, 
                                   ml240_autophagy_fullscale_avg, ml240_misc_fullscale_avg)


# Define palette 
colors_scaled <- c(seq(-1,0, length=10),0, seq(0,0.33,length = 30), seq(0.33, 0.66, length =30), seq(0.66,1, length = 30))
my_palette <- colorRampPalette(c("gray25","gray50","gray75","gray100", "#4B2991", "#EA4F88", "#EDD9A3"))


#Plot heatmaps
ggplot(data = wt_total_fullscale_avg, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                                    y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradientn(colors = my_palette(99), breaks = colors_scaled, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg")

ggplot(data = a2234d_total_fullscale_avg, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                                        y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradientn(colors = my_palette(99), breaks = colors_scaled2, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg")

ggplot(data = c1264r_total_fullscale_avg, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                                        y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradientn(colors = my_palette(99), breaks = colors_scaled, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg")

ggplot(data = ml240_total_fullscale_avg, mapping = aes(x = factor(Sample, level = c("X0Hr_ratio_Avg", "X05Hr_ratio_Avg", "X10Hr_ratio_Avg", "X15Hr_ratio_Avg", "X20Hr_ratio_Avg", "X30Hr_ratio_Avg")),
                                                       y = factor(Gene, level = rev(ordered_interactors)), fill = Enrichment)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_gradientn(colors = my_palette(99), breaks = colors_scaled, na.value = "black") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg")
