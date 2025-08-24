library(ggplot2)
library(dplyr)
library(tidyverse)

library(extrafont)
pdf.options(family = "Arial")

setwd("/data/yangyu/Project/skinAtlas/70samples/Result/4.composition_ratio/plot_final/")
# female ####
female_sig <- read_csv('../Female_ALL_cell.Percentage_Ratio_ave.csv')

female_files <- list.files(path = '../', pattern = "^Female_.*\\_ave.csv$")
female_files <- female_files[-1]
data_list <- list()
for (file in female_files) {
  file <- paste0('../',file)
  tmp <- read_csv(file)
  data_list[[file]] <- tmp
}
female_df <- bind_rows(data_list)

## female NF ####
female_sig_NF <- female_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NF') %>%
  #arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
female_sig_NF_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
                         "Spinous.Ker","Spinous.Ker.ADGRL3",
                         "Mel_1","Mel_2",
                         "Papillary.Fb","Mesenchymal.Fb.INHBA","Mesenchymal.Fb.COCH",
                         "Pc1","arteriole.ECs","LEC1",
                         "M2","cDC2","LC2","MigLC",
                         "helper T cells" ,
                         "cytotoxic T cells","ILC3"
) 
female_df_NF <- female_df %>% filter(celltype_Granular %in% female_sig_NF_order) %>% 
  filter(sample_group == 'NF') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NF %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NF_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NF_sig_subtype_boxplot.pdf',dpi = 1200,w=6,h=5)  



female_sig_NF <- female_sig %>% filter(p_value <= 0.01) %>% 
  filter(sample_group == 'NF') %>%
  #arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
female_sig_NF_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
                         "Spinous.Ker","Spinous.Ker.ADGRL3",
                         "Mel_1","Mel_2",
                         "Papillary.Fb","Mesenchymal.Fb.INHBA","Mesenchymal.Fb.COCH",
                         "Pc1","arteriole.ECs","LEC1",
                         "cDC2",#"LC2","MigLC",
                         "helper T cells" ,
                         "cytotoxic T cells"
) 
female_df_NF <- female_df %>% filter(celltype_Granular %in% female_sig_NF_order) %>% 
  filter(sample_group == 'NF') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NF %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NF_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NF_sig_P0.01_subtype_boxplot.pdf',dpi = 1200,w=6,h=5)  

female_sig_NF <- female_sig %>% filter(corrected_p_value <= 0.05) %>% 
  filter(sample_group == 'NF') %>%
  #arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
female_sig_NF_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
                         "Spinous.Ker","Spinous.Ker.ADGRL3",
                         "Mel_1","Mel_2",
                         "Papillary.Fb","Mesenchymal.Fb.INHBA","Mesenchymal.Fb.COCH",
                         #"Pc1",
                         "arteriole.ECs",#"LEC1",
                         #"M1",
                        # "M2",#"cDC1","MigDC","Mast cells",
                        # "cDC2",#"LC2","MigLC",
                         "helper T cells" ,
                        "cytotoxic T cells" #,#"ILC3" 
) 
female_df_NF <- female_df %>% filter(celltype_Granular %in% female_sig_NF) %>% 
  filter(sample_group == 'NF') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NF %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NF_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NF_sig_adjP_subtype_boxplot.pdf',dpi = 1200,w=6,h=4)  


## female NP ####
female_sig_NP <- female_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NP') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
# female_sig_NP_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
#                          "Spinous.Ker.ADGRL3", "Granular.Ker",
#                          "KC.Mel",
#                          "Mesenchymal.Fb.FGFBP2","Mesenchymal.Fb.INHBA",
#                          "Mesenchymal.Fb.COCH",
#                          "Pc1","Pc2","vSMC1","vSMC2",
#                          "arteriole.ECs","venule.ECs","capillary.ECs",
#                          "M2",
#                          "cDC1","MigDC","LC2",
#                          "helper T cells","cytotoxic T cells","NKT cells")
female_sig_NP_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
                         "Spinous.Ker.ADGRL3", "Granular.Ker",
                         "Mesenchymal.Fb.FGFBP2","Mesenchymal.Fb.INHBA",
                         "Pc1","Pc2","vSMC1","vSMC2",
                         "arteriole.ECs","venule.ECs","capillary.ECs","EMT.ECs",
                         "M2",
                         "LC2",
                         "helper T cells","cytotoxic T cells")
female_df_NP <- female_df %>% filter(celltype_Granular %in% female_sig_NP) %>% 
  filter(sample_group == 'NP') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NP %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NP_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NP_sig_subtype_boxplot.pdf',dpi = 1200,w=6,h=4)  


female_sig_NP <- female_sig %>% filter(p_value <= 0.01) %>% 
  filter(sample_group == 'NP') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
# female_sig_NP_order <- c("Basal.Ker.COL17A1","Basal.Ker.ADGRL3","Proliferating.Ker",
#                          "Spinous.Ker.ADGRL3", "Granular.Ker",
#                          "KC.Mel",
#                          "Mesenchymal.Fb.FGFBP2","Mesenchymal.Fb.INHBA",
#                          "Mesenchymal.Fb.COCH",
#                          "Pc1","Pc2","vSMC1","vSMC2",
#                          "arteriole.ECs","venule.ECs","capillary.ECs",
#                          "M2",
#                          "cDC1","MigDC","LC2",
#                          "helper T cells","cytotoxic T cells","NKT cells")
female_sig_NP_order <- c("Basal.Ker.COL17A1",
                         "Spinous.Ker.ADGRL3", "Granular.Ker",
                         "Mesenchymal.Fb.FGFBP2","Mesenchymal.Fb.INHBA",
                         "Pc2","vSMC1","vSMC2",
                         "venule.ECs","capillary.ECs",
                         "M2")
female_df_NP <- female_df %>% filter(celltype_Granular %in% female_sig_NP) %>% 
  filter(sample_group == 'NP') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NP %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NP_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NP_sig_P0.01.subtype_boxplot.pdf',dpi = 1200,w=6,h=4) 


## female NE ####
female_sig_NE <- female_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NE') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
female_sig_NE_order <- c("Basal.Ker.ADGRL3","Proliferating.Ker","Spinous.Ker","Granular.Ker",
                         "Channel-gap_1",
                         "HF-PCs", "Mesenchymal.Fb","Mesenchymal.Fb.COCH","capillary.ECs"
)
female_df_NE <- female_df %>% filter(celltype_Granular %in% female_sig_NE) %>% 
  filter(sample_group == 'NE') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NE %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NE_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NE_sig_subtype_boxplot.pdf',dpi = 1200,w=5,h=3)  
## female NT ####
female_sig_NT <- female_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NT') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
female_sig_NT_order <- c( "Proliferating.Ker","Spinous.Ker",
                          "Papillary.Fb","Mesenchymal.Fb.COCH","Mel_1",
                          "Mel_2",
                          
                          "Pc1","vSMC1","vSMC2",
                         "cDC2"
)
female_df_NT <- female_df %>% filter(celltype_Granular %in% female_sig_NT) %>% 
  filter(sample_group == 'NT') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
female_df_NT %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(female_sig_NT_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('NT_sig_subtype_boxplot.pdf',dpi = 1200,w=5,h=3)  

# Male ####
male_sig <- read_csv('../Male_ALL_cell.Percentage_Ratio_ave.csv')

male_files <- list.files(path = '../', pattern = "^Male_.*\\_ave.csv$")
male_files <- male_files[-1]
data_list <- list()
for (file in male_files) {
  file <- paste0('../',file)
  tmp <- read_csv(file)
  data_list[[file]] <- tmp
}
male_df <- bind_rows(data_list)

## male PS ####
male_sig_PS <- male_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'PS') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
male_sig_PS_order <- c("Basal.Ker.ADGRL3",
                       "Proliferating.Ker",
                       "Mel_2",
                       "Mesenchymal.Fb","Mesenchymal.Fb.INHBA",
                       "vSMC2","post-capillary-venule.ECs",
                       "M2" ,"cDC2", 
                       "helper T cells","cytotoxic T cells","ILC3")
male_df_PS <- male_df %>% filter(celltype_Granular %in% male_sig_PS) %>% 
  filter(sample_group == 'PS') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
male_df_PS %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(male_sig_PS_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('PS_sig_subtype_boxplot.pdf',dpi = 1200,w=6,h=4) 


male_sig_PS <- male_sig %>% filter(p_value <= 0.01) %>% 
  filter(sample_group == 'PS') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
male_sig_PS_order <- c(
                       "Proliferating.Ker",
                       
                       "Mesenchymal.Fb","Mesenchymal.Fb.INHBA",
                       "vSMC2","post-capillary-venule.ECs",
                       "M2" ,"cDC2", 
                       "helper T cells","cytotoxic T cells")
male_df_PS <- male_df %>% filter(celltype_Granular %in% male_sig_PS) %>% 
  filter(sample_group == 'PS') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
male_df_PS %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(male_sig_PS_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('PS_sig_P0.01_subtype_boxplot.pdf',dpi = 1200,w=6,h=4) 
## male NE ####
male_sig_NE <- male_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NE') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
male_sig_NE_order <- c("Proliferating.Ker","Mesenchymal.Fb",
                       "post-capillary-venule.ECs",
                       "EMT.ECs",
                       "cDC2",
                       "ILC3")
male_df_NE <- male_df %>% filter(celltype_Granular %in% male_sig_NE) %>% 
  filter(sample_group == 'NE') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
male_df_NE %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(male_sig_NE_order)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('male_NE_sig_subtype_boxplot.pdf',dpi = 1200,w=6,h=4) 

## male NT ####
male_sig_NT <- male_sig %>% filter(p_value <= 0.05) %>% 
  filter(sample_group == 'NT') %>%
  arrange(desc(Log_Percentage_Ratio_Ave))%>%
  pull(celltype_Granular)
male_sig_NT_order <- c( "arteriole.ECs",
                        "vSMC2","cDC2",
                        "helper T cells","cytotoxic T cells" )
male_df_NT <- male_df %>% filter(celltype_Granular %in% male_sig_NT) %>% 
  filter(sample_group == 'NT') %>%
  mutate(color = ifelse(Log_Percentage_Ratio_Ave > 0, "increase", "decrease")) #%>% distinct()
male_df_NT %>% 
  ggplot(aes(x = log2(Percentage_Ratio), y = celltype_Granular)) +
  geom_boxplot(aes(fill = color),outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("increase" = "red", "decrease" = "blue"))  +
  scale_y_discrete(limits = rev(male_sig_NT)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) + 
  theme_classic()+
  labs(x = "Log2(Percentage Ratio)")
ggsave('male_NT_sig_subtype_boxplot.pdf',dpi = 1200,w=6,h=4) 
