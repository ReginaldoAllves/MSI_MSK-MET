rm(list = ls())  ## clean the environment

### Libraries

library(readxl)        # read excel files easily
library(tidyverse)     # data wrangling
library(ggplot2)
library(ggpubr)

################################################################################
############################# FUNCTIONS ########################################
################################################################################

### Generates density plots for metastasis frequency among MSI groups

density.plot = function(x.axis, x.lab){
  plot = ggdensity(df.2, x = x.axis,
            add = "mean", rug = TRUE,
            color = "msi", fill = "msi",
            palette = c("#00AFBB", "#E7B800","#E7B801")) +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = x.lab) +
    theme (text=element_text(size = 13),
           axis.text.x=element_text(colour="black", size = 11),
           axis.text.y=element_text(colour="black", size = 11)) + 
    theme(axis.line = element_line(size = 1, linetype = "solid"),
          axis.text = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12))
  return(plot)
}

### Merges metastasis burden and metastasis number from different analyses


merge.plots = function(plot.A, plot.B) {
     ggarrange(
        plot.A, plot.B, labels = c("A", "B"), 
        common.legend = TRUE, legend = "bottom")
}

### Generates correlation plots for different metastasis against TMB

correlation.plot = function(axis.x, axis.y, tumor_type, m.x, m.y, x.title){
  df.2 %>% 
    ggscatter(x = axis.x, y = axis.y, 
              add = "reg.line", conf.int = TRUE,
              title = paste(axis.x, " VS ", axis.y, "(", tumor_type, ")", sep = " "),
              xlab = x.title, ylab = "Metastasis (n)",
              shape = 21, fill = "lightgray", color = "black", size = 2.2) +
    stat_cor(method = "spearman", label.x = 1, label.y = m.y) +
    coord_cartesian(ylim = c(0,m.y), xlim = c(0,m.x)) +
    theme (text=element_text(size = 13),
           axis.text.x=element_text(colour="black", size = 11),
           axis.text.y=element_text(colour="black", size = 11)) + 
    theme(axis.line = element_line(size = 1, linetype = "solid"),
          axis.text = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 14))
}

### Generates correlation plots for MSIscores against TMB

corr.tmb.msi = function(tumor_type, m.x, m.y){
  df.2 %>% 
    ggscatter(x = "tmb", y = "msi_score", 
            add = "reg.line", conf.int = TRUE,
            title = paste("TMB vs MSI score (", tumor_type, ")", sep = ""),
            xlab = "TMB", ylab = "MSI score",
            shape = 21, fill = "lightgray", color = "black", size = 2.2) +
  stat_cor(method = "spearman", label.x = 1, label.y = m.y - 1) +
  coord_cartesian(ylim = c(0,m.y), xlim = c(0,m.x)) +
  theme (text=element_text(size = 13),
         axis.text.x=element_text(colour="black", size = 11),
         axis.text.y=element_text(colour="black", size = 11)) + 
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 14))
}

### Generates violin plots for Metastasis metrics between MSI and MSS tumors

violin.plot = function(axis.x, axis.y, title.x, title.y, m.y){
   ggviolin(df.2, x = axis.x, y = axis.y, fill = axis.x,
                    palette = c("#00AFBB", "#E7B800", "#E7B801"),
                    add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(method = "wilcox.test", label.y = m.y + 5) + # Adds significance levels
    scale_y_continuous(name = title.y) +
    scale_x_discrete(name = title.x) +
    theme (text=element_text(size = 13),
           axis.text.x=element_text(colour="black", size = 11),
           axis.text.y=element_text(colour="black", size = 11)) + 
    theme(axis.line = element_line(size = 1, linetype = "solid"),
          axis.text = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 14))
}

### Generates violin plots for Metastasis metrics between MSI and MSS tumors

box.plot = function(axis.x, axis.y, title.x, title.y, m.y){
  ggboxplot(df.2, x = axis.x, y = axis.y,
            color = axis.x, palette =c("#00AFBB", "#E7B800", "#E7B801"),
            add = "jitter", shape = axis.x) + 
    stat_compare_means(method = "wilcox.test", label.y = m.y + 5) + # Add significance levels
    stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
    scale_y_continuous(name = title.y) +
    scale_x_discrete(name = title.x) +
    theme (text=element_text(size = 13),
           axis.text.x=element_text(colour="black", size = 11),
           axis.text.y=element_text(colour="black", size = 11)) + 
    theme(axis.line = element_line(size = 1, linetype = "solid"),
          axis.text = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 14)) 
}

####################################################################
############################# MAIN #################################
####################################################################

### Data Wrangling

# Reading table

df = read_excel("input/MSK-MET.Suppl.Table.1.xlsx")

# Generates a data frame with variables to be used (table.desc)

table.desc = df %>% select(curated_organ_system,
                           cancer_type,
                           oncotree_code,curated_subtype_display,
                           sample_type,
                           met_count,
                           met_site_count,
                           msi_score,
                           msi_type,
                           tmb)

# Renaming variables

table.desc = table.desc %>% rename(system = 1,
                                   tumor_type = 2,
                                   code = 3,
                                   subtype = 4,
                                   number_metastasis = 5,
                                   metastasis_burden = 6,
                                   msi_status=8)

# Creating variables with classification of MSI status

table.desc = table.desc %>% mutate(msi_status = replace(msi_status, 
                                                        msi_score <= 10, "MSS"),
                                   msi_status = replace(msi_status, 
                                                        msi_score > 10, "MSI"))

# Creating variables with classification of TMB status

table.desc = table.desc %>%
  mutate(tmb_status = tmb) %>%
  mutate(tmb_status = replace(tmb_status, tmb <= 10, "low"),
         tmb_status = replace(tmb_status, tmb > 10, "high"))

# Creating a table with tumor code info for searching when needed

table.code = df %>% select(cancer_type,
                           oncotree_code,
                           curated_subtype_display)

table.code = unique(table.code)

write.table(table.code, "output/Table.Code.xls", 
            row.names = FALSE,
            quote = FALSE, 
            sep = "\t")

# Grouping variables for statistic analyses

table.desc.2 = table.desc %>% 
  select(3,6,7,8,9,10,11)

grouped.table.1 = table.desc.2 %>% group_by(code) ## to obtain frequencies of each tumor code
grouped.table.2 = table.desc.2 %>% group_by(code, msi_status) ## to obtain frequencies of each code, by MSI status
grouped.table.3 = table.desc.2 %>% group_by(code, tmb_status) ## to obtain frequencies of each code, by TMB status

# Creating tables with descriptive statistics

description.1 = grouped.table.1 %>% summarize(
  observations = n(),
  mean_msi_score = mean(msi_score),
  st_msi = sd(msi_score),
  mean_tmb = mean(tmb),
  sd_tmb = sd(tmb)
) %>% ungroup () %>% droplevels(.)

description.2 = grouped.table.2 %>% summarize(
  observations = n(),
  mean_msi_score = mean(msi_score),
  sd_msi = sd(msi_score),
) %>% ungroup () %>% droplevels(.)

description.3 = grouped.table.3 %>% summarize(
  observations = n(),
  mean_tmb = mean(tmb),
  st_tmb = sd(tmb)
) %>% ungroup () %>% droplevels(.)

# Frequency of MSS and MSI tumors

n.msi.mss = description.2 %>% filter(msi_status == "MSS") %>%
  select(1,3) %>%
  rename(n.mss = 2)

final.table = description.1 %>%
  left_join(n.msi.mss, by = "code") %>%
  mutate(n.msi = observations - n.mss) %>%
  relocate(observations, .before = n.mss)

# Frequency of TMB-high and TMB-low tumors

n.tmb.low = description.3 %>% filter(tmb_status == "low") %>%
  select(1,3) %>%
  rename(n.low = 2)

final.table = final.table %>%
  left_join(n.tmb.low, by = "code") %>%
  mutate(n.high = observations - n.low)


# Creating columns to save statistics values from Wilcox test

final.table = final.table %>%
  mutate(Wilcox.MSI.burden = 0) %>%
  mutate(Wilcox.TMB.burden = 0)

final.table = final.table %>%
  mutate(Spearman.MSI.burden = 0) %>%
  mutate(Spearman.TMB.burden = 0)

# Vector with tumor codes

tumor.type = unique(final.table$code)

# Iteration to save p-values from Wilcox test

for(i in 1:length(tumor.type)){
  if(tumor.type[i] == final.table$code[i]){
    met.number.msi = table.desc %>% 
      filter(code == tumor.type[i] & msi_status == "MSI") %>%
      pull(met_site_count)
  
    met.number.mss = table.desc %>% 
      filter(code == tumor.type[i] & msi_status == "MSS") %>%
      pull(met_site_count)
  
    met.number.high = table.desc %>% 
      filter(code == tumor.type[i] & tmb_status == "high") %>%
      pull(met_site_count)
  
    met.number.low = table.desc %>% 
      filter(code == tumor.type[i] & tmb_status == "low") %>%
      pull(met_site_count)
    
    if(length(met.number.msi) > 0){
      wt.msi = wilcox.test(met.number.msi, met.number.mss)
      final.table$Wilcox.MSI.burden[i] = wt.msi$p.value
    } else{
      final.table$Wilcox.MSI.burden[i] = NA
    }
    
    if(length(met.number.high > 0)){
      wt.tmb = wilcox.test(met.number.high, met.number.low)
      final.table$Wilcox.TMB.burden[i] = wt.tmb$p.value
    } else{
      final.table$Wilcox.TMB.burden[i] = NA
    }
  }
  else{
    final.table$Wilcox.MSI.burden[i] = NA
    final.table$Wilcox.TMB.burden[i] = NA
  }
}

# Spearman Correlation

final.table = final.table %>%
  rename(Spearman.MSI.burden.r = Spearman.MSI.burden,
         Spearman.MSI.burden.p = Spearman.TMB.burden) %>%
  mutate(Spearman.TMB.burden.r = 0,
         Spearman.TMB.burden.p = 0)

for(i in 1:length(tumor.type)){
  if(tumor.type[i] == final.table$code[i]){
    data = table.desc %>% 
      filter(code == tumor.type[i]) %>%
      select(metastasis_burden, msi_score, tmb)
    if(!is.na(final.table$Wilcox.MSI.burden[i])){
      correlation.msi = cor.test(data$metastasis_burden, data$msi_score, 
                                 method = "spearman", use = "complete.obs")
      final.table$Spearman.MSI.burden.r[i] = correlation.msi$estimate
      final.table$Spearman.MSI.burden.p[i] = correlation.msi$p.value
    }else{
      final.table$Spearman.MSI.burden.r[i] = NA
      final.table$Spearman.MSI.burden.p[i] = NA
    }
    
    if(!is.na(final.table$Wilcox.TMB.burden[i])){
      correlation.tmb = cor.test(data$metastasis_burden, data$tmb, 
                                 method = "spearman", use = "complete.obs")
      final.table$Spearman.TMB.burden.r[i] = correlation.tmb$estimate
      final.table$Spearman.TMB.burden.p[i] = correlation.tmb$p.value
    }else{
      final.table$Spearman.TMB.burden.r[i] = NA
      final.table$Spearman.TMB.burden.p[i] = NA
    }
  } else{
    final.table$Spearman.MSI.burden.r[i] = NA
    final.table$Spearman.MSI.burden.p[i] = NA
    final.table$Spearman.TMB.burden.r[i] = NA
    final.table$Spearman.TMB.burden.p[i] = NA
  }
}

# Writing the final table with statistics 

write.table(final.table, "output/Statistics.values.xls", 
            row.names = FALSE,
            quote = FALSE, 
            sep = "\t")


#### Plots

# Tumor types

tumors = unique(df$oncotree_code)

# Loop to generate and save graphs per tumor type

for(i in 1:length(tumors)){

  df.2 = df %>%
    filter(oncotree_code == tumors[i]) %>% 
    select(sample = patient_id,
           organ = cancer_type,
           tumor = oncotree_code,
           meta_burden = met_site_count,
           meta_number = met_count,
           msi = msi_type,
           msi_score = msi_score,
           tmb = tmb,
           type = sample_type)
  
  df.2 = df.2 %>%
    mutate(msi = replace(msi, msi_score <= 10, "MSS"),
           msi = replace(msi, msi_score > 10, "MSI"))
  
  df.2$msi = as.factor(df.2$msi)
  
  df.2 = df.2 %>%
    mutate(tmb.2 = tmb) %>%
    mutate(tmb.2 = replace(tmb.2, tmb <= 10, "low"),
           tmb.2 = replace(tmb.2, tmb > 10, "high"))
  
  df.2$tmb.2 = as.factor(df.2$tmb.2)
  
  # Creating an output directory
  
  directory = paste("Output/", tumors[i], sep = "")
  dir.create(directory)
  
  # Normality
  
  qplot = paste(directory,"/norm.qqplot.tmb.tiff", sep = "")
  tiff(filename = qplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggqqplot(df.2$tmb))
  dev.off()
  
  qplot = paste(directory,"/norm.qqplot.meta_burden.tiff", sep = "")
  tiff(filename = qplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggqqplot(df.2$meta_burden))
  dev.off()
  
  qplot = paste(directory,"/norm.qqplot.meta_number.tiff", sep = "")
  tiff(filename = qplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggqqplot(df.2$meta_number))
  dev.off()
  
  denplot = paste(directory,"/norm.density.plot.tmb.tiff", sep = "")
  tiff(filename = denplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggdensity(df.2$tmb))
  dev.off()
  
  denplot = paste(directory,"/norm.density.plot.meta.number.tiff", sep = "")
  tiff(filename = denplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggdensity(df.2$meta_number))
  dev.off()
  
  denplot = paste(directory,"/norm.density.plot.meta.burden.tiff", sep = "")
  tiff(filename = denplot, width = 846, height = 703, res = 300, pointsize = 12)
    print(ggdensity(df.2$meta_burden))
  dev.off()

  ## Removing MSS Hypermutated tumors
  
  df.2 = filter(df.2, !(tmb >= 10 & msi == "MSS"))
  
  ###  Density Plot

  meta.number = density.plot(x.axis = "meta_number", 
                             x.lab = "Number of metastatic tumors per patient")
  meta.burden = density.plot(x.axis = "meta_burden", 
                             x.lab = "Metastasis Burden")
  merged = merge.plots(plot.A = meta.number, plot.B = meta.burden)
  
  title.plot = paste(directory, "/Density.Plot.tiff",sep = "")
  tiff(filename = title.plot, width = 1464, height = 763, res = 150, pointsize = 12)
    print(merge.plots(plot.A = meta.number, plot.B = meta.burden))
  dev.off()
  
  
  ## Correlation
  
  max.msi = round(max(df.2$msi_score)) + 2
  max.met = round(max(df.2$meta_number)) + 2
  max.burden = round(max(df.2$meta_burden)) + 2
  max.tmb = round(max(df.2$tmb)) + 2
  
  cor.met.number = correlation.plot(axis.x = "msi_score", axis.y = "meta_number",
                                    tumor_type = tumors[i], 
                                    m.x = max.msi, 
                                    m.y = max.met, 
                                    x.title = "MSI score")
  
  cor.met.burden = correlation.plot(axis.x = "msi_score", axis.y = "meta_burden",
                                    tumor_type = tumors[i], 
                                    m.x = max.msi, 
                                    m.y = max.burden, 
                                    x.title = "MSI score")
  
  cor.met.number.tmb = correlation.plot(axis.x = "tmb", 
                                        axis.y = "meta_number",
                                        tumor_type = tumors[i], 
                                        m.x = max.tmb, 
                                        m.y = max.met, 
                                        x.title = "TMB")
  
  cor.met.burden.tmb = correlation.plot(axis.x = "tmb", 
                                        axis.y = "meta_burden",
                                        tumor_type = tumors[i], 
                                        m.x = max.tmb, 
                                        m.y = max.burden, 
                                        x.title = "TMB")
  
  title.plot = paste(directory, "/Correlation.MSI.Metastasis.tiff",sep = "")
  tiff(filename = title.plot, width = 1464, height = 763, res = 150, pointsize = 12)
    print(merge.plots(plot.A = cor.met.burden, plot.B = cor.met.number))
  dev.off()
  
  title.plot = paste(directory, "/Correlation.TMB.Metastasis.tiff",sep = "")
  tiff(filename = title.plot, width = 1464, height = 763, res = 150, pointsize = 12)
    print(merge.plots(plot.A = cor.met.burden.tmb, plot.B = cor.met.number.tmb))
  dev.off()
  
  
   title.plot = paste(directory, "/Correlation.MSI.TMB.tiff",sep = "")
   tiff(filename = title.plot, width = 1200, height = 720, res = 200, pointsize = 12)
     print(corr.tmb.msi(tumor_type = tumors[i], m.x = max.tmb, m.y = max.msi))
   dev.off()

  # Violin Plot

  violin.met.burden = violin.plot(axis.x = "msi", axis.y = "meta_burden", 
                                  title.x = "", title.y = "Metastasis burden",
                                   m.y = max.burden)
  
  violin.met.number = violin.plot(axis.x = "msi", axis.y = "meta_number", 
                                  title.x = "", title.y = "Number of metastatic tumors per patient",
                                  m.y = max.met)
  
  violin.met.burden.tmb = violin.plot(axis.x = "tmb.2", axis.y = "meta_burden", 
                                  title.x = "TMB", title.y = "Metastasis burden",
                                  m.y = max.burden)
  
  violin.met.number.tmb = violin.plot(axis.x = "tmb.2", axis.y = "meta_number", 
                                  title.x = "TMB", title.y = "Number of metastatic tumors per patient",
                                  m.y = max.met)
  
  title.plot = paste(directory, "/Violin.MSI.met.tiff",sep = "")
  tiff(filename = title.plot, width = 1500, height = 850, res = 150, pointsize = 12)
    print(merge.plots(plot.A = violin.met.burden, plot.B = violin.met.number))
  dev.off()
  
  title.plot = paste(directory, "/Violin.TMB.met.tiff",sep = "")
  tiff(filename = title.plot, width = 1500, height = 850, res = 150, pointsize = 12)
    print(merge.plots(plot.A = violin.met.burden.tmb, plot.B = violin.met.number.tmb))
  dev.off()

  # BoxPlot
  
  boxplot.met.burden = box.plot(axis.x = "msi", axis.y = "meta_burden", 
                                  title.x = "", title.y = "Metastasis burden",
                                  m.y = max.burden)
  
  boxplot.met.number = box.plot(axis.x = "msi", axis.y = "meta_number", 
                                  title.x = "", title.y = "Number of metastatic tumors per patient",
                                  m.y = max.met)
  
  boxplot.met.burden.tmb = box.plot(axis.x = "tmb.2", axis.y = "meta_burden", 
                                      title.x = "TMB", title.y = "Metastasis burden",
                                      m.y = max.burden)
  
  boxplot.met.number.tmb = box.plot(axis.x = "tmb.2", axis.y = "meta_number", 
                                      title.x = "TMB", title.y = "Number of metastatic tumors per patient",
                                      m.y = max.met)
  
  title.plot = paste(directory, "/Boxplot.MSI.met.tiff",sep = "")
  tiff(filename = title.plot, width = 1500, height = 850, res = 150, pointsize = 12)
    print(merge.plots(plot.A = boxplot.met.burden, plot.B = boxplot.met.number))
  dev.off()
  
  title.plot = paste(directory, "/Boxplot.TMB.met.tiff",sep = "")
  tiff(filename = title.plot, width = 1500, height = 850, res = 150, pointsize = 12)
    print(merge.plots(plot.A = boxplot.met.burden.tmb, plot.B = boxplot.met.number.tmb))
  dev.off()
}
