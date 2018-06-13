

```R
library(phyloseq)
library(ggplot2)
library(vegan)
library(readr)
library(data.table)
library(dplyr)
library(picante)
library(car)
library(dunn.test)
library(ggpubr)
```

### Pre-processing


```R
#CREATE PHYLOSEQ OBJECT

biom <- import_biom("~/Data/data.biom")
map <- import_qiime_sample_data("~/Data/metadata.txt")
tree <- read_tree("~/Data/phylo.tre")
batrun_nt <- merge_phyloseq(biom,tree,map)
batrun_nt <- prune_samples((sample_sums(batrun_nt) > 1000), batrun_nt)

#CHANGE TAXA COLUMN NAMES FROM (E.G. FROM"RANK" TO "KINGDOM")
colnames(tax_table(batrun_nt)) <- c("kingdom", "phylum", "class", "order", "family",  "genus", "species")
```

#### Summarize sequencing depth by sample type


```R
batrun_sdt = data.table(as(sample_data(batrun_nt), "data.frame"),
                 TotalReads = sample_sums(batrun_nt), keep.rownames = TRUE)
setnames(batrun_sdt, "rn", "SampleName")
pSeqDepth = ggplot(batrun_sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.





![png](output_4_2.png)



```R
batrun.summary = batrun_sdt %>%
    group_by(SampleType) %>%
summarise(avg = mean(TotalReads), 
              sd = sd(TotalReads), 
              sem = sd/sqrt(length(TotalReads))) %>%
    arrange((SampleType))
batrun.summary[complete.cases(batrun.summary),]
```


<table>
<thead><tr><th scope=col>SampleType</th><th scope=col>avg</th><th scope=col>sd</th><th scope=col>sem</th></tr></thead>
<tbody>
	<tr><td>FecalFTA  </td><td>41963.77  </td><td>21697.85  </td><td>1082.1905 </td></tr>
	<tr><td>TongueEtOH</td><td>28650.17  </td><td>22365.17  </td><td>1153.3956 </td></tr>
	<tr><td>SkinEtOH  </td><td>26992.98  </td><td>15483.25  </td><td> 723.4847 </td></tr>
</tbody>
</table>



### Alpha diversity

#### Estimate Shannon diversity index


```R
#sample_tab = read.csv("~/Data/metadata.csv")
sample_tab = read.csv("/group/gilbert-lab/Lutz/Batmicrobiome/batrun_finaldata_for_phyloseq_Oct2017/Metadata_files/FinalMappingFile_noswabs.csv")
colnames(sample_tab)[1] <- "SampleID"

```


```R
#ESTIMATE SHANNON DIVERSITY
sha = estimate_richness(batrun_nt, split = TRUE, measures = "Shannon")
sha$SampleID = rownames(sha)

dfsha = merge(sha, sample_tab, on="SampleID")

ordsha = sha %>% group_by("SampleType") %>%
    summarize(med.sha = median(Shannon)) %>%
    arrange(desc(med.sha)) %>%
    magrittr::extract2("SampleType")

stats = dfsha %>%
    group_by(SampleType) %>%
    summarize(median = median(Shannon), 
              mean = mean(Shannon),
              sd = sd(Shannon),
              sem = sd(Shannon/sqrt(length(Shannon)))) %>%
    arrange(-desc(SampleType))
stats
```


<table>
<thead><tr><th scope=col>SampleType</th><th scope=col>median</th><th scope=col>mean</th><th scope=col>sd</th><th scope=col>sem</th></tr></thead>
<tbody>
	<tr><td>FecalFTA  </td><td>1.411972  </td><td>1.524695  </td><td>0.8681219 </td><td>0.04181596</td></tr>
	<tr><td>SkinEtOH  </td><td>4.563169  </td><td>4.423806  </td><td>0.9183671 </td><td>0.04263412</td></tr>
	<tr><td>TongueEtOH</td><td>1.395439  </td><td>1.463988  </td><td>0.9202091 </td><td>0.04647758</td></tr>
</tbody>
</table>




```R
#CHECK FOR NORMAL/NON-NORMAL DISTRIBUTION
p = ggqqplot(dfsha, x = "Shannon")
p = p + facet_wrap("SampleType")
p

#CHECK FOR SIGNIFICANCE
stats.krusk = kruskal.test(Shannon~SampleType, data=dfsha)
stats.krusk

#POST-HOC ANALYSIS
dunn = dunn.test(dfsha$Shannon, g=dfsha$SampleType, method="bonferroni")
dunn
```




    
    	Kruskal-Wallis rank sum test
    
    data:  Shannon by SampleType
    Kruskal-Wallis chi-squared = 796.02, df = 2, p-value < 2.2e-16



      Kruskal-Wallis rank sum test
    
    data: x and group
    Kruskal-Wallis chi-squared = 796.0165, df = 2, p-value = 0
    
    
                               Comparison of x by group                            
                                     (Bonferroni)                                  
    Col Mean-|
    Row Mean |   FecalFTA   SkinEtOH
    ---------+----------------------
    SkinEtOH |  -24.09624
             |    0.0000*
             |
    TongueEt |   0.761793   24.27289
             |     0.6693    0.0000*
    
    alpha = 0.05
    Reject Ho if p <= alpha/2



<dl>
	<dt>$chi2</dt>
		<dd>796.016530788115</dd>
	<dt>$Z</dt>
		<dd><ol class=list-inline>
	<li>-24.0962458584237</li>
	<li>0.761793527707069</li>
	<li>24.272899664308</li>
</ol>
</dd>
	<dt>$P</dt>
		<dd><ol class=list-inline>
	<li>1.36841093298805e-128</li>
	<li>0.223091621948333</li>
	<li>1.89504741884996e-130</li>
</ol>
</dd>
	<dt>$P.adjusted</dt>
		<dd><ol class=list-inline>
	<li>4.10523279896414e-128</li>
	<li>0.669274865844999</li>
	<li>5.68514225654989e-130</li>
</ol>
</dd>
	<dt>$comparisons</dt>
		<dd><ol class=list-inline>
	<li>'FecalFTA - SkinEtOH'</li>
	<li>'FecalFTA - TongueEtOH'</li>
	<li>'SkinEtOH - TongueEtOH'</li>
</ol>
</dd>
</dl>




![png](output_10_4.png)



```R
#PLOT (DENSITY PLOT) SHANNON DIVERSITY BY SAMPLE TYPE
dfsha$SampleType <- factor(
     dfsha$SampleType, 
     labels=c('Gut','Skin', 'Tongue'))
p2 = ggdensity(dfsha, x = "Shannon", color = "SampleType", fill= "SampleType", rug=TRUE)
p2
```




![png](output_11_1.png)


#### Estimate Observed ESV richness


```R
sample_tab = read.csv("~/Data/metadata.csv")
```


```R
#ESTIMATE OBSERVED ESV RICHNESS
obs = estimate_richness(batrun_nt, split = TRUE, measures = "Observed")
obs$SampleID = rownames(obs)

dfobs = merge(obs, sample_tab, on="SampleID")

ordobs = obs %>% group_by("SampleType") %>%
    summarize(med.obs = median(Observed)) %>%
    arrange(desc(med.obs)) %>%
    magrittr::extract2("SampleType")

stats = dfobs %>%
    group_by(SampleType) %>%
    summarize(median = median(Observed), 
              mean = mean(Observed),
              sd = sd(Observed),
              sem = sd(Observed/sqrt(length(Observed)))) %>%
    arrange(-desc(SampleType))
stats
```


<table>
<thead><tr><th scope=col>SampleType</th><th scope=col>median</th><th scope=col>mean</th><th scope=col>sd</th><th scope=col>sem</th></tr></thead>
<tbody>
	<tr><td>FecalFTA  </td><td> 49       </td><td> 68.80510 </td><td> 68.71589 </td><td> 3.309928 </td></tr>
	<tr><td>SkinEtOH  </td><td>413       </td><td>553.16164 </td><td>388.47839 </td><td>18.034657 </td></tr>
	<tr><td>TongueEtOH</td><td> 50       </td><td> 77.91837 </td><td> 92.12656 </td><td> 4.653094 </td></tr>
</tbody>
</table>




```R
#CHECK FOR NORMAL/NON-NORMAL DISTRIBUTION
p3 = ggqqplot(dfobs, x = "Observed")
p3 = p3 + facet_wrap("SampleType")
p3

#CHECK FOR SIGNIFICANCE
stats.krusk = kruskal.test(Observed~SampleType, data=dfobs)
stats.krusk

#POST-HOC ANALYSIS
dunn = dunn.test(dfobs$Observed, g=dfobs$SampleType, method="bonferroni")
dunn

#PLOT (BOXPLOT) OBSERVED ESV RICHNESS BY SAMPLE TYPE
p4 = ggplot(na.omit(dfobs), aes(x = SampleType, y = Observed)) + 
        ggtitle("ESV Richness by Sample Type")

p4 = p4 + geom_boxplot() + ylab("Observed ESVs")

p4 = p4 + geom_point(size = .75, aes(color = SampleType),  
                   position = position_jitter(h = 0, w = 0.2))

p4 = p4 + scale_color_manual(values=c("#9999CC","#66CC99","#CC6666"))

p4 = p4 + theme_bw() + theme(text=element_text(size=12, color="black",family="Arial"),
                           legend.key = element_blank(),
                           legend.position = "none",
                           axis.title.y=element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())
p4
```




    
    	Kruskal-Wallis rank sum test
    
    data:  Observed by SampleType
    Kruskal-Wallis chi-squared = 780.3, df = 2, p-value < 2.2e-16



      Kruskal-Wallis rank sum test
    
    data: x and group
    Kruskal-Wallis chi-squared = 780.2963, df = 2, p-value = 0
    
    
                               Comparison of x by group                            
                                     (Bonferroni)                                  
    Col Mean-|
    Row Mean |   FecalFTA   SkinEtOH
    ---------+----------------------
    SkinEtOH |  -24.10057
             |    0.0000*
             |
    TongueEt |   0.279667   23.78661
             |     1.0000    0.0000*
    
    alpha = 0.05
    Reject Ho if p <= alpha/2



<dl>
	<dt>$chi2</dt>
		<dd>780.296278447862</dd>
	<dt>$Z</dt>
		<dd><ol class=list-inline>
	<li>-24.1005720699503</li>
	<li>0.279667159291494</li>
	<li>23.7866124989105</li>
</ol>
</dd>
	<dt>$P</dt>
		<dd><ol class=list-inline>
	<li>1.23271177702263e-128</li>
	<li>0.389866438181682</li>
	<li>2.29756168634141e-125</li>
</ol>
</dd>
	<dt>$P.adjusted</dt>
		<dd><ol class=list-inline>
	<li>3.69813533106789e-128</li>
	<li>1</li>
	<li>6.89268505902424e-125</li>
</ol>
</dd>
	<dt>$comparisons</dt>
		<dd><ol class=list-inline>
	<li>'FecalFTA - SkinEtOH'</li>
	<li>'FecalFTA - TongueEtOH'</li>
	<li>'SkinEtOH - TongueEtOH'</li>
</ol>
</dd>
</dl>






![png](output_15_5.png)



![png](output_15_6.png)


#### Examine ESV richness by elevation


```R
#PLOT ESV RICHNESS BY ELEVATION
p5 = ggplot(dfobs, aes(x=Elevation, y=Observed)) + geom_point() + 
        ggtitle("Observed ESV Richness by Elevation")
p5 = p5 + geom_smooth(method = "lm", se=TRUE, alpha=.2)
p5 = p5 + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
p5 = p5 + facet_wrap(~SampleType)
p5
```




![png](output_17_1.png)



```R
#SUMMARIZE REGRESSION COEFFICIENTS - Break down by sample type before running
fit <- lm(Elevation ~ Observed, data = <data>)
summary(fit)
```

### Beta diversity

#### Transform data


```R
batrun <- batrun_nt %>%
    transform_sample_counts(function(x) x / sum(x))
```

#### Subset transformed data by sample type (for beta diversity metrics)


```R
#GUT
batrun_gut = subset_samples(batrun, SampleType=="FecalFTA")
batrun_gut = prune_taxa(taxa_sums(batrun_gut)>0, batrun_gut)
#ORAL
batrun_oral = subset_samples(batrun, SampleType=="TongueEtOH")
batrun_oral = prune_taxa(taxa_sums(batrun_oral)>0, batrun_oral)
#SKIN
batrun_skin = subset_samples(batrun, SampleType=="SkinEtOH")
batrun_skin = prune_taxa(taxa_sums(batrun_skin)>0, batrun_skin)
```

#### Calculate beta diversity based on various distance metrics


```R
# WEIGHTED UNIFRAC
batgut.dist.wunif = distance(batrun_gut, method = "wunifrac", type = "samples")
batoral.dist.wunif = distance(batrun_oral, method = "wunifrac", type = "samples")
batskin.dist.wunif = distance(batrun_skin, method = "wunifrac", type = "samples")

# UNWEIGHTED UNIFRAC
batgut.dist.unif = distance(batrun_gut, method = "unifrac", type = "samples")
batoral.dist.unif = distance(batrun_oral, method = "unifrac", type = "samples")
batskin.dist.unif = distance(batrun_skin, method = "unifrac", type = "samples")

# BRAY-CURTIS
batgut.dist.bc = distance(batrun_gut, method = "bray", type = "samples")
batoral.dist.bc = distance(batrun_oral, method = "bray", type = "samples")
batskin.dist.bc = distance(batrun_skin, method = "bray", type = "samples")
```

#### Perform ADONIS test(s) to check for significance of variables


```R
#GUT
batgut.wunif.adtest = adonis2(batgut.dist.wunif ~ Host_Genus/Host_species + Ecto_status + Haem_status, strata= Locality, by = "margin", data=batgut.df)
batgut.wunif.adtest
#ORAL
batoral.wunif.adtest = adonis2(batoral.dist.wunif ~ Host_Genus/Host_species + Ecto_status + Haem_status, strata= Locality, by = "margin", data=batoral.df)
batoral.wunif.adtest
#SKIN
batskin.wunif.adtest = adonis2(batskin.dist.wunif ~ Host_Genus/Host_species + Ecto_status + Haem_status, strata= Locality, by = "margin", data=batskin.df)
batskin.wunif.adtest
```

#### Create PCoA plots


```R
#E.G. CREATE PCOA PLOT FOR GUT MICROBIOME / WEIGHTED UNIFRAC
batgut.wunif.log = ordinate(batrun_gut, method = "PCoA", distance = "wunifrac")
batgut.wunif.evals = batgut.wunif.log$values$Eigenvalues
batgut.wunif.pcoa = plot_ordination(batsall_.00, batsall.00.wunif.log, color = "Locality") + ggtitle("Bat Gut Microbiome - Weighted Unifrac PCoA") + geom_point(size=2) + coord_fixed(sqrt(batgut.wunif.evals[2] / batgut.wunif.evals[1]))
batgut.wunif.pcoa
```

### Plotting Rank Abundance

#### Convert phyloseq objects to data.frame


```R
batgut.mdf = psmelt(batrun_gut)
batoral.mdf = psmelt(batrun_oral)
batskin.mdf = psmelt(batrun_skin)
```

#### Establish top six most abundant phyla for each anatomical site (can plot any rank)


```R
#GUT

N = 6

topN = batgut.mdf %>%
    group_by(Rank2) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank2) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank2")

mdfgut.other = batgut.mdf %>%
    mutate(Rank2.other = ifelse(Rank2 %in% topN, as.character(Rank2), "Other")) %>%
    group_by(Sample, Rank2.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdfgut.summary = mdfgut.other %>%
    group_by(Rank2.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank2.other))
```

    Selecting by m



```R
#ORAL

N = 6

topN = batoral.mdf %>%
    group_by(Rank2) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank2) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank2")

mdforal.other = batoral.mdf %>%
    mutate(Rank2.other = ifelse(Rank2 %in% topN, as.character(Rank2), "Other")) %>%
    group_by(Sample, Rank2.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdforal.summary = mdforal.other %>%
    group_by(Rank2.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank2.other))
```

    Selecting by m



```R
#SKIN

N = 6

topN = batskin.mdf %>%
    group_by(Rank2) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank2) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank2")

mdfskin.other = batskin.mdf %>%
    mutate(Rank2.other = ifelse(Rank2 %in% topN, as.character(Rank2), "Other")) %>%
    group_by(Sample, Rank2.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdfskin.summary = mdfskin.other %>%
    group_by(Rank2.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank2.other))
```

    Selecting by m



```R
#ADD SAMPLETYPE COLUMN
mdfgut.other$SampleType <- "gut"
mdforal.other$SampleType <- "oral"
mdfskin.other$SampleType <- "skin"

#MERGE RESULTS FOR PLOTTING
mdfall.other = rbind(mdfgut.other, mdforal.other, mdfskin.other)
```

#### Create stacked barcharts faceted by anatomical site


```R
mdfall.other$SampleType <- factor(
     mdfall.other$SampleType, 
     labels=c('Gut','Oral','Skin'))

mdfall.other$Rank2.other <- factor(
     mdfall.other$Rank2.other, 
     labels=c('Other   ', 'Actinobacteria   ', 'Bacteriodetes   ', 'Chlamydiae   ', 'Cyanobacteria   ', 'Euryarcheota   ', 'Firmicutes   ', 'Fusobacteria   ', 'Proteobacteria   ', 'Tenericutes   '))


p6 = ggplot(mdfall.other[order(mdfall.other$Rank2.other,decreasing=T),], aes(x=Sample, y = Abundance, fill = Rank2.other))
                        
p6 = p6 + facet_wrap(~SampleType, shrink = TRUE, scales = "free_x", ncol = 10) #labeller = as_labeller("Phylum" = top))

p6 = p6 + ylab("Relative Abundance") 

p6 = p6 + geom_bar(stat="identity")

#p6 = p6 + scale_fill_brewer(palette = "Spectral")

p6 = p6 + scale_y_continuous(expand = c(0,0)) + scale_fill_brewer(palette = "Spectral", na.value= "#0E024C") 

p6 = p6 + guides(fill = guide_legend(title = "Phylum", nrow = 4, keywidth = 1, keyheight = 1))

p6 = p6 + theme_bw() + theme(panel.grid.major = element_blank(),
              axis.title.y=element_text(margin=margin(0,20,0,0)),
              text=element_text(size=10, color="black",family="Arial"),
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              legend.key.size = unit(1, "mm"),
              legend.spacing.x = unit(1,"in"),
              axis.title.x=element_blank())
              

p6 = p6 + theme(legend.position="bottom") 

p6 = p6 + theme(panel.spacing = unit(.5, "lines"))

p6 = p6 + ggtitle("Bat Microbiome - Top 6 Phyla by Anatomical Site") + theme(plot.title = element_text(hjust = 0.5))

p6
```




![png](output_39_1.png)


#### Plot ten most abundant bacterial families by host family


```R
#GUT

N = 10

topN = batgut.mdf %>%
    group_by(Rank5) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank5) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank5")

mdfgut.other = batgut.mdf %>%
    mutate(Rank5.other = ifelse(Rank5 %in% topN, as.character(Rank5), "Other")) %>%
    group_by(Host_Family,Rank5.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdfgut.summary = mdfgut.other %>%
    group_by(Host_Family,Rank5.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank5.other))
```

    Selecting by m



```R
#PLOT RELATIVE ABUNDANCE OF BACTERIAL FAMILY FACETED BY HOST FAMILY

p7 <- ggplot(mdfgut.other, aes(x="", y=Abundance, fill=Rank5.other)) +
geom_bar(stat = "identity", width = 1, position = "fill") +
coord_polar("y", start=0) + 
labs(title="Gut Microbiome - Top 10 Families by Host Family")

p7 = p7 + theme_minimal() + theme(legend.position="bottom", 
                                  plot.title = element_text(hjust = 0.5),
                                  axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  panel.grid  = element_blank())

p7  = p7 + facet_wrap( ~Host_Family, shrink = TRUE) + 
        scale_fill_brewer(palette="Spectral", na.value="#0E024C")


p7
```




![png](output_42_1.png)



```R
#ORAL

N = 10

topN = batoral.mdf %>%
    group_by(Rank5) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank5) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank5")

mdforal.other = batoral.mdf %>%
    mutate(Rank5.other = ifelse(Rank5 %in% topN, as.character(Rank5), "Other")) %>%
    group_by(Host_Family, Rank5.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdforal.summary = mdforal.other %>%
    group_by(Host_Family,Rank5.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank5.other))
```

    Selecting by m



```R
#PLOT RELATIVE ABUNDANCE OF BACTERIAL FAMILY FACETED BY HOST FAMILY

p8 <- ggplot(mdforal.other, aes(x="", y=Abundance, fill=Rank5.other)) +
geom_bar(stat = "identity", width = 1, position = "fill") +
coord_polar("y", start=0) + 
labs(title="Oral Microbiome - Top 10 Families by Host Family")

p8 = p8 + theme_minimal() + theme(legend.position="bottom", 
                                  plot.title = element_text(hjust = 0.5),
                                  axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  panel.grid  = element_blank())

p8  = p8 + facet_wrap( ~Host_Family, shrink = TRUE) + 
        scale_fill_brewer(palette="Spectral", na.value="#0E024C")


p8
```




![png](output_44_1.png)



```R
#SKIN

N = 10

topN = batskin.mdf %>%
    group_by(Rank5) %>%
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Rank5) %>%
    summarize(m = mean(Abundance)) %>%
    top_n(N) %>%
    arrange(desc(m)) %>%
    magrittr::extract2("Rank5")

mdfskin.other = batskin.mdf %>%
    mutate(Rank5.other = ifelse(Rank5 %in% topN, as.character(Rank5), "Other")) %>%
    group_by(Host_Family, Rank5.other) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Abundance)

mdfskin.summary = mdfskin.other %>%
    group_by(Host_Family,Rank5.other) %>%
    summarise(avg = mean(Abundance), 
              sd = sd(Abundance), 
              sem = sd/sqrt(length(Abundance)),
              proportion = round(avg*100,2)) %>%
    arrange((Rank5.other))
```

    Selecting by m



```R
#PLOT RELATIVE ABUNDANCE OF BACTERIAL FAMILY FACETED BY HOST FAMILY

p9 <- ggplot(mdfskin.other, aes(x="", y=Abundance, fill=Rank5.other)) +
geom_bar(stat = "identity", width = 1, position = "fill") +
coord_polar("y", start=0) + 
labs(title="Skin Microbiome - Top 10 Families by Host Family")

p9 = p9 + theme_minimal() + theme(legend.position="bottom", 
                                  plot.title = element_text(hjust = 0.5),
                                  axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  panel.grid  = element_blank())

p9  = p9 + facet_wrap( ~Host_Family, shrink = TRUE) + 
        scale_fill_brewer(palette="Spectral", na.value="#0E024C")


p9
```




![png](output_46_1.png)

