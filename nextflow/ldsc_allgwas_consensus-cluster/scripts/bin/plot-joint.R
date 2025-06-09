#!/usr/bin/env Rscript

library(ggplot2)
library(glue)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
library(tidyr)

args = commandArgs(TRUE)

savep = function(p, f, h=8, w=4){
    png(f, height=h, width=w, units="in", res = 150)
    print(p)
    dev.off()
}
  

d = read.csv(args[1], sep='\t', header=T, as.is=T)
d$cluster = gsub(".broadPeaks_fdr1perc_noBlacklist", "", d$cluster)
d$cluster = gsub("_broadpeaks", "", d$cluster)
d = d[! d$cluster %in% c('static_annotations.Conserved.LindbladToh', 'static_annotations.Coding.UCSC'),]

head(d)

height = nrow(d) * 0.5
maxval = max(d$Coefficient*1e+07 + 1.96*d$Coefficient_std_error*1e+07)
if (maxval < 20 ){
    maxx = maxval + 1
} else {
    maxx = 20
}

d$min_coeff = d$Coefficient*1e+07 - 1.96*d$Coefficient_std_error*1e+07
d$coeff_sig = ifelse( (d$Coefficient - 1.96*d$Coefficient_std_error) > 0 | (d$Coefficient + 1.96*d$Coefficient_std_error) < 0, "s", "ns"  )

p = ggplot(d, aes(y=reorder(cluster, Coefficient), x=Coefficient*1e+07)) +
    geom_point(aes(color = coeff_sig), size=2.5) +
    geom_errorbarh(aes(xmax=Coefficient*1e+07 + 1.96*Coefficient_std_error*1e+07,
                       xmin=Coefficient*1e+07 - 1.96*Coefficient_std_error*1e+07)) +
    facet_wrap(~trait) +
    coord_cartesian(xlim=c(-10, maxx)) +
                                        #xlim(-10, 10) +
    geom_vline(xintercept = 0)

savep(p, glue("{args[2]}.coeff.png"), height, 10)

d$enrich_sig = ifelse( (d$Enrichment - 1.96*d$Enrichment_std_error) > 1 | (d$Enrichment + 1.96*d$Enrichment_std_error) < 1, "s", "ns"  )
p = ggplot(d, aes(y=reorder(cluster, Enrichment), x=Enrichment)) +
    geom_point(aes(color = enrich_sig), size=2.5) +
    geom_errorbarh(aes(xmax=Enrichment + 1.96*Enrichment_std_error,
                       xmin=Enrichment - 1.96*Enrichment_std_error)) +
    facet_wrap(~trait) +
    geom_vline(xintercept = 1)

savep(p, glue("{args[2]}.enrichment.png"), height, 10)


