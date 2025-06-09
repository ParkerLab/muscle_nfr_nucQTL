#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(glue))
ggplot2::theme_set(theme_cowplot())
suppressMessages(library(data.table))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
suppressMessages(library(ggrepel))

args = commandArgs(TRUE)

prefix = args[2]

savep = function(p, f, h=8, w=4){
    print(p)
    png(glue("{f}.png"), height=h, width=w, units="in", res = 150)
    print(p)
    dev.off()
}

d = read.csv(args[1], sep='\t', header=T, as.is=T)
d = d[d$Category == "L2_0",]
d$annot = gsub("proinsulin.annot.", "", d$trait)
d$annot = gsub(".results", "", d$annot)
d$trait = "proinsulin"

d$annot = gsub(".broadPeaks_fdr1perc_noBlacklist", "", d$annot)
d$annot = gsub(".33samp_broadpeaks", ".33samp", d$annot)

# Coding and Conserved were in the baseline duh
d = d[! d$annot %in% c('static_annotations.Conserved.LindbladToh', 'static_annotations.Coding.UCSC'),]

head(d)


d$min = d$Coefficient*1e+07 - 1.96*d$Coefficient_std_error*1e+07
p = ggplot(d, aes(y=reorder(annot, Coefficient), x=Coefficient*1e+07)) +
geom_point() +
geom_errorbarh(aes(xmax=Coefficient*1e+07+1.96*Coefficient_std_error*1e+07,
                   xmin=Coefficient*1e+07-1.96*Coefficient_std_error*1e+07)) +
facet_wrap(~trait) +
coord_cartesian(xlim=c(-10,40)) +
#xlim(-10, 10) +
geom_vline(xintercept = 0) + 
theme_bw() +
theme(panel.background=element_rect(fill="white", colour="black")) 
savep(p, glue("{prefix}.single_annot_coeff"), 12, 8)

d$min = d$Enrichment - 1.96*d$Enrichment_std_error
p = ggplot(d, aes(y=reorder(annot, Enrichment), x=Enrichment)) +
geom_point() +
geom_errorbarh(aes(xmax=Enrichment + 1.96*Enrichment_std_error,
                   xmin=Enrichment - 1.96*Enrichment_std_error)) +
facet_wrap(~trait) +
#coord_cartesian(xlim=c(-10,40)) +
#xlim(-10, 10) +
geom_vline(xintercept = 1) +
theme_bw() +
theme(panel.background=element_rect(fill="white", colour="black"))
savep(p, glue("{prefix}.single_annot_enrich"), 12, 8)

#SUBSET plots
keep = read.csv("/lab/work/arushiv/proinsulin/subset_plot", sep='\t', header=T, as.is=T)$Annotation
d = d[d$annot %in% keep,]

d$min = d$Coefficient*1e+07 - 1.96*d$Coefficient_std_error*1e+07
p = ggplot(d, aes(y=reorder(annot, Coefficient), x=Coefficient*1e+07)) +
geom_point() +
geom_errorbarh(aes(xmax=Coefficient*1e+07+1.96*Coefficient_std_error*1e+07,
                   xmin=Coefficient*1e+07-1.96*Coefficient_std_error*1e+07)) +
facet_wrap(~trait) +
coord_cartesian(xlim=c(-10,40)) +
#xlim(-10, 10) +
geom_vline(xintercept = 0) + 
theme_bw() +
theme(panel.background=element_rect(fill="white", colour="black")) 
savep(p, glue("{prefix}.single_annot_coeff_subset"), 8, 8)

d$min = d$Enrichment - 1.96*d$Enrichment_std_error
p = ggplot(d, aes(y=reorder(annot, Enrichment), x=Enrichment)) +
geom_point() +
geom_errorbarh(aes(xmax=Enrichment + 1.96*Enrichment_std_error,
                   xmin=Enrichment - 1.96*Enrichment_std_error)) +
facet_wrap(~trait) +
#coord_cartesian(xlim=c(-10,40)) +
#xlim(-10, 10) +
geom_vline(xintercept = 1) +
theme_bw() +
theme(panel.background=element_rect(fill="white", colour="black"))
savep(p, glue("{prefix}.single_annot_enrich_subset"), 8, 8)


