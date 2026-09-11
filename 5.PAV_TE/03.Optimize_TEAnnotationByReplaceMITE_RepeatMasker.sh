#!/usr/bin/bash

bedtools intersect -a Sitalica_Yugu1.EDTA.TEanno.bed -b Sitalica_Yugu1.genomeV2.RMSK.Harbinger.200-400.bed -v -wa > Sitalica_Yugu1.EDTA.tmp.bed
cat Sitalica_Yugu1.EDTA.tmp.bed Sitalica_Yugu1.genomeV2.RMSK.Harbinger.200-400.bed > Sitalica_Yugu1.EDTA.Replaced.bed 
