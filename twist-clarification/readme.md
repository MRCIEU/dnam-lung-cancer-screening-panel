From Twist:

> Using our default filters it looks like we would be able to cover 80.39% 
> of your targets, with the following regions filtered out: 
> [bed file](twist-clarification/Target_bases_not_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_1X_MTE-93452736_hg19_230901090027.bed)
>
> We can recover a few more targets if we allow shifting of probes but the following would still remain uncovered: 
> [bed file](twist-clarification/all_target_segments_not_covered_by_probes_withshifting.bed)

We assessed the sites the could not be reliably targeted (without 'shifting'). 
Below we refer to these as 'uncovered'.

```
Rscript src/clarify.r \
	../panel.csv \
	../ancestry/output/sites.csv \
	../episcores/episcore-sites.csv \
	all_target_segments_not_covered_by_probes_withshifting.bed \
	output/stats.rmd
```

Output is here: [output/stats.md](output/stats.md)

