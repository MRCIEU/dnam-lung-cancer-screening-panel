Quality control analyses and outputs
for the probes used in the final panel.

```
Rscript src/check.r \
	../panel.csv \
	data/all_target_segments_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed \
    data/all_target_segments_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed \
    data/all_target_segments_not_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed \
	data/merged_probe_file_shareable_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed \
	output/stats.rmd
```

Output is here: [output/stats.md](output/stats.md)

