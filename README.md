# vALK
#### _A bioinformatic pipeline optimized for the processing and assessment of variants at the ALK gene locus_

>This project is the result of a [master's thesis](https://github.com/ivadym/masters_thesis) entitled: 
>
>_Development of a Bioinformatic Pipeline for the Dectection of Somatic Mutations in Liquid Bipsy Samples from Non-Small Cell Lung Cancer Patients_

### Pipeline characterization

To detect variants at the ALK gene locus, specific conditions for each of them were established from a set of previously validated samples. The following flowchart shows the structure of the developed pipeline and the selection criteria based on certain variables such as:
- **MOL_COUNT**: molecular coverage
- **READ_COUNT**: read coverage
- **MAPD**: median of the absolute values of all pair-wise differences
- **FILTER**: Ion Reporter™ internal filter (`Oncomine™ Variants v5.12`)
- **CLN_SIG**: clinical significance
- **VAR_CLASS**: Oncomine™ Variant Class
- **POS**: variant position
- **AF**: allele frequency

<img src="https://github.com/ivadym/vALK/blob/main/img/flowchart.png" width="600">

### Graphical User Interface (GUI)

The implemented user interface was developed to facilitate data input and interpretation of results.

<img src="https://github.com/ivadym/vALK/blob/main/img/GUI.png" width="600">

On the one hand, it is possible to select the filter/s to apply to the non-filtered-oncomine.tsv file, as well as the gene to study. Currently only the ALK gene has been addressed. On the other hand, the user can select both the source and output files through a pop-up screen.

Finally, to fully characterize the variants that have passed a particular filter, the results are displayed specifying the mutation type, the row number of the variants in the non-filtered-oncomine.tsv file, and the identifier of the alternate allele. Simultaneously, each of the identified variants is appended to the `.csv` output file along with its main characteristics.

### Performance

The initial sequenced samples were used to study the performance of the developed algorithm regarding the ALK mutations confirmed by dPCR, which was considered the gold standard. To assess this, sensitivity, specificity, and positive and negative predictive values, shown in the following table, were calculated based on 30 patients.

| **Statistic**             | **Value**     | **95% CI**        |
| -------------             | ------------- | -------------     |
| Sensitivity               |  87.50%       | 47.35% to 99.68%  | 
| Specificity               | 81.82%        | 59.72% to 94.81%  |
| Disease prevalence        | 20.00%        | -                 |
| Positive Predictive Value | 54.61%        | 32.31% to 75.20%  |
| Negative Predictive Value | 96.32%        | 80.55% to 99.40%  |
| Accuracy                  |  82.95%       | 64.83% to 94.13%  |

In this context, the algorithm managed to successfully identify 7 of the 8 mutations confirmed by dPCR, reaching a specificity of 87.50%. On the other hand, 4 patients were incorrectly classified as carriers of ALK mutations, thus a specificity of 81.82% was obtained. Finally, the algorithm reported an accuracy of 82.95%, with a remarkable performance in terms of discarding samples without any mutation since a negative predictive value of 96.32% was achieved. Noteworthy, the Oncomine™ Variants v5.12 filter only detected ALK mutations in 3 patients.

### Technical requirements
The programming language used in this study was `R v3.6.3`, which capabilities were extended through additional packages such as `Tcl/Tk` and `Scales v1.1.1`.
