# keller-hotspots
analyzing hotspots from keller et al 2018 data


The goal for this project is to analyze, with mediation analysis and pairwise pleiotropy tests, four of the five expression trait hotspots from the Keller et al. 2018 GENETICS data set. Pairwise pleiotropy tests will all use pairs that have one local expression trait and one nonlocal trait.

Note that we've already done [the analysis for Keller's Chr2 hotspot](https://github.com/fboehm/keller2018-chr2-hotspot-chtc
), so, in this repo, I consider only the four expression trait hotspots on Chromosomes 5, 7, 11, and 13. Note that Keller defined hotspots as having at least 100 associated transcripts.

Here is Keller's Table S1 (from the Supplemental Documents).

Chr   position (Mb)   number of transcripts   candidate mediator   number of genes with LOD difference > 1.5
---   -----------     ---------------------   ------------------   -----------------------------------------
2     165.5           147                     Hnf4a                88     
5     146             182                     Pdx1                 77
7     46              123                     Fam83e               96
11    71              126                     Sat2                 115
13    112.5           104                     Il6st                82

Table: From Keller et al. 2018 GENETICS, Supplementary table 1

