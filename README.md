# keller-hotspots
analyzing hotspots from keller et al 2018 data


The goal for this project is to analyze, with mediation analysis and pairwise pleiotropy tests, the five expression trait hotspots from the Keller et al. 2018 GENETICS data set. Pairwise pleiotropy tests will all use pairs that have one local expression trait and one nonlocal trait.

Note that we've already done [the analysis for Keller's Chr2 hotspot](https://github.com/fboehm/keller2018-chr2-hotspot-chtc
). However, the early analysis of Chr 2 was incomplete. Thus, in this repo, I consider all five expression trait hotspots on Chromosomes 2, 5, 7, 11, and 13. Note that Keller defined hotspots as having at least 100 associated transcripts.

Here is Keller's Table S1 (from the Supplemental Documents).

| Chr | position (Mb) |  number of transcripts  | candidate mediator  | number of genes with LOD difference > 1.5  |     
| --- | -----------   |  ---------------------  |  ------------------ |  ----------------------------------------- |   
| 2    | 165.5        |   147                   |  Hnf4a              |  88                                        | 
| 5    | 146          |   182                   |  Pdx1               |  77                                        |
| 7    | 46           |   123                   |  Fam83e             |  96                                        |
| 11   | 71           |   126                   |  Sat2               |  115                                       |
| 13   | 112.5        |   104                   |  Il6st              |  82                                        |

Table: From [Keller et al. 2018 GENETICS](https://www.genetics.org/content/209/1/335), [Supplementary table 1](https://figshare.com/articles/Supplemental_Material_for_Attie_et_al_2018_in_review_/5977459)

## Data

[Here](https://datadryad.org/resource/doi:10.5061/dryad.pj105) is the Data Dryad link for the data.


## Pairwise pleiotropy tests

For each of the four hotspots, we need to identify strong local expression QTL that will be counted as the set of local traits for that hotspot. We then pair each nonlocal trait with each of the local traits and test for pleiotropy. 

We also do all pairwise tests for pairs of nonlocal traits. 


## Mediation analysis

We assess whether each of the local traits is a mediator for the QTL - nonlocal trait association. We determine both LOD difference and LOD difference proportion statistics.




## Questions

What does it mean if multiple local traits show evidence of mediating the association between QTL and a single nonlocal trait? Is it enough to say that the local trait with the greatest LOD difference is the true mediator? What assumptions go into such a conclusion?









