
HIPPIE (Human Integrated Protein-Protein Interaction rEference) PPI data
downloaded from http://cbdm.uni-mainz.de/hippie/ on 12 June 2021
with 
  wget http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt
  wget http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-current.mitab.txt
for resprectively HIPPIE tsv format and PSI-MITAB 2.5 format
(http://psicquic.github.io/MITAB25Format.html).

As per download page http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
current version is v2.2 last updated 02/14/19

See also HIPPIE_download_page.pdf (saved from 
http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php) and
HIPPIE_informatino_page.pge (saved from
http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/information.php).

List of scores assigned to experimental techniques (used to score the 
interactions in HIPIE) downloaded from
http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/RS/experimental_scores.tsv

Citations for HIPPIE:


    Schaefer MH, Fontaine J-F, Vinayagam A, Porras P, Wanker EE, et al. (2012) HIPPIE: Integrating Protein Interaction Networks with Experiment Based Quality Scores. PLoS ONE, 7(2): e31826 [link]

    Schaefer MH, Lopes TJS, Mah N, Shoemaker JE, Matsuoka Y, et al. (2013) Adding Protein Context to the Human Protein-Protein Interaction Network to Reveal Meaningful Interactions. PLoS Computational Biology, 9(1): e1002860 [link]

    Suratanee A, Schaefer MH, Betts M, Soons Z, Mannsperger H, et al. (2014) Characterizing Protein Interactions Employing a Genome-Wide siRNA Cellular Phenotyping Screen. PLoS Computational Biology, 10.9: e1003814 [link]

    Alanis-Lobato G, Andrade-Navarro MA, & Schaefer MH. (2016). HIPPIE v2.0: enhancing meaningfulness and reliability of protein-protein interaction networks. Nucleic Acids Research, gkw985 [link]


PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification
files downloaded from

http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR16.0_human

(human, Panther version 16.0) on 21 June 2021.

See: http://pantherdb.org/about.jsp

and also the downloaded README file downloaded from 
http://data.pantherdb.org/ftp/sequence_classifications/current_release/README
as README_panther.txt

Citations for PANTHER are:

Thomas, P. D., Campbell, M. J., Kejariwal, A., Mi, H., Karlak, B., Daverman, R., ... & Narechania, A. (2003). PANTHER: a library of protein families and subfamilies indexed by function. Genome research, 13(9), 2129-2141.

Mi, H., Muruganujan, A., & Thomas, P. D. (2013). PANTHER in 2013: modeling the evolution of gene function, and other gene attributes, in the context of phylogenetic trees. Nucleic acids research, 41(Database issue), D377–D386. https://doi.org/10.1093/nar/gks1118

Mi, H., Muruganujan, A., Huang, X., Ebert, D., Mills, C., Guo, X., & Thomas, P. D. (2019). Protocol Update for large-scale genome and gene function analysis with the PANTHER classification system (v. 14.0). Nature protocols, 14(3), 703-721.

Mi, H., Lazareva-Ulitsky, B., Loo, R., Kejariwal, A., Vandergriff, J., Rabkin, S., ... & Thomas, P. D. (2005). The PANTHER database of protein families, subfamilies, functions and pathways. Nucleic acids research, 33(suppl_1), D284-D288.

Mi, H., Ebert, D., Muruganujan, A., Mills, C., Albou, L. P., Mushayamaha, T., & Thomas, P. D. (2021). PANTHER version 16: a revised family classification, tree-based classification tool, enhancer regions and extensive API. Nucleic Acids Research, 49(D1), D394-D403.

"How to cite" says to always cite Mi et al. (2021); and for classification
system, Mi et al. (2019).

Citations for Gene Ontology (GO) are:

Ashburner, Michael, Catherine A. Ball, Judith A. Blake, David Botstein, Heather Butler, J. Michael Cherry, Allan P. Davis et al. "Gene ontology: tool for the unification of biology." Nature genetics 25, no. 1 (2000): 25-29.

The Gene Ontology Consortium, "The Gene Ontology resource: enriching a GOld mine." Nucleic Acids Research 49, no. D1 (2021): D325-D334.


The extraction, conversion and matching of the HIPPIE data and the PANTHER
data are described in:

Stivala, A., Lomi, A. Testing biological network motif significance with exponential random graph models. Appl Netw Sci 6, 91 (2021). https://doi.org/10.1007/s41109-021-00434-y

and the scripts and data are available from:

https://github.com/stivalaa/bionetworks_estimations

For more details see the 
scripts/getHippiePPIhighConfidenceEdgelistWithPantherAttributes.R
R script in that repository.

The hippie_ppi_high_edgelist.txt is copied directly from 
undirected/hippie/model2/hippie_ppi_high_edgelist.txt
file in that repository, and the ategorial attribute for cellular component 
extracted from the file
undirected/hippie/model2/hippie_ppi_high_actors.txt
with:

cat ~/bionetworks_estimations/undirected/hippie/model2/hippie_ppi_high_actors.txt | sed -n '/^Categorical Attributes:/,$p' | tail -n+2 | awk '{print $2}' > hippie_catattr.txt

ADS
Wed Feb  2 03:36:12 CET 2022
