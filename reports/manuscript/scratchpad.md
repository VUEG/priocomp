### Introduction

+ Why are ESs relevant for BD conservation in general?
  + Trade-offs important
  + More studies are explicitly incorporating ESs into the analyses
+ What is spatial conservation prioritization (SCP)?
  + SCP-methods have been primarily developed for biodiversity
+ Multiple SCP methods are available
  + This is also reflected in what types of datasets are used
  + More data covering broad geographical ranges and species is becoming available (RS)
+ How is SP for ESs different?
  + What has already been done?
  + Why is application of SCP methods nevertheless useful?
+ Why can direct application of SCP methods to ESs create issues? (**MISSING**)
+ In this manuscript, we...




### Discussion

##### Rank priority patterns

- While our objective here was not primarily to do a directly	policy-relevant prioritization, the results already as they are informative.
- The general priority rank patterns are generally explained and in line with...
  - Kukkala et al. (2016)
  - Schulp et al. (2014)
  - Mouchet et al. (2017)
- Priority rank patterns
  - ES: These regions had high values for erosion prevention in particular. Carbon sequestration and wood production had a minor, but still noticeable effect concentrating priorities in Southern Fennoscandia and South-Western France. The effects of other ES features were more diffuse. 		
  - BD: biomodal distribution
- Stats
  - Not surprisingly, the highest priorities overlapped more than the lowest
- Performance and trade-offs
  - Cimon-Morin et al. (2016): “ Achieving all conservation targets starting from a network that was primarily built for either ES or biodiversity features alone was two to five times less efficient than considering both ES and biodiversity simultaneously in conservation assessment. A better framework is required to translate these spatial synergies into effective joint conservation actions.”
- **Bottom line**: the data and purpose of the prioritization matters

##### Method performance

+ We show that performance-wise, all methods are very similar
  + Only for ES, Zonation seemed to produce a more balanced solution. This is probably because more spatially aggregated distributional patterns in the BD features, and the higher number of features. In any case, the result is related to the iterative cell removal process
  + Previous studies (e.g. Albuquerque & Beier, 2015) have suggested similar performance for RWR and ZON (though with CAZ) and suggested that they are close to optimal. Here, we show similar results.
+ The practical aspect of the prioritization analysis are different too.
  + RWR: 
    + Easy to set up
    + Runtimes:
      + ALL: 600 s (~10 min), ES: 12 s, BD: 573 s (~9.5 min)
  + ZON:
    + Building analysis variants takes expertise
    + Excellent reporting functionalities
    + Runtimes:
      + ALL: 220024 s (~2.55 days), ES: 1679 s (~30 min), BD: 206420 s (2.39 days)
  + ILP
    + Data preprocessing almost the same to RWR
    + Solver if proprietary but available with an academic license
    + Runtimes:
      + ALL: 3640 s (~60 min), ES: 2162 s (~36 min), BD: 3476 s (~58 min) 
    + Major advances lately as demonstrated by Beyer et al. (2016)
  + NOTE! All the above are heavily dependant on the implementation details. 
+ However, it is important to note that the prioritization analyses we have used are simple; For any real-life application several complicating factors would probably need to be included. 




##### Integrated prioritization of biodiversity conservation andecosystem services supply

+ Trade-offs between ecosystem service provision and biodiversity conservation are most likely common, ***as we have shown***. It does not necessarily follow that priority areas for the provision of ecosystem services are automatically priority areas also for biodiversity (Anderson et al., 2009; Thomas et al., 2012). 

+ It still remains unclear how exactly ecosystem services should best be incorporated into prioritization schemes that have been developed with biodiversity conservation in mind.

  + Complicating factors
    + Spatiotemporal scales
    + Supply/demand
    + Places/actions
    + Connectivity (Kukkala & Moilanen 2016)
  + Results would probably change a lot if ES flows and demands would be accounted for through e.g. connectivity (Kukkala & Moilanen, 2016).

+ Many prioritization methods designed with biodiversity features in mind give high emphasis on relative rarity of features. Furthermore, ways of accounting for spatial arrangement of features are predominantly done with the ecology of species mind (e.g. ecological connectivity).

+ Furthermore, asymptotic benefit-functions often used in SCP (Arponen et al., 2005; Wilson et al., 2009)⁠ are not suitable for all ecosystem services for which either linear or more complex relationships are more appropriate (Barbier et al., 2008; Luck et al., 2012)⁠. Emphasizing the relative importance of rare features over more common features (Moilanen et al., 2005; Williams et al., 1996)⁠ is another principle which may be more suitable for biodiversity rather than ecosystem services features.  

+ Ecosystem services require additional considerations (e.g. the availability of alternative means of providing benefits supplied by services, human demand, and the scale of, and site dependency in, the delivery of services), which can be considerably different to those of species. In terms of selecting the suitable prioritization tool for the job, we recommend the following:

  1. Study the assumptions behind the tools you are about to use. Is it geared more towards species or ES?
  2. Embrace flexibility, but avoid complexity.

  ​




##### Selecting the right tool

+ Often, the tool is selected first and other considerations follow
+ Operationalizing ecosystem services requires institutional adaptation, case-specific tailoring of methods, and deliberation among practitioners and stakeholders (Rinne and Primmer, 2015)
+ Qualitative method comparison. Bagstad et al. (2013) provide a useful set of criteria against which methods and software tools can be assessed.
  1. The fit of the underlying assumption/concepts with topic and research question
  2. Quantification and uncertainty
  3. Time requirements
  4. Capacity for independent application
  5. Level of development and documentation
  6. Scalability
  7. Generalizability
  8. Nonmonetary and cultural perspectives
  9. Affordability, insights, integration with existing environmental assessment
+ A further practical advantage of adapting SCP methods into spatial prioritization of ES capacity is that SCP has already seen wide operationalization in real-life decision-making (Kareksela et al., 2013; Knight et al., 2009; Lehtomäki and Moilanen, 2013; Pressey et al., 2013; Whitehead et al., 2016)⁠. Based on experiences from a broad array of applied projects and the existing literature on the applicability of methods and tools to practice, it is possible to assess the potential of SCP methods in the context of ecosystem services. In the broader context of providing decision-support tools capable of dealing with ecosystem services, there are multiple good reviews available assessing the technical and practical aspects of different software tools (Bagstad et al., 2013; Langemeyer et al., 2016)⁠. However, only few assessments explicitly consider spatial methods and combining multiple ecosystem provision and biodiversity conservation features simultaneously. 
+ In addition to quantifying the differences between the different spatial prioritization methods, we also present the full implementation of the analysis that can be adjusted to other types of data. 