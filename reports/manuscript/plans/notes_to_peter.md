## Version 0.8.0 (2017-04-05)



### Main text



#### General

+ After changing the scope slightly, there is probably a lot redundant text. The length can be cut down.

+ Introduction and Discussion are still a bit of a mess, I would appreciate comments/focus especially on these.

+ I'm still working on putting the Supplementary Material together, but you can start with the main MS.

  ​



#### Introduction

+ Fig. 1 is still unfinished, tell me what you think about it (some reference already in Introduction and Discussion).

  ​



#### Material

+ Quite long already. Including the cost description here would't necessarily make the section that much longer.
+ ILP: I'm having serious doubts that there's something wrong with my implementation of ILP:
  1. ILP does not have the highest performance in any measure I could think of. IF the solutions produced by ILP would be truly optimal, then it should also have highest feature representation conditioned on the area target. This could be the same as RWR, though.
  2. RWR_ES / ILP_ES and RWR_BD / RWR BD are almost exactly the same, whereas RWR_ALL / ILP_ALL are not. In ILP_ALL (as in RWR_ALL and ZON_ALL), I give both feature groups (ES and BD) the same aggregate weight (1.0). I'm starting to suspect that this aggregate weighting in ILP is not working as intended.
  3. We also need a contingency plan: if RWR and ILP turn out to be practically the same, then the main conclusions about balance and performance don't hold and we're again back to square one. I can see at least the following solutions:
     1. Use costs in the main line of work. This introduces differences, but makes everything else less clear. 

     2. Drop RWR or ILP. I don't think contending "they are exactly the same" is very good as there's no reason to expect they wouldn't be.

     3. Come up with alternatives for either RWR or ILP (implementation). 

        ​



#### Results

+ This section is quite short. In fact, it's not completely impossible to include costs here (now in the supplements

  ​



#### Discussion

+ There's still quite some interpretation of the priority pattern within Europe. Some is needed for sure, but I don't want to concentrate on it too much.

  ​



#### Conclusions

+ To be added when the content is fixed

  ​

### Appendix

+ Note that there currently are both Appendix (at the end of the main text) and Electronic Supplementary Material (separate file). Appendix has miscellaneous content relevant for the main text, Supplement is costs only.

  ​

### Supplementary Material

#### Costs

+ All costs related analyses are now in the Supplement, which I'm still working on. However, 1) why we chose to do the cost analyses, 2) and why then do we place them in the supplement needs more work.

