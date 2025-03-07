# EndothermEcoEvo
R Code for the simulations of the paper:
Rubalcaba, J. G. (under review). "The evolution of homeothermic endothermy via life-history optimization"

Juan G. Rubalcaba
Department of Biodiversity, Ecology and Evolution; Faculty of Biological Sciences; Complutense University of Madrid
jg.rubalcaba@gmail.com

Endothermy is an energetically expensive trait and yet it posed an evolutionary advantage across different lineages – a paradox that remain puzzling biologists. Here I investigate whether endothermy can evolve through life-history optimization using a model of the balance between energy assimilation vs energy allocation to either somatic maintenance, thermoregulation, growth, or reproduction. The model displays bistable strategies when assimilation rates and thermoregulatory costs increase, respectively, exponentially and linearly with body temperature: the ‘heterothermic strategy’ consists of minimizing the costs of thermoregulation by maintaining body temperature close to ambient temperature; and the ‘homeothermic strategy’ consists of increasing body temperature until the costs of thermoregulation are fully compensated by the increased assimilation capacity at higher temperatures. These strategies produce similar fitness outcomes and thus emerge as alternative stable states of the system, maintained by strong stabilizing selection preventing transitions between them. Using quantitative genetics simulations, I show that a drop in ambient temperature may push populations towards an evolutionary branching point enabling the rapid radiation of homeothermic lineages coupled with body size reductions. I thus propose that life-history optimization of energy balance can explain the radiation of homeothermic endothermy associated with either climate cooling or migration to colder regions by early endothermic lineages.

Code developer: Juan G. Rubalcaba

# R version and required packages:
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

RColorBrewer_1.1-3 
ggpubr_0.6.0       
ggplot2_3.5.1    

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       stringi_1.8.4     purrr_1.0.2       car_3.1-2        
 [7] generics_0.1.3    glue_1.7.0        backports_1.5.0   colorspace_2.1-1  scales_1.3.0      fansi_1.0.6      
[13] abind_1.4-8       carData_3.0-5     munsell_0.5.1     tibble_3.2.1      rstatix_0.7.2     lifecycle_1.0.4  
[19] ggsignif_0.6.4    compiler_4.4.1    dplyr_1.1.4       pkgconfig_2.0.3   tidyr_1.3.1       rstudioapi_0.16.0
[25] R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3    tools_4.4.1      
[31] withr_3.0.0       gtable_0.3.5      broom_1.0.6  

# Code description
The main R code computes the energy-balance model described in the paper and generates Figures 2, 3 and 4. 
The additional code "functions_multiple_scales.R" loads functions used for Figure 3 and 4.









