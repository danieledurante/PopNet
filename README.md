# PopNet: Nonparametric Bayes modeling of populations of networks

This repository is associated with the article: [Durante, Dunson and Vogelstein (2017) *Nonparametric Bayes Modeling of Populations of Networks*](https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1219260?journalCode=uasa20) (see also discussions and rejoinder [here](https://www.tandfonline.com/toc/uasa20/112/520?nav=tocList)), and aims at providing **materials and codes** to perform posterior inference via the **Gibbs sampler** described in the article (see file [`code`]() in the folder  [`Implementation`]()), with a focus on the **simulation study** in Section 5 (see data [`simulated_data.RData`]() in the folder [`Implementation`]()).

The analyses are performed with a **MacBook Pro (OS X El Capitan, version 10.11.6)**, using a `R` version **3.6.3**.

**NOTE**: All the materials are meant to provide a **basic `R implementation**, mostly meant to offer a clear understanding of the computational routines and steps associated with the proposed model. Optimized computational routines relying on C++ coding can be easily considered and may be required in the focus is on much larger networks.
