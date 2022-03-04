
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/emosca-cnr/margheRita">
    <img src="vignettes/images/logo.png" alt="Logo" width="500" height="300">
  </a>

  <h3 align="center"> margheRita </h3>

  <p align="center">
    an R package for analyzing the entire workflow of mass spectrometry-based metabolic profiles !
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#Getting Started">Getting Started</a>
      <ul>
        <li><a href="#Prerequisites">Prerequisites</a></li>
        <li><a href="#R installation">R Installation</a></li>
        <li><a href="#installation">Git Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

margheRita is intended for any kind of mass spectrometer raw data, including both $MS$ and $MS/MS$ data. It takes as input `.txt/.csv` files that contain the metabolic profiles generated by **[MS-Dial](http://prime.psc.riken.jp/compms/msdial/main.html)** and metadata for sample processing. The first pre-processing step performs the integration of the metabolic profiles in a unique data structure and generates plots for quality control about outliers, drifts or batch effects. Then, metabolites and samples can be filtered based on the quantification of mass defect values, missing values and coefficient of variation of m/z features along QC and study samples. Moreover, margheRita provides various methods for missing value imputation and data normalization, including those particular recommended for metabolomic profiles, like the Probabilistic Quotient Normalization (PQN) or the normalization to a standard factor (e.g. protein concentration, urine osmolality). Subsequent analyses include: parametric and non-parametric statistical tests for the identification of significant metabolite signatures; sample clustering; metabolite correlation analysis and pathway analysis (using databases like **[KEGG](https://www.genome.jp/kegg/)** and **[Biocyc](https://biocyc.org)**). Additional advantage derives from the $MS/MS$ data management and annotations with the possibility of data browsing and searching.


<!-- GETTING STARTED -->
## Getting Started


### Prerequisites

To run margheRita you need [R](https://www.r-project.org/) version $\ge$ 3.8 and [Bioconductor](https://www.bioconductor.org/install/) 3.13 

### R Installation

1. You have first to install [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

    ```r
    if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")}
    ```
2. Install the repository and its dependencies from [Github](https://github.com/emosca-cnr/margheRita)

   ```r
    devtools::install_github("emosca-cnr/margheRita", dependencies = T)
    install_dependencies("margheRita")
   ```
3. load the package

   ```r
    load_all("margheRita")
   ```

### Git Installation
1. From bash

    ```sh
    cd path/to/R/library
    git clone https://github.com/emosca-cnr/margheRita.git
    ```
2. In R with [devtools](https://cran.r-project.org/web/packages/devtools/index.html) installed
   ```r
    load_all("margheRita")
   ```



<!-- USAGE EXAMPLES -->
## Usage

Here is presented a possible main workflow:   
```r
# 1) Read input files and creation of margheRita object
mRlist <- read_input_file(input = "inst/extdata/example1.xlsx", 
  metadata = "inst/extdata/example1_meta.xlsx", 
  split_QC = TRUE, rt_col = 2, mz_col =3, 
  data_start_col = 4)

mRlist

# 2) Filtering m/z features and NA values 
mRlist <- m_z_filtering(mRlist = mRlist, 
  lower_quality_mass_acc = 0.4, 
  upper_quality_mass_acc = 0.8, 
  do_plot = F, 
  color="black")

mRlist <- filter_NA(mRlist)

# 3) Normalization step with probabilistic quotient normalization
mRlist <- normalize_profiles(mRlist = mRlist, method = "pqn")

# 4) Principal component analysis
mRlist <- pca_gen(mRlist, dirout,
  col_by="class",
  scaling=c("none", "Pareto", "uv"),
  include_QC=TRUE,
  type=c("component", "coordinate"),
  dist.method="euclidean",
  top=Inf)

# 5) Calculation of mean median average and standard deviation of all the samples and FC for all comparisons 
mRlist <- mean_media_stdev_samples(mRlist, write_output = F)

mRlist <- calculate_lfc_all(mRlist = mRlist, lfc_theshold = 0.25)

# Annotation of compound
mRlist <- metabolite_annotation(feature_data = mRlist, reference = NULL,
  feature_spectra = NULL, reference_spectra= NULL, rt_err_thr=1,
  unaccept_flag=15, accept_flag=5, suffer_flag=10, acceptable_RI = 10,
  n_peaks=1, acceptable_PPM_err = 10)

```                                                                                               
_For more examples, please refer to the [Vignette](https://github.com/emosca-cnr/margheRita.git)_






<!-- LICENSE -->
## License

Distributed under the GPL License.



<!-- CONTACT -->
## Contact

Edoardo Bellini - [@Edoardo Bellini](https://www.researchgate.net/profile/Edordo-Bellini) - bellini.edoardo@hsr.it

Valeria Mannella - [@Valeria Mannella](https://www.researchgate.net/profile/Valeria-Mannella) - mannella.valeria@hsr.it

Maria Ulaszewska - [@Maria Ulaszewska](https://www.researchgate.net/profile/Maria-Ulaszewska) - ulaszewska.maria@hsr.it

Ettore Mosca - [@Ettore Mosca](https://www.researchgate.net/profile/Ettore-Mosca) - ettoremosca@itb.cnr.it

Marco J. Morelli - [@Marco J. Morelli](https://www.researchgate.net/profile/Marco-Morelli) - morelli.marco@hsr.it

Annapaola Andolfo - [@Annapaola Andolfo](https://www.researchgate.net/profile/Annapaola-Andolfo) - andolfo.annapaola@hsr.it

Denise Drago - [@Annapaola Andolfo](https://www.researchgate.net/profile/Annapaola-Andolfo) - drago.denise@hsr.it

Project Link: [https://github.com/your_username/repo_name](https://github.com/your_username/repo_name)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[product-screenshot]: images/screenshot.png
