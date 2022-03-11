
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
        <li><a href="#Installation">Installation</a></li>
        <li><a href="#Installation">Documentation</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

margheRita is intended for any kind of mass spectrometer raw data, including both MS and MS/MS data. It takes as input `.xls/.txt/.csv` files that contain the metabolic profiles generated by **[MS-Dial](http://prime.psc.riken.jp/compms/msdial/main.html)** and metadata for sample processing. The first pre-processing step performs the integration of the metabolic profiles in a unique data structure and generates plots for quality control about outliers, drifts or batch effects. Then, metabolites and samples can be filtered based on the quantification of mass defect values, missing values and coefficient of variation of m/z features along QC and study samples. Moreover, margheRita provides various methods for missing value imputation and data normalization, including those particular recommended for metabolomic profiles, like the Probabilistic Quotient Normalization (PQN) or the normalization to a standard factor (e.g. protein concentration, urine osmolality). Subsequent analyses include: parametric and non-parametric statistical tests for the identification of significant metabolite signatures; sample clustering; metabolite correlation analysis and pathway analysis (using databases like **[KEGG](https://www.genome.jp/kegg/)** and **[Biocyc](https://biocyc.org)**). Additional advantage derives from the MS/MS data management and annotations with the possibility of data browsing and searching.


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

To run margheRita you need [R](https://www.r-project.org/) version $\ge$ 3.8

### Installation

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
3. Load the package

  ```r
    load_all("margheRita")
    i Loading margheRita

         ##############################################
         #                                            #
         #          Welcome to margheRita             #
         #                 from                       #
         #              all our team                  #
         #     we hope you to enjoy our tool          #
         #                                            #
         ##############################################
  ```

### Documentation

For more examples, please refer to the [Vignette](https://github.com/emosca-cnr/margheRita.git)_


<!-- LICENSE -->
## License

Distributed under the GPL License.


<!-- CONTACT -->
## Contacts

Ettore Mosca - CNR-ITB, Bioinformatics Lab  [@Ettore Mosca](https://www.researchgate.net/profile/Ettore-Mosca) - ettore.mosca@itb.cnr.it

Annapaola Andolfo - HSR, Proteomics And Metabolomcis Facility [@Annapaola Andolfo](https://www.researchgate.net/profile/Annapaola-Andolfo) - andolfo.annapaola@hsr.it

Project Link: [https://github.com/emosca-cnr/margheRita](https://github.com/emosca-cnr/margheRita)


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
