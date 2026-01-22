# Genomics underlying facultative CAM variability in a homoploid hybrid, *Yucca gloriosa*

Project abstract:
Photosynthesis is central to plant carbon acquisition, but the typical C3 pathway is inhibited by abiotic stresses such as water limitation due to the dual function of stomata in carbon intake and transpiration. Crassulacean acid metabolism (CAM) is a form of photosynthesis in which stomata are opened mostly at night, when temperature is reduced, thus reducing transpiration and increasing the plant’s water use efficiency. CAM exists across a spectrum in many lineages, including in the genus Yucca; some species, such as Y. filamentosa, are strictly C3, whereas others, such as Y. aloifolia, use mainly the CAM pathway, and others use CAM only under certain conditions (“facultative CAM” or “C3+CAM”). One such species is Y. gloriosa, a homoploid hybrid of CAM Y. aloifolia and C3 Y. filamentosa, which exhibits intraspecific variation in the degree of facultative CAM under drought.
	The combination of homoploidy and unique CAM biology makes Y. gloriosa an excellent system for the study of both CAM evolution and homoploid speciation. To this end, a 24-hour timecourse RNA-seq dataset was collected for a controlled drought experiment run on both parent species and 26 naturally occurring genotypes of Y. gloriosa. These data are being leveraged to understand both the expression variability arising from homoploid hybridization and that underlying intraspecific CAM variation in this species.

#### Subdirectories:
- `camgene_modeling`: Contains scripts related to polynomial modeling of CAM gene expression over time.
- `data_wrangling`: Contains various notebooks/scripts for general data wrangling.
- `differential_expression`: Contains scripts for DESeq2 (both general & for homeolog expression bias analysis).
- `machine_learning`: Contains scripts & notebooks relevant to XGBoost machine learning modeling.
- `masigpro`: Contains scripts for running maSigPro (time-ordered gene coexpression).
- `pca`: Contains notebooks used to create expression principal component analysis (PCA) plots.
- `physiology`: Contains notebooks/code relevant to making calculations on physiology data (see [Karolina Heyduk's GitHub repo](https://github.com/kheyduk/Yucca_physiology))
- `rnaseq_processing`: Contains scripts relevant to RNA-seq data processing (i.e. creating TPM & counts matrices).
