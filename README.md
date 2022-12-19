# Consistent gene expression of Oesophageal squamous cell carcinoma, Glioblastoma and Melanoma.
## I took 7 cell lines which are related to Oesophageal squamous cell carcinoma, Glioblastoma and Melanoma. I am uploading my Master's project here in this repostory.
I have used Rstudio, Linux, Jupiter notebook, Visual studio here with different types of library in this project.
I am uploading my snakemakefile which i used in this project to download my SRA files, the fastqc files to check the quality control of my data, to build genome index using HISAT2, for the aloignment (single and pairwise), for the conversion of my reads into BAM files, feature counts and normalizing the feature counts.
I used Multiple array veiwer to visualize my data with the help of Rank Products. Moreover i have done also done the Functional analyses of my data using Rstudio.

Consistent gene expression of Oesophageal squamous cell carcinoma, Glioblastoma and Melanoma.

Table of Contents
Abstract:	2
Aim:	2
Methods and Materials:	3
1. Introduction	3
1.1 LnCaP (Lymph Node Carcinoma of the Prostate):	5
1.2 KYSE150, KYSE70, KYSE450:	6
1.3 U-CH11	7
1.4 D19	8
2. Methods and Materials.	9
2.1 Materials	9
2.2.1 Method:	9
2.2.2 Collecting the Data:	9
2.2.3 Extracting the data.	10
2.2.4 Quality control of the data.	11
2.2.5 Building the genome index and Mapping:	16
2.2.6 Alignment:	18
2.2.7 Conversion of reads:	20
2.2.8 Feature Counts:	22
2.2.9 Normalizing the FeatureCounts.	24
2.2.10 Multiple Array Viewer:	25
2.2.11 Functional analyses:	29
2.2.12 Packages	29
3. Results:	30
3.1 MultiQC:	30
3.2 UpsetR analysis:	33
    Extracted Gene ID with link:	37
3.4 Functional analyses of Consistent genes:	37
4.Discussion:	49
5.References:	51



Abstract:
Harmful development is a principal wellspring of death all over the planet. Despite numerous investigation movements in the field, the genetic changes dealing with the difference in common oral cells into perilous cells have not been totally explained. A couple of assessments have evaluated carcinogenesis at the sub-nuclear level. Harmful development cell lines are typically used in biomedical assessment since they give a boundless wellspring of cells and address various periods of beginning and development of carcinogenesis. The multistep cycle of oral carcinogenesis requires the amassing of a few hereditary changes that influence proto-oncogene and cancer silencer quality exercises. The examination of synthetic alterations known as epigenetic adjustments, which control the openness and clarity of DNA, is important to recognize a few changes that can't be found essentially by sequencing DNA. Neoplastic change results from the aggregation of epigenetic adjustments. The Stupp routine, which incorporates radiation and chemotherapy following a medical procedure, is as of now utilized as first line treatment for glioblastoma (GBM), the most widely recognized and forceful glial growth. The strong threat glioblastoma generally returns and constantly kills the patient notwithstanding forceful medical procedure, radiation, and simultaneous and adjuvant temozolomide therapy. Notwithstanding progress nearby, the guess for people with high-risk or high level metastatic melanoma stays grim. Melanoma occurrence is rising worldwide. Medical procedure is trailed by adjuvant treatment or clinical preliminary enrolment for patients with thick (2.0 mm) essential melanoma regardless of territorial lymph hub metastases. Cell lines are utilized by melanoma scientists to rein act an assortment of cancer peculiarities. In this manner, it is vital to grasp the similitudes and contrasts between cell lines and the cancers. Finding that melanoma cell lines transcriptional and mutational characteristics equal those of ordinary cell lines.
Aim
To find the reliable qualities between the given cell lines and to do the utilitarian investigations with the goal that we can track down a medication.
Methods and Materials:
Several databases, including the NCBI, Uniport, Science direct, PubMed, Google Scholar, Expasy, Array-express using the alternative arrangements of the essential words below like Homo-sapiens, cancer cell lines, U-CH11, LnCaP, Mel-888, KYSE450, KYSE70, KYSE150, D19, series, high throughput sequencing.

Keywords: U-CH11, LnCaP, Mel-888, KYSE450, KYSE70, KYSE150, D19, consistent genes, functional analyses, Cell lines, Python, RStudio, Humans, Homo-sapiens, UpsetR, TDM, pre-processor, Euler, Quantile, Packages, Linux, Data, Research, Biological, Cancer,
NCBI, Linux, Alignment, Fastqc, Package.

1. Introduction:
Cell lines are in vitro model frameworks that are broadly utilized in a few areas of clinical exploration, including fundamental disease examination and drug improvement. Their utilization stems for the most part from their ability to deliver a boundless stockpile of organic material for logical purposes. Under the right settings and with the right controls, validated disease cell lines protect most of the hereditary highlights of the first growth. Contrasting genomic information from most malignant growth cell lines throughout recent years has upheld this declaration and those found while researching tumoral tissue reciprocals in the Cancer Genome Atlas (TCGA) data set. We presently have unmatched admittance to point by point open access cell line information bases reporting their atomic and cell changes. For different malignant growth types, including glioma, bosom, colorectal, and ovarian disease, contrasts in genomic and transcriptional designs between disease cell lines and cancer tests have been analyzed. The most frequently utilized cell lines, for high-grade serous ovarian malignant growth, for example, showed that they are not exceptionally illustrative of the cancer partners.
Current disease research actually depends vigorously on human malignant growth cell lines. As a matter of fact, they are habitually used as preclinical model frameworks to comprehend how systems work and how to treat patients. Prominently, the improvement of - omics innovations has prompted the production of broad information bases for the portrayal of most of realized cell lines. Moreover, the data got from these examinations that was made accessible web-based filled in as an important asset for the investigation of malignant growth cell lines and made it simpler for scientists to pick the best in vitro model framework for their exploration projects. It is urgent to consider various imperative distributions that have been delivered in under 10 years in this specific situation.Transcriptomics has undergone a revolution thanks to the discovery of high-throughput next-generation sequencing (NGS), which makes it possible to analyse RNA through the sequencing of complementary DNA (cDNA). Our knowledge of the intricate and dynamic nature of the transcriptome has been completely transformed by this technique, known as RNA sequencing (RNA-Seq). Gene expression, alternative splicing, and allele-specific expression are all better understood and quantified by RNA-Seq. Deep profiling of the transcriptome and the possibility to shed light on many physiological and pathological situations have been made possible by recent improvements in the RNA-Seq methodology, which span sample preparation, sequencing platforms, and bioinformatic data interpretation.
Carcinoma is the name given to a gathering of tumours that beginning in epithelial cells. These cells make up the tissue that lines the surfaces inside and outside your body. This tissue, called the epithelium, is tracked down outwardly surfaces of your skin and inner organs. It additionally incorporates the internal parts of empty organs, similar to your gastrointestinal system and veins.
Carcinoma is the most ordinarily analysed sort of disease. It's characterized into subtypes in light of the area and kind of the epithelial cell it begins from. Squamous cells, found in the epithelium's outermost layer, are where squamous cell carcinoma occurs. Squamous cell carcinoma most frequently relates to skin cancer, however it may as frequently afflict other sections of the body. A category of malignancies known as adenocarcinomas begins in glandular cells, which are specialised epithelial cells. Most organ linings contain glandular cells, which release something akin to mucus. Adenocarcinomas that are most prevalent include Colorectal cancer, Lung Cancer, Pancreatic Cancer, Prostate Cancer, Breast Cancer. The overall risk factors for carcinoma are the same as the risk factors for all other malignancies. These risk elements consist of: bad eating habits, with ageing, abuse of drugs like alcohol and tobacco, genetics exposure to certain chemicals. Exposure to UV radiation (including tanning beds) is a substantial risk factor for both basal cell carcinomas and squamous cell carcinomas, two types of skin cancer. Melanoma is a type of cancer that develops in pigment-producing cells called melanocytes. It is the most dangerous form of skin cancer and can lead to death if left untreated. One of the most effective ways to prevent melanoma is by performing regular self-examinations and identifying new symptoms early. Sadly, many people are not aware that they have this disease until it has already spread, making it difficult to reverse the damage.


1.1 LnCaP (Lymph Node Carcinoma of the Prostate):
High-fondness explicit androgen and estrogen receptors are tracked down in the cytosol and atomic portions. The LNCaP line is hormonally touchy, as shown by 5 alpha-in dihydrotestosterone's vitro guideline of cell development and corrosive phosphatase creation. LNCaP cells additionally express Prostate Specific Antigen (PSA). The most predominant malignant growth kind and a central point in men's disease related mortality is prostate disease. Nowadays, most of recently analyzed cases are locally forceful and nearly sluggish, requiring no dynamic treatment. The circumstance is more serious for guys with cutting edge prostate malignant growth. Androgen-hardship therapy (ADT) is the main line of safeguard since it assumes an essential part in all components of ordinary prostate disease capability as well as all periods of prostate malignant growth development.  
Figure 1: Prostate malignant growth improvement graph showing the many periods of the ailment and the job of LNCaP cells and their subsidiaries in grasping the different stages.
Maiming safe prostate disease (CRPC), so named in light of the fact that it keeps on depending on AR even without any androgens, in the end creates from these cancers when they return in most of cases. With the bone as the primary area, CRPC is ordinarily lethal and normally metastatic (mCRPC). Extra treatment, for example, the utilization of second-age antiandrogens, may bring about considerably more serious illness types that are much of the time AR-negative and may foster qualities of neuroendocrine prostate disease (NEPC). These forceful varieties have altogether less fortunate infection results than mCRPC because of their high bone and delicate tissue metastasis rates.Some prostate tumours are considered to progress more rapidly as a result of mutations that alter the selectivity of steroid-binding. In actuality, the LNCaP cell line included the first mutation in the androgen receptor discovered in a prostate cancer. This cell line is an androgen-responsive cell line derived from a human prostate cancer that was metastatic to bone. The functional consequence of this mutation in exon 8 of the androgen receptor is to change the steroid-binding specificity of the androgen receptor, resulting in oestradiol and the anti-androgens hydroxy-flutamide, nilutamide, and cyproterone acetate becoming agonists. Since the discovery of the first androgen receptor mutation in prostate cancer, mutations in the androgen receptor have been found in additional prostate cancer specimens and cell lines. These mutations alter or expand the androgen receptor's binding specificity for glucocorticoids, some adrenal androgens, weaker androgens, anti-androgens, estrogens, and progesterone. This sort of mutation has clear consequences in human prostate cancer. An anti-androgen or an oestrogen administered to the patient would be harmful. The antagonist has changed from being a tumour growth inhibitor to an agonist for the altered receptor, which may promote tumour development. All of these mutations would be categorised as gain of function mutations.
1.2 KYSE150, KYSE70, KYSE450:
One of the six most pervasive dangerous growths on the planet is esophageal disease (EC). A cell subpopulation in the growth tissue might be vital to carcinogenesis, neoplasia, and metastasis, as per mounting information. Malignant growth undifferentiated organisms (CSCs) are a gathering of cells that are profoundly metastatic and impervious to chemotherapy and radiation treatment. Besides, research has shown that CSCs are impervious to both radiation and chemotherapy. In spite of the way that researchers have removed CSCs from an assortment of cancer tissues and cell lines.
Esophageal SSC (ESSC) is more normal in non-industrial countries, particularly China. From Northeast China to the Middle East, there is a noticeable belt of esophageal disease occurrence, especially SCC. Smoking, stoutness, utilization of hot beverages and red meat, a high liquor consumption, and an unfortunate admission of new leafy foods are undeniably connected to the high commonness of ESCC . Furthermore, there are no dependable screening strategies and early esophageal disease has no conspicuous clinical signs, making early discovery troublesome in clinical practice. Since most of patients have progressed sickness when they get a conclusion and the disease has spread, this is valid.
1.3 U-CH111.3 U-CH11
sicknesses. As the drawn out articulation of the undeveloped record factor brachyury is the most intermittent trademark, chordoma cells are accepted to have created from the notochord's cell remains. Elevated degrees of atomic brachyury, communicated by the TBXT quality, appear to essentially add to carcinogenesis as well as being a delicate chordoma marker. Subsequently, it has been shown that TBXT concealment brought about by shRNA and sgRNA in chordoma cell lines represses cell expansion and triggers apoptosis.An uncommon form of bone cancer called chordoma often affects the spine or the skull. It often develops near the base of the spine or where the skull rests on top of the spine (sacrum)Chordoma starts in cells that were originally a group of cells in the growing embryo that would later develop into the discs of the spine. By the time you are born or shortly after, the majority of these cells disappear. However, a small number of these cells can occasionally survive and, in rare cases, become cancerous. Although it can occur at any age, chordoma most frequently affects individuals between the ages of 40 and 60.
1.4 D19
Females incorporate D-19. In clinical examinations, tamoxifen and radiation therapy are utilized to treat glioblastoma multiforme (GBM). The way that tamoxifen is a radiosensitizer is the legitimization for this treatment. Be that as it may, there isn't a lot of help for this. Consequently, the creators examined the effect of blend radiation and tamoxifen therapy in three GBM cell lines got from people. 
Figure 2: Cycle of hypoxia and neovascularization in glioblastoma.
1. Glioma cells use the oxygen that the healthy vasculature provides. 2. Vaso-occlusion and necrosis are brought on by endothelial damage, prothrombotic factors, and increased mechanical pressure in areas with a high density of glioma cells. 3. In response to hypoxia, perivascular glioma cells adopt a "go" phenotype. 4. Pro-angiogenic substances are secreted by pseudo palisading glioma cells. 5. The development of aberrant, highly permeable neo vasculature is induced by pro-angiogenic agents, leading to increased hypoxia and rapid progression. 6. The cycle restarts when pseudo palisading cells move to a new vasculature.

2. Methods and Materials.
2.1 Materials 
There were 7 cell lines. The raw data was taken from various databases such as NCBI, Array-express. Furthermore, the databases like Cellosaurus, Uniport, PubMed, European nucleotide archive. There were various packages and libraries were used in Ubuntu, RStudio and Python. Packages such as UpSetR, preprocessCore, 
Libraries like UpSetR, TDM, preprocessCore, euler, Quantile

2.2.1 Method:
We will utilize information from different data sets like NCBI GEO and Array-express. We will utilize the provided promotion number in the data set to look through the Gene Expression Omnibus (GEO) to find the crude sequencing information. The GEO page for promotion number might be found at https://www.ncbi.nlm.nih.gov/gds/?term= . There is data on this page concerning the areas of the stockpiling for the sequencing information for each exploration test.
.2.2 Collecting the Data:
Succession Read Archive (SRA) records are the organization utilized by the NCBI to store sequencing information. Every one of the examples is connected to a bunch of the previously mentioned SRA promotion numbers. Downloading the SRA runs for each example is the initial step. The SRA documents will then be utilized to make FASTQ records. The information design expected for mass RNA-sequencing examination is the FASTQ record. We will get the crude readings of every informational index from the SRA pursue choice assembling the information from the data sets. Different strategies exist for downloading the information freely of NCBI.
The European Nucleotide Archive (https://www.ebi.ac.uk/ena/program/home ) or SRA tool stash (utilizing the prefetch choice) may likewise be utilized to get the information through the Ubuntu terminal. We'll utilize the SRA tool compartment given by the NCBI to download the SRA documents to our PC. We can introduce the tool compartment utilizing the order line in the event that we are showing Linux to composing "sudo able introduce sra-tool stash". The prefetch order is utilized by the tool compartment to initially download the SRA document connected to the provided SRA ID. The sequencing information connected to the SRA ID can be downloaded from NCBI utilizing the "guidelines" contained in the SRA record.
We might make a message document with all the promotion numbers for every cell line since we have a ton of information to download. By doing this, we can rapidly download every one of the readings without a moment's delay.
2.2.3 Extracting the data.
We will have a record with the expansion ". sra" subsequent to downloading the records from the SRA tool compartment. We will utilize "fastq-dump" to change over the ". sra" documents into ". fastq "document" since this record must be transformed into a FASTQ document. The ". sra" records might be changed over into ". fastq" records or ". fastq.gz" records, which are compacted compress documents and won't occupy as much room. To change over the records into ". fastq" records, there's another choice we can download the "". fastq.gz" records straightforwardly from "The European Nucleotide Archive".The ENA offers a thorough and well-rounded resource for this essential source of biological knowledge by combining databases for raw sequence data, assembly details, and functional annotation. The provision of submission services—including interactive and programmatic submission tools—search services—including text and sequence similarity search tools—as well as data display and retrieval services—are essential components of the ENA.
To identify sequencing mistakes, PCR artefacts, or contamination, quality control for the raw reads entails analysing the sequence quality, GC content, the presence of adaptors, overrepresented k-mers, and duplicated reads. Acceptable amounts of duplication, k-mer, or GC content vary according on the experiment and the organism, but these values must be consistent among samples from the same research. We advise excluding outliers that have more than 30% disagreement. For these studies on Illumina reads, FastQC is a well-liked tool, whereas NGSQC may be used on any platform. As a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability. Creating FASTQ-format files with reads from an NGS platform, matching these reads to an annotated reference genome, and measuring gene expression are all part of the traditional process for RNA-Seq data . RNA-Seq analysis provides a special computational difficulty even though fundamental sequencing analysis tools are more widely available than before.

2.2.4 Quality control of the data.
FastQC is a tool intended to identify possible issues in datasets from high-throughput sequencing. One or more raw sequence files in fastq or bam format are subjected to a series of analysis, and a report summarising the results is generated. To use FastQC we will use the command “fastqc -- (“. fastq or. fastq.gz”)”.  We may utilise the Loop to execute on fastq files as we cannot do FastQC on every sample. When the job is finished,.html files will be present. On the HPC, you may access them by opening them in the  browser . The primary metrics to examine are: * Adapter content * Per base N content * Per base sequence quality. A QC report is often produced by most sequencers as part of their analysis pipeline, however this report is typically primarily concerned with finding issues that were caused by the sequencer. FastQC intends to deliver a QC report that can identify issues that either stem from the sequencer or the initial library material.
One of two modes can be used to run FastQC. It may be executed in two different ways: interactively as a standalone programme for the quick examination of a few FastQ files, or in a non-interactive mode that would be appropriate for incorporating into a broader analytic pipeline for the methodical processing of many files. Here are some fastqc reports of my cell lines.


 

The distribution of quality ratings across all positions in the read of “SRR12454227_1.fastq.gz”. This plot can alert us to whether there were any problems occuring during sequencing and whether we might need to contact the sequencing facility. The read's location is shown on the x-axis, and the quality ratings are shown on the y-axis. The plot's colour coding designates the high, medium, and low quality ratings.
 
Figure 4: Assessment of quality
Here we can notice a quality decline at the endpoints of the readings, which may be caused by phasing or signal decay. There are no more alarming symptoms, hence the facility's sequencing data is of high quality.
 
Figure 5: Duplicate Sequences
Duplicate sequences found in the library in large quantities. This plot can aid in identifying a low complexity library, which may be the consequence of insufficient starting material or excessive PCR amplification cycles. If this were a pilot project, we may change the number of PCR cycles, quantity of input, or amount of sequencing for future libraries. Normally, we don't do anything to address this in the analysis for RNA-seq. Due to the fact that the subset of data we are dealing with has the over-expression, there appear to be a lot of duplicated sequences in this study, but this is to be anticipated.
 
Figure 5: Per base sequence content
The distribution of GC over all sequences is seen in the "Per sequence GC content" graphic. In general, it's a good idea to note whether the centre peak's GC content matches the predicted GC percentage for the organism. Additionally, barring overrepresented sequences (sharp peaks on a normal distribution) or contamination with another organism, the distribution should be normal (broad peak). The steep peaks in this figure would suggest some sort of overrepresented sequence, which might be contamination or a gene that is overexpressed.
   
Figure 6: The GC content of another sample which “SRR12454229_1.fastq.gz”
The FastQC report gives us the report about the quality scores of the read which we want. 
We examine the data quality before applying mapping methods to identify which gene or transcript the reads came from. However, the mapping methods we utilise (salmon and STAR) are able to take into account for adaptor contamination, vector contamination, and low-quality bases at the ends of reads. The quality of the data is critical for identifying where it matches to on the genome or transcriptome. In order to align or map our raw reads to the reference genome or transcriptome, we must first identify any QC errors. Plot of readings per sample and GC percentage. Theoretical distributions are shown assuming a constant GC content across all reads. The GC content of each read should form a normal distribution for whole genome shotgun sequencing, with the apex of the curve occurring at the mean GC content for the sequenced organism. FastQC will declare a Fail if the observed distribution departs too much from the theoretical. This is likely to happen in a variety of circumstances, thus the assignment can be disregarded. In RNA sequencing, for instance, the mean GC content of the transcripts may be distributed more or less, leading the observed plot to be broader or narrower than an assumed normal distribution.
 
Another crucial thing is the "Overrepresented sequences" table, which lists the sequences (with at least 20 bases) that make up more than 0.1% of all sequences. This table makes it easier to spot contaminants like adapter or vector sequences. This table can assist in locating the source if the% GC content in the module above was incorrect. If the sequence is not identified as a recognised adapter or vector, BLAST can be used to assist identify it. List of sequences in the file that are more common than predicted. Only the first 50bp are taken into account. A sequence is considered overrepresented if it accounts for ≥ 0.1% of the total reads. To try to identify each overrepresented sequence, it is matched to a list of typical contaminants. While a small number of adaptor reads is not unusual, no one sequence should appear in DNA-Seq data with a frequency high enough to be reported. It's likely that certain transcripts in RNA-Seq data will be so numerous that they will appear as overrepresented sequences.

2.2.5 Building the genome index and Mapping:
It is important to build genome index because using Indices, the aligner can quickly and efficiently determine where a query sequence could have originated in the genome. In Linux, we can utilise this site (https://www.ensembl.org/Homo_sapiens/Info/Index) with the "wget" option to get the genome file which will be in “. fa" form, or we can just download the genome file directly via link. We'll get the human GRCh38 downloaded (Genome Reference Consortium Human Build 38). The human GRCh38 stand for alternate haplotypes and significantly affect how well we are able to identify and analyse genomic variation that is unique to populations that have these alternate haplotypes. Alternate or ALT contigs are used to represent common complicated variation, including variation at HLA loci. The Read Alignment Because many RNA-Seq reads map across splice junctions, mapping reads from DNA sequencing to the genome is much easier than mapping reads from RNA-Seq. Because they can't handle spliced transcripts, traditional read mapping methods like Bowti and BWA are really not advised for mapping RNA-Seq reads to the reference genome. A solution to this issue is to add sequences from exon-exon splice junctions obtained from well-known gene annotations to the reference genome. A "splicing-aware" aligner that can distinguish between reads with a brief insertion and reads aligning across an exon-intron border is the recommended method for mapping reads. Several splicing-aware mapping methods have emerged as RNA-Seq data have been more commonly utilised. A variety of splicing-aware mapping tools have been created especially for mapping transcriptome data as RNA-Seq data have been more often utilised. GSNAP is one of the most popular RNA-Seq alignment programmes. Performance, speed, and memory use benefits vary depending on the aligner. Based on these parameters and the general goals of the RNA-Seq investigation, the optimal aligner should be chosen. GENCODE's RNA-Seq Genome Annotation has spearheaded efforts to systematically assess the performance of RNA-Seq aligners and has discovered significant performance differences between alignment tools on a number of benchmarks, including alignment yield, basewise accuracy, mismatch and gap placement, and exon junction discovery.
Before building the genome it is necessary to download the Homo-sapiens annotation file
rule all:
    input:
        #get fa and gtf files
        "genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        "genome/Homo_sapiens.GRCh38.106.gtf.gz",
        #Get annotation GTF

rule get_genome_gtf:
    "Downloading Genome annotation file from Ensemble, Homo sapiens primary assembly (GRCh38)"
    output:
        gtf = "genome/Homo_sapiens.GRCh38.106.gtf.gz"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
        " && gunzip -k Homo_sapiens.GRCh38.106.gtf.gz 

To build the genome we used this pipeline:
rule hisat2_indexing:
    input:
        ref="genome.fa"
    output:
        touch("hisat2/makeidx.done")
    params:
        threads=20,
        idx="hisat2/genome_hisat2.idx"
    shell:
        """
        hisat2-build -p {params.threads} {input.ref} {params.idx}
        """


2.2.6 Alignment:
Peruses are for the most part planned to either a genome or a transcriptome. The extent of planned peruses, a general proportion of the by and large sequencing precision and the presence of sullying DNA, is an essential measurement for planning quality. Contingent upon the read mapper used, we guess that somewhere in the range of 70 and 90% of ordinary RNA-seq peruses will plan onto the human genome, with a sizeable piece of peruses relating to few comparative regions uniformly. We expect rather lower by and large planning rates when peruses are planned against the transcriptome on the grounds that peruses from unannotated records will be lost, and extensively more multi-planning peruses on the grounds that readings fall into exons that are shared by numerous record isoforms of a similar quality. Subsequent to investigating the nature of the information, we decide from which quality or record the peruses started from utilizing planning apparatuses. The nature of the information is significant while figuring out where it adjusts to on the genome or transcriptome, yet for the planning instruments we will utilize HISAT2 (progressive ordering for grafted arrangement of records) can represent bad quality bases at the finishes of readsIf our information were created utilizing a coordinated procedure, it would be useful to be aware for arrangement (HISAT2). This infers that utilizing specific strategies, it is plausible to decide if the readings are from similar strand as the qualities or records or from the contrary strand. Here qualities from the two strands are available at a similar site. You might know about the interaction used to create your information, yet assuming you are uncertain or uncertain, you might actually take a look at it. Accordingly, in the wake of taking note of any QC issues, we can involve our crude peruses for the arrangement or planning to the reference genome or transcriptome. We should "map" each read in your FASTQ document to a reference genome. Finding the area in the reference genome that most intently looks like the part of the mRNA record trapped in the read is alluded to as "planning." It is important to "map" each read in our FASTQ document to a reference genome. Finding the area in the reference genome that most intently looks like the piece of the mRNA record trapped in the read is alluded to as "planning". 
 We will use HISAT2 because the most accurate approach available, HISAT is the quickest system with little memory use. Even though HISAT has a lot of indexes, it only needs 4.3 GB of RAM. Any genome size, even ones with more than 4 billion bases, is supported by HISAT.
For aligning Pair-end reads we used:
hisat2 -p {params. threads} -x {params.idx} -1 {input.trim1} -2 {input.trim2} -S {output}
r1="Reads/{sample}_1.fastq.gz"
r2="Reads/{sample}_2.fastq.gz"

(SAMPLES,)=glob_wildcards("Reads/{sample}_1.fastq.gz")

rule all:
    input:
        expand("hisat2/{sample}.sam",sample=SAMPLES),
       

rule hisat2_Alignment:
    input:
        idxdone="hisat2/makeidx.done",
        trim1="Reads/{sample}_1.fastq.gz",
        trim2="Reads/{sample}_2.fastq.gz",
    output:
        "hisat2/{sample}.sam"
    params:
        idx="hisat2/genome_hisat2.idx",
        threads=5
    shell:
        """
        hisat2 -p {params.threads} -x {params.idx}  -1 {input.trim1} -2 {input.trim2} -S {output}
        """
        


For aligning Single-end reads we used:
hisat2 -p {params.threads} -x {params.idx}  -U {input.trim1} -S {output}
Here in - x {params.idx}. The basename of the file for the reference genome. The basename is the name of any of the record documents up to however excluding the last. Hisat searches for the predetermined list first in the flow registry, then, at that point, in the catalog determined in the HISAT_INDEXES climate variable.
-U: A collection of files containing unpaired reads that need to be aligned, separated by commas.
-S: To get the outputs in “. SAM” form
(SAMPLE,)=glob_wildcards("trimmedReads/{sample}.fastq")

rule all:
    input:
        expand("hisat2/{sample}.sam",sample=SAMPLE),

rule hisat2_Alignment:
    input:
        idxdone="hisat2/makeidx.done",
        trim1="trimmedReads/{sample}.fastq",
    output:
        "hisat2/{sample}.sam"
    params:
        idx="hisat2/genome_hisat2.idx",
        threads=20
    shell:
        """
        hisat2 -p {params.threads} -x {params.idx}  -U {input.trim1} -S {output}
        """

2.2.7 Conversion of reads: 
When utilising BAM files, the majority of functionality may be summed up as follows:
1.	BAM files are created by converting SAM files (samtools view)
samtools view -@ {threads} -b -o {output} {input}
(SAMPLE,)=glob_wildcards("hisat2/{sample}.sam")

rule all:
    input:
        expand("hisat2/{sample}_hisat2_sorted.bam",sample=SAMPLE),
        
                
rule convert_samtobam:
    input:
        "hisat2/{sample}.sam"
    output:
        "hisat2/{sample}_hisat2_sorted.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -b -o {output} {input}
        """
        

1.	By reference coordinates, BAM files are organised (samtools sort)
2.	Indexing of sorted BAM files (samtools index)
3.	Filters are applied to sorted, index BAM files according to location, flags, and mapping quality (samtools view with filtering options)
(SAMPLE,) = glob_wildcards("aligned/{sample}.sorted.bam")

rule all:
    input: 
        expand ("aligned/{sample}_hisat2_sorted.bam.bai"),
    

rule samtools_index:
    input:
        "aligned/{sample}_hisat2_sorted.bam",       
    output:
        "aligned/{sample}sorted.bam.bai",
    params:
        threads=20,
        idx="hisat2/genome_hisat2.idx"
    shell:
        "samtools index {input} {output}"

The main distinction among BAM and SAM documents is that BAM records are in paired document design, which can't be perused by people. Instead of SAM records, BAM documents are more modest and more viable for programming to work with, saving time and bringing down computation and capacity costs. Most programming that examinations adjusted readings hopes to ingest information in BAM design since arrangement information is typically consistently put away in BAM documents (frequently with a BAM record, to be examined later here).
2.2.8 Feature Counts:
FeatureCounts, a read summary tool that can count reads produced by either RNA or genomic DNA sequencing operations. FeatureCounts uses extremely effective feature blocking and chromosomal hashing algorithms. It uses a lot less computer memory and is noticeably faster than current techniques (by an order of magnitude for gene-level summarization). Read counts give a broad overview of the coverage for the relevant genetic feature. In particular, when several transcripts are being expressed from a single gene, gene-level counts from RNA-seq offer an overall summary of the gene's level of expression but do not differentiate between isoforms. Reads can normally be assigned to genes with a high degree of confidence, but determining the expression levels of specific isoforms is inherently more challenging due to the significant genomic overlap that exists between various isoforms of the same gene. Using data from reads clearly allocated to locations where isoforms differ, a variety of model-based techniques have been developed that aim to deconvolve the expression levels of separate transcripts for each gene using RNA-seq data. It offers a variety of settings suitable for various sequencing applications and works with single-end or paired-end reads. Genes, exons, promoters, gene bodies, genomic bins, and chromosomal positions are among the genomic characteristics that featureCounts counts mapped reads for. It is a highly effective general-purpose read summarization application. You can have paired or unpaired reads. When using paired reads, each read pair specifies a DNA or RNA fragment that is bookended by the two reads. Instead of counting reads in this situation, featureCounts will count fragments. In the event that paired readings do not appear in consecutive spots in the SAM or BAM file, featureCounts automatically sorts reads by name. SAM/BAM files and an annotation file with the chromosomal coordinates of features are inputs for featureCounts. Numbers of readings allocated to features are output (or meta-features). It also provides stat info for the overall summarization outcomes, including number of successfully assigned reads and number of reads that failed to be assigned owing to different reasons (these reasons are provided in the stat info) (these reasons are included in the stat info). This may be executed simultaneously on all SAM/BAM files.
A list of genomic features in Gene Transfer Format (GTF), General Feature Format (GFF), or Simplified Annotation Format (SAF), plus one or more files containing aligned reads (short or long reads) in either SAM or BAM format make up the data input to featureCounts. The input readings' format is automatically determined (SAM or BAM). Either GFF or SAF format can be used to specify the genetic characteristics. With only five mandatory fields for each feature (feature identification, chromosomal name, start position, end position, and strand), the SAF format is the easier to use. The bare minimum of data required for read quantification is provided by these five columns. According to the widely used Gene Transfer Format (GTF) revision of GFF, the feature IDs are considered to be unique in both formats.
FeatureCounts will naturally reorder the peruses to organize the peruses from a similar pair contiguous each other prior to counting them on the off chance that the information contains area arranged matched end peruses. At the point when name-arranged matched end input peruses are contradictory with featureCounts (for instance, due to the detailing of multi-planning results), featureCounts will furthermore reorder the peruses on a case by case basis. Prior to taking care of the peruses to featureCounts, clients can couple up the peruses utilizing the utility program fix that we offer.
A similar reference genome is utilized for both read arrangement and read counting. The BAM/SAM document contains data about each perused, remembering its beginning situation for the chromosome or contig/framework, the name of the reference chromosome or contig the read planned to, which contains point by point arrangement data about inclusions and cancellations comparative with the beginning position.
GTF considers the determination of hereditary qualities. As per the Gene Transfer Format (GTF) amendment of GFF, which is much of the time utilized, the element IDs are viewed as interesting. If strand-explicit data is given, featureCounts gives strand-explicit read counting. The aftereffects of perused planning normally give planning quality appraisals to planned peruses.
Here is the code we involved it for featureCounts:
featureCounts - T 8 - a Single/Homo_sapiens.GRCh38.106.gtf - g 'gene_id' - o featureCounts/Single.txt/*.bam
The featureCounts document would seem to be this: 
Figure 7: FeatureCounts file of cell line KYSE70.
2.2.9 Normalizing the FeatureCounts.
The constraints of restricted input material and different kinds of predisposition or clamor intrinsic in the sequencing system are tended to by information standardization, which is fundamental for single-cell sequencing. There are a few such standardization strategies, some of which depend on spike-in qualities, particles provided in referred to sums to go about as a standardization model's establishment. Explicit advantages of one methodology over another might be communicated relying upon the data that is accessible and the sort of information. We utilize two genuine informational collections with spike-in qualities and one mimicked examination to look at the viability of seven normalizing procedures accessible for single-cell sequencing.
We involved Quantile standardization in RStudio for standardization. Making two circulations genuinely indistinguishable is an interaction called quantile standardization. Without a reference dispersion, sort the conveyances as in the past, and afterward set the normal (frequently the arithmetical mean) of the circulations to quantile standardize them to each other. In this way, in all conditions, the most elevated esteem turns into the mean of the greatest qualities, and the second-most noteworthy worth turns into the mean of the second-most elevated values. 
Figure 8: The picture depicts the feature counts after normalizing.
2.2.10 Multiple Array Viewer:
Multiple Array Viewer (MeV) is a desktop Java programme that offers an easy-to-use graphical user interface and sophisticated analysis of gene expression data. Here, we present a substantial improvement to MeV that makes it possible to analyse RNA-Seq data using these reliable, strong techniques. We further mention that many RNA-Seq-specific capabilities have been added to MeV, addressing the different analytic needs for this data format compared to conventional gene expression data. In addition to differential expression discovery and functional annotation enrichment detection based on published methodologies, these tools also feature automated conversion functions from raw count data to processed RPKM or FPKM values.  
First, we load for our data into MeV, the data is basically in text delimited form.
 
Figure 9: Loading the Data

A unique test for identifying differentially expressed genes in several replicates is called Rank Products. This approach varies from many other methods in that it calculates rank products, a quicker and easier way, as opposed to using a complex statistical model. We used Rank Prod analysis because it can dramatically minimise the number of duplicate tests necessary to provide trustworthy findings and is effective with excessively noisy data. 
(A) Three experimental designs are supported by MeV's RP.
1.	One-class, finds the genes that are highly up- or down-regulated within the included group and is commonly conducted on two-colour data. Uncheck the box next to a sample's name in the left pane of the one-class screen to remove it from the analysis.
2.	When samples belong to one of two groups and the participants vary between the two groups, the two-class unpaired design (analogous to a between subjects t-test). The t-test dialogue box resembles the startup dialogue box.
3.	When samples are two-class paired, they are not only divided into two groups but also paired one-to-one (for example, when measuring the gene expression of a group of patients, measurements are made before (Group A) and after (Group B) each person receives a pharmacological treatment).
4.	Specified Genes, MeV supports two methods for sample selection, as do the majority of other modules: the choice of each button individually the sample analysis tabs on the cluster selection tabs to allocate your samples. Unassigned or unchecked samples will not be considered in the study.
(B) Data
 
Figure 10: The result after loading the Data.
Only clusters that are "redder" are worthy of attention in a big dataset, which is likely to include many clusters that fall within the size range. To save time, it is crucial to give the proper colour gradient limits. A subsequent part will discuss how to change the colour gradient. Simply place the mouse over the cluster's root to display a pop-up window and examine information about each cluster.
  (C) Parameters 
 
Figure 12: Setting the parameters
T-values are determined for each gene, and p-values are generated either from permutations of the data for each gene or from the theoretical t-distribution. A gene's p-value is compared to the user-specified critical p-value or alpha to see if there is a statistically significant difference between the two groups. Alternatively, the p-values are adjusted to account for multiple testing.
(D) Extracting the up-regulated, Down-regulated and Non-significant genes.
To extract Up and Down Regulated Genes we used the software MEV Multiple Experiment Viewer). With the help of Rank prod analysis, we found the genes. The p-Value cut-off was set to 0.001. The results of the analysis will look this:
 
Figure 13: Cluster Information of all the regulated genes.

2.2.11 Functional analyses:
To analyze and break down a gathering of qualities to recognize and concentrate on the qualities that are engaged with the guideline of specific natural pathways is known as practical investigation or natural pathways examination.
For instance, a quality name, overlap change, and various testing-remedied P-an incentive for quality articulation research; or a SNP, chances proportion, and P-an incentive for an association with a sickness aggregate from a hereditary affiliation review — are the normal results of conventional bioinformatics investigation. We utilized Package UpsetR is RStudio to examine our UP and DOWN controlled qualities. The utilitarian examination is displayed in results.
2.2.12 Packages
 UpSet analysis presents the data in a tabular style with bar graphs. You may include as many items as you like in each category and quickly identify the points where they meet. It is implemented via a R package called UpSetR.The greatest tool for identifying patterns in large, complicated datasets with many attributes is UpSetR. It accomplishes this by combining data points with a lot of the same values across many attributes. In other words, UpSetR identifies sets with the biggest intersection.It works best with variables that are already binary (i.e., have two categories) or that can be made binary. It is significant to remember that every category feature may be converted into a binary feature using one-hot encoding.
Subread is a group of high-performance software tools for processing data from next-generation sequencing are included in the Subread/Rsubread packages. These packages contain the utility programmes exact SNP calling, Subread aligner, Subtunic aligner, Subindel long indel identification, featureCounts read quantification, and others. The featureCounts tool is made to associate genomic characteristics like genes, exons, and promoters with mapped reads or fragments (paired-end data). This little read counting application can count both gDNA-seq and RNA-seq reads for genomic characteristics.
Euler creates area-proportional euler diagrams that use circles or ellipses to represent set relationships (intersections, unions, and disjoints). The absence of any set interactions is a condition for Euler diagrams, unlike Venn diagrams (whether they are empty or not). In other words, depending on the input, eulerr will occasionally create Venn diagrams and occasionally not.
Euler was written from scratch, employs alternative optimizers, returns statistics seen in veneer and Euler APE, and accepts a variety of inputs and conditioning on extra variables. It is built on the enhancements to vene eulerr that Ben Fredrickson offered in venn.js. Additionally, it can represent set relationships for any quantity of sets with ellipses.


3. Results:
3.1 MultiQC:
Many popular bioinformatics tools are supported by MultiQC, but research groups will inevitably need their own custom scripts or other customizations. To account for this, MultiQC is designed in a way that makes it simple to integrate new code into its operations. Code hooks provide third-party plugins access to and control over a program's internal operations. Modules, templates, and plugins may be retained in a distinct code base and nevertheless run as a part of the main MultiQC application thanks to the usage of Python setup tools entry points.
A single self-contained HTML report produced by MultiQC may be shared and accessed in any current web browser. Using the JavaScript charting package HighCharts, reports produce charts. Plots may be zoomed in and out of, and some are interactive. Using a report toolkit, samples may be renamed, hidden, and highlighted. A variety of publication-ready formats are available for plot export.
3.1.1 MultiQC report of LnCaP:
 
We performed MultiQC report of every cell line and this is the report of LnCaP cell line. As we can see that the report gives us the percentage of Duplicates, the percentage of GC content, the length of sequences.

Sequence Counts:
 
 
This is the sequence counts plot where the blue colour represents the Unique Reads and the black colour represents the Duplicate Reads. It also shows how many percentages of reads are unique and how many are duplicate.
3.1.2 MultiQC report of KYSE150:
 
3.1.3 MultiQC report of KYSE70:
 
3.1.4 MultiQC report of KYSE450:
 



3.2 UpsetR analysis:  
1. (A) KYSE70 (Pair-end):

 
No genes are intersecting in the figure, indicating that there are no shared genes among up-regulated, down-regulated, and non-significant genes. While there are 421 down-regulated genes and 550 up-regulated genes, we can observe that there are a large number of 60584 non-significant genes.

(B) KYSE70 (Single-end):

 
The given figure represents the number of Up-regulated, Down-regulated and non-significant genes.
(C) KYSE150 Single-end
     
(D) KYSE450 Single-end:
   
(E) LnCaP Pair-end
 
(F) LnCaP Single-end:
 
3.3 UpsetR analysis of Single-end Up-regulated genes for all cell lines.

  
As we can see, LnCap has 7091 up-regulated genes, which is the highest number, and KYSE70 has 97, which is the lowest number. The picture also demonstrates the 161 genes that are shared by all cell lines. Additionally, it displays the junction of two, three, or even more cell lines. For instance, there are just 7 genes shared by KYSE70 and KYSE150.
 
The single-end up-regulated genes' euler plot. LnCaP, KYSE150, and KYSE450 are intersecting at one point, but KYSE70 is not intersecting, indicating that these cell lines share 155 genes in common.
Extracted Gene ID:
The results are in html form too with the extracted gene ids available. I am pasting the link here, please kindly check it. Thank you.
Results
3.4 Functional analyses of Consistent genes:
The extraction of consistent genes was done using the Python. The analyses were done using http://www.webgestalt.org/ 
3.4.1 KYSE450 Single-end
(A) 1. Up regulated Genes
 
After mapping the up-regulated genes in Gene otology (Functional Database). We get three types of bar chart where it shows the function of genes. There are three categories namely:
1. Biological process.
2. Cellular Component.
3. Molecular Function.

 

As we can see we there are two genes which were able to mapped out of 70.

In order to repair DNA double-strand breaks, this gene's protein collaborates with DNA ligase IV and DNA-dependent protein kinase. This protein is essential for both the completion of V(D)J recombination and non-homologous ends joining. Endocrine dysfunction, microcephaly, and low stature can all result from mutations in this gene.
2. Down-regulated Genes: 
The gene set is 926 but the mapping rate is very low. In addition, the mapping input was 57 but it mapped only 11 genes.

 
3. Gene Signature:
 
3.4.2 KYSE70 Single-end.
(A) 1. KYSE70 Down-regulated Genes: -
 
 
As we can see the biological processes are very low in KYSE70 Down-regulated Genes. The mapping of genes is comparatively very low.

2. KYSE70 Single-end Up-regulated Genes:
 
The biological processes of KYSE70 Up-regulated Genes is very high comparatively to the down-regulated genes. In addition, more than half of the genes possess the metabolic genes and the protein bindation is also very high amongst the genes. We have 3 bar graphs. We can see there are 317 genes which are functioning and the least functioning gene is in Cellular component which is microbody
3. Gene Signatures of both regulated genes (Single-end):
 
The enrichment results of both up and down regulated genes shows that the translation initiation of all the genes is highest while the lowest is the cellular protein localization. Meanwhile the overlapping of all genes is of total 33 genes. The overlapping genes are as follows:
  
This are the gene-signature when we combine both regulated genes fromKYSE70. Most of them are catalysing proteins. For example: The ribosomal protein L11P family includes the protein in question. In the cytoplasm is where it is found. The 26S rRNA is directly bound by the protein. The U65 snoRNA, which is situated in the fourth intron of this gene, co-transcribes with it. This gene is distributed throughout the genome as many processed pseudogenes, as is normal for genes encoding ribosomal proteins.
3.4.2 KYSE70 Pair-end
(A) 1.Up-regulated Genes:
 
(B) 2.Down-regulated Genes:
 
The FDR is less than 0.05. Meanwhile the genes mapped were only 2.
  
If we look at the bar chart the process are way lower, keeping the affinity in mind we can say the the processes are way lower 
 
(C) Gene Signature Pair-end:
 
The overlapping genes are very less if we compare it to the gene sets.



3.4.3 LnCap:
 (A) 1. Down-regulated genes Pair-end:
 
 
The genes function least in Cellular protein localization and highest in translational initiation.
 

As compare to other Cell lines the overlapping of LnCap is higher which is 106. The mapping rate of this cell lines is also more than other cell lines.
2. Up-regulated Genes Pair-end:
 
 
3. Gene Signature of both regulated Genes Pair-end: 
 
(B) Single-end 
1. Up-regulated Genes:
 
2.Down-Regulated Genes:
 
3. Gene Signature Single-end:
 



3.4.4 KYSE150:
(A) Single-end Up-regulated genes:
 
The genes mapped are
 
(B) Single-end Down-regulated genes:
 
 
(C) Gene-signature:
 

4.Discussion:
In KYSE70, in Single Down-regulated genes we found that the genes were family of olfactory receptor family 6 subfamily T member 1. The genes were related to a neural response that results in the experience of a scent, odorant molecules in the nose interact with olfactory receptors. The olfactory receptor proteins are derived from a single coding-exon gene and belong to the broad family of G-protein-coupled receptors (GPCR). Olfactory receptors are in charge of odorant signal detection and G protein-mediated transduction and share a 7-transmembrane domain structure with several neurotransmitter and hormone receptors. The biggest gene family in the genome is that of the olfactory receptor. Meanwhile, in Pair-end the extracellular matrix is broken down by proteins in this family during both healthy physiological processes including embryonic development, reproduction, and tissue remodelling as well as pathological activities like arthritis and metastasis. To create the mature protease, the encoded preproprotein undergoes proteolytic processing. Both soluble and insoluble elastin are destroyed by this protease. The gene are been linked to lung function and chronic obstructive pulmonary disease, and it may contribute to the development of aneurysms (COPD).KYSE450, Single-end genes cell line produces a DNA ligase protein, which binds single-strand breaks to double-stranded polydeoxynucleotides in an ATP-dependent process. This protein is necessary for DNA double-strand break (DSB) repair via nonhomologous end joining and V(D)J recombination (NHEJ). This protein interacts with DNA-dependent protein kinase in addition to forming a complex with X-ray repair cross complementing protein 4 (XRCC4) (DNA-PK). It is understood that XRCC4 and DNA-PK are both necessary for NHEJ.KYSE150, single-end, a member of the protein family known as receptor tyrosine kinases and a by-product of the proto-oncogene MET. Alpha and beta subunits of the encoded preproprotein are produced during proteolytic processing and joined by disulphide bonds to create the mature receptor. The M10 peptide, which has been found to lessen lung fibrosis, is created by further processing of the beta subunit. Hepatocyte growth factor, the receptor's ligand, binding results in dimerization and activation, which affects cellular survival, embryogenesis, and cellular migration and invasion,a cell migration regulator. The encoded protein seems to participate in the RhoA-Dia1 signal transduction pathway (Ras homolog gene family, member A). Differentially spliced transcripts.LnCaP, A serine/threonine protein kinase called casein kinase II phosphorylates acidic proteins like casein. It is implicated in a number of biological functions, including circadian rhythm, apoptosis, and cell cycle regulation. An alpha, an alpha-prime, and two beta subunits make up the tetramer form of the kinase. While the beta subunits are autophosphorylated, the alpha subunits have the catalytic activity. The alpha subunit is represented by the protein that this gene encodes. For this gene, many transcript variants have been discovered that encode various protein isoforms. Overall, we found consistent genes from four cell lines. The Highest number of genes which are found was from LnCap and the Lowest number of genes were found from KYSE450 which was 2. However, this represents the genes we were able to extract from the cell lines.

5.References:
1.	https://thejns.org/view/journals/j-neurosurg/90/3/article-p533.xml
2.	Tawfik, G., Dila, K., Mohamed, M., Tam, D., Kien, N., Ahmed, A. and Huy, N., 2019. A step by step guide for conducting a systematic review and meta-analysis with simulation data. Tropical Medicine and Health, 47(1).
3.	(https://iiste.org/Journals/index.php/ADS/article/view/52906, 2020)
4.	Webgestalt.org. 2022. WebGestalt (WEB-based GEne SeT AnaLysis Toolkit). [online] Available at: <http://www.webgestalt.org/> [Accessed 29 July 2022].
5.	Ncbi.nlm.nih.gov. 2022. XRCC4 X-ray repair cross complementing 4 [Homo sapiens (human)] - Gene - NCBI. [online] Available at: <https://www.ncbi.nlm.nih.gov/gene/?term=7518> [Accessed 29 July 2022].
6.	Clayden, J., Vallverdú, L. and Helliwell, M., 2006. Conformational communication between the Ar–CO and Ar–N axes in 2,2′-disubstituted benzanilides and their derivatives. Org. Biomol. Chem., 4(11), pp.2106-2118.
7.	Kirkuk University Journal-Scientific Studies, 2020. https://kujss.iraqjournals.com/pdf_166170_8dd024058ce4abb6c364bec514cecef8.html. 15(2), pp.1-16.
8.	HISAT2. 2022. HowTo | HISAT2. [online] Available at: <http://daehwankimlab.github.io/hisat2/howto/> [Accessed 29 July 2022].
9.	2022. [online] Available at: <https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19> [Accessed 29 July 2022].
10.	Chung, M., Bruno, V., Rasko, D., Cuomo, C., Muñoz, J., Livny, J., Shetty, A., Mahurkar, A. and Dunning Hotopp, J., 2021. Best practices on the differential expression analysis of multi-species RNA-seq. Genome Biology, 22(1).
11.	Gaming Law Review and Economics, 2016. LAS VEGAS SANDS CORP., a Nevada corporation, Plaintiff, v. UKNOWN REGISTRANTS OF www.wn0000.com, www.wn1111.com, www.wn2222.com, www.wn3333.com, www.wn4444.com, www.wn5555.com, www.wn6666.com, www.wn7777.com, www.wn8888.com, www.wn9999.com, www.112211.com, www.4456888.com, www.4489888.com, www.001148.com, and www.2289888.com, Defendants. 20(10), pp.859-868.
12.	Cross, J., 2006. MEDLINE, PubMed, PubMed Central, and the NLM. Editors' Bulletin, 2(1), pp.1-5.
13.	https://www.uniprot.org/
14.	Luiz de Freitas Vieira, J. and Almeida Có, M., 1997. https://sobraep.org.br/artigo/sistema-retificador-inversor-com-corrente-pulsada-no-barramento-cc-para-acionamento-de-motores-de-inducao/. Eletrônica de Potência, 2(1), pp.35-42.
15.	Lonstein, J., 2016. Update on NLM/PubMed Status. Spine Deformity, 4(3), p.165.
16.	https://pubmed.ncbi.nlm.nih.gov/
17.	Yang Liao, Gordon K. Smyth, Wei Shi, featureCounts: an efficient general purpose program for assigning sequence reads to genomic features, Bioinformatics, Volume 30, Issue 7, 1 April 2014, Pages 923–930, https://doi.org/10.1093/bioinformatics/btt656
18.	Kukurba KR, Montgomery SB. RNA Sequencing and Analysis. Cold Spring Harb Protoc. 2015 Apr 13;2015(11):951-69. doi: 10.1101/pdb.top084970. PMID: 25870306; PMCID: PMC4863231.
19.	Koch CM, Chiu SF, Akbarpour M, Bharat A, Ridge KM, Bartom ET, Winter DR. A Beginner's Guide to Analysis of RNA Sequencing Data. Am J Respir Cell Mol Biol. 2018 Aug;59(2):145-157. doi: 10.1165/rcmb.2017-0430TR. PMID: 29624415; PMCID: PMC6096346.
20.	Wikis.utexas.edu. 2022. Mapping with HISAT2 - Bioinformatics Team (BioITeam) at the University of Texas - UT Austin Wikis. [online] Available at: <https://wikis.utexas.edu/display/bioiteam/Mapping+with+HISAT2> [Accessed 19 August 2022].
21.	Tonelli C, Morelli MJ, Bianchi S, Rotta L, Capra T, Sabò A, Campaner S, Amati B. Genome-wide analysis of p53 transcriptional programs in B cells upon exposure to genotoxic stress in vivo. Oncotarget. 2015 Sep 22;6(28):24611-26. doi: 10.18632/oncotarget.5232. PMID: 26372730; PMCID: PMC4694782.
22.	Wikis.utexas.edu. 2022. Mapping with HISAT2 - Bioinformatics Team (BioITeam) at the University of Texas - UT Austin Wikis. [online] Available at: <https://wikis.utexas.edu/display/bioiteam/Mapping+with+HISAT2> [Accessed 19 August 2022].
23.	Chung, M., Bruno, V.M., Rasko, D.A. et al. Best practices on the differential expression analysis of multi-species RNA-seq. Genome Biol 22, 121 (2021). https://doi.org/10.1186/s13059-021-02337-8.
24.	Mary Piper, R., 2022. Quality control: Assessing FASTQC results. [online] Introduction to RNA-Seq using high-performance computing - ARCHIVED. Available at: <https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html> [Accessed 19 August 2022].
25.	Conesa, A., Madrigal, P., Tarazona, S. et al. A survey of best practices for RNA-seq data analysis. Genome Biol 17, 13 (2016). https://doi.org/10.1186/s13059-016-0881-8
26.	Anders, S., Huber, W. Differential expression analysis for sequence count data. Genome Biol 11, R106 (2010). https://doi.org/10.1186/gb-2010-11-10-r106
27.	Yalamanchili HK, Wan YW, Liu Z. Data Analysis Pipeline for RNA-seq Experiments: From Differential Expression to Cryptic Splicing. Curr Protoc Bioinformatics. 2017 Sep 13;59:11.15.1-11.15.21. doi: 10.1002/cpbi.33. PMID: 28902396; PMCID: PMC6373869.
28.	Cuiping Li, Dongmei Tian, Bixia Tang, Xiaonan Liu, Xufei Teng, Wenming Zhao, Zhang Zhang, Shuhui Song, Genome Variation Map: a worldwide collection of genome variations across multiple species, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D1186–D1191, https://doi.org/10.1093/nar/gkaa1005
29.	 Decodon.com. 2022. [online] Available at: <https://www.decodon.com/files/doc/MeV_Manual_4_0.pdf> [Accessed 19 August 2022].
30.	: Documentation.partek.com. 2022. Normalize to baseline - Flow Documentation - Partek® Documentation. [online] Available at: <https://documentation.partek.com/display/FLOWDOC/Normalize+to+baseline> [Accessed 19 August 2022].
31.	: Ncbi.nlm.nih.gov. 2022. ATP5MF ATP synthase membrane subunit f [Homo sapiens (human)] - Gene - NCBI. [online] Available at: <https://www.ncbi.nlm.nih.gov/gene/?term=9551> [Accessed 19 August 2022].
32.	F., docs, R., browser, R. and Ellipses, e., 2022. eulerr: Area-Proportional Euler and Venn Diagrams with Ellipses version 6.1.1 from CRAN. [online] Rdrr.io. Available at: <https://rdrr.io/cran/eulerr/> [Accessed 19 August 2022].
33.	Metagenomics-workshop.readthedocs.io. 2022. Determine gene coverage in metagenomic samples — Metagenomics Workshop SciLifeLab 1.0 documentation. [online] Available at: <https://metagenomics-workshop.readthedocs.io/en/stable/comparative-functional-analysis/genecoverage.html> [Accessed 19 August 2022].
34.	Metagenomics-workshop.readthedocs.io. 2022. Functional annotation — Metagenomics Workshop SciLifeLab 1.0 documentation. [online] Available at: <https://metagenomics-workshop.readthedocs.io/en/stable/comparative-functional-analysis/annotation.html> [Accessed 19 August 2022].
35.	Disc.marseille.inserm.fr. 2022. Next Generation Sequencing data analysis with Subread and IGV. Example of a real-world SNP analysis.. [online] Available at: <https://disc.marseille.inserm.fr/_static/techdays_20180220_subread_GB.html#slide2> [Accessed 19 August 2022].
36.	Liao Y, Smyth GK, Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res. 2013 May 1;41(10):e108. doi: 10.1093/nar/gkt214. Epub 2013 Apr 4. PMID: 23558742; PMCID: PMC3664803.
37.	 Medium. 2022. UpSetR is the Greatest Set Visualization Since the Venn Diagram. [online] Available at: <https://towardsdatascience.com/upsetr-is-the-greatest-set-visualization-since-the-venn-diagram-8bccbaef698a> [Accessed 19 August 2022].
38.	Notes, T. and FAQ, F., 2022. FastQC Tutorial & FAQ. [online] Rtsf.natsci.msu.edu. Available at: <https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/> [Accessed 19 August 2022].


