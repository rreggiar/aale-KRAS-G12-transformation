# Tutorial for DNAnexus Platform (web)

All test samples and genome data are shared on our public DNAnexus project. You don't have to download any data for testing our pipeline on DNAnexus platform.

There are two methods to run our pipeline on DNAnexus.

1) [Building your own DX workflow from `atac.wdl` with dxWDL (CLI)](tutorial_dx_cli.md)
2) Using a pre-built DX workflow on our public DX project (Web UI)

This document describes instruction for the item 2).

1. Sign up for a [DNAnexus account](https://platform.DNAnexus.com/register).

2. Create a new [DX project](https://platform.DNAnexus.com/projects) by clicking on "+New Project" on the top left.

3. Move to one of the following workflow directories according to the platform you have chosen for your project (AWS or Azure). These DX workflows are pre-built with all parameters defined.

* [AWS test workflow](https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/ATAC-seq/workflows): Use `[LATEST_VER]/test_ENCSR356KRQ_subsampled`.
* [Azure test workflow](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/ATAC-seq/workflows): Use `[LATEST_VER]/test_ENCSR356KRQ_subsampled`.

4. Copy it to your project by right-clicking on the DX workflow `atac` and choose "Copy". 

5. Choose your project and create a folder for the test run by clicking on the "Folder+" icon.

6. Click on "Copy into this folder" on the bottom left.

7. Move to the target folder and click on the DX workflow `atac`.

9. Specify an output directory by clicking "Workflow Actions" on the top right. Click on "Set output folder" and choose an output folder.

10. Click on "Run as Analysis..." and you will be automatically redirected to the "Monitor" tab.

11. It will take about an hour. You will be able to find all outputs on your output folder. Final QC report (`qc.html`)/JSON (`qc.json`) will be found on it.

11. See full specification for [input JSON file](input.md).


## Extras for advanced users

1. DNAnexus allows only one copy of a workflow per project. The example workflow in the previous section is pre-built for the subsampled test sample [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/) with all parameters defined already.

2. Choose your main platform (AWS or Azure). Move to [ENCODE ATAC-seq pipeline repository for AWS](https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/ATAC-seq/workflows) or [ENCODE ATAC-seq pipeline repository for Azure](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/ATAC-seq/workflows).

3. Choose a folder with the latest available version.

4. Copy one of the following workflows according to the platform you have chosen for your project.
> **IMPORTANT**: Make sure that you have chosen a correct platform (AWS or Azure) for your project.

  * general: General workflow without pre-defined reference genome.
  * hg38: Worfklow with pre-defined hg38 reference genome.
  * hg19: Worfklow with pre-defined hg19 reference genome.

5. Click on the DX workflow `atac`.

6. Specify your input files (FASTQs, BAMs, TAG-ALIGNs, ...) on the top left. For example, click on the item "fastqs_rep1_R1" and choose your R1 FASTQ file for replicate 1. See details [here](input.md) for other input types.

7. Choose a reference genome. See details [here](input.md).

8. Click on "Run as Analysis..." and you will be automatically redirected to the "Monitor" tab.

