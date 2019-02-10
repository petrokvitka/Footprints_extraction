# Footprints extraction

There is a number of tools to find peaks from different file formats. Still there was no suitable tool to extract the peaks from a file with uncontinuous signal. Therefore the _footprints_extraction.py_ was developed. The signal in the input file in .bigWig format is used to estimate peaks or in other words footprints. A sliding window algorithm is used to look through the signal from the input file. A region with significant high signal will be saved as a footprint. To estimate the significance of a footprint the signal is compared with the mean of the signals within the window. After the run the footprints are saved in a .bed file.

# Installation

Please mention that [CONDA](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is needed for the installation of the Footprints extraction working environment.

1. Clone the directory
```bash
git clone https://github.com/petrokvitka/Footprints_extraction
```
2. Switch to the directory
```bash
cd Footprints_extraction
```
3. Create the needed environment from the file fp_extract.yaml
```bash
conda env create --file fp_extract.yaml
```
4. Activate the environment
```bash
source activate fp_extract
```

## Required input
There are two required input files for the peak calling. They are:
* --bigwig a bigWig file with the signal;
* --bed a corresponding .bed file with peaks of interest.

## Output files
The user can set the name for the output .bed file using the parameter --output_file. By default the output file is called _footprints_extraction_output.bed_ and is saved to the working directory. While using the whole pipeline, the default name will be used for the output file to save it as an intermediate result.

For each footprint an individual score is found. The score shows the quality of the footprint and is found as the mean of the signals from bigWig file for the corresponding positions. The output .bed file has 8 columns and a header marked with #. The header identifies the names for all columns:
```
#chr start end name score len max_pos bonus_info
```
* _chr_ for the chromosom;
* _start_ for the start position of the footprint;
* _end_ for the end position of the footprint;
* _name_ for the unique name of the footprint;
* _score_ for the quality of the footprint;
* _strand_ for the strand of the footprint;
* _len_ for the length of the footprint;
* _max_pos_ for the relative position within the footprint, where the highest signal from the bigWig file is (if there are several positions with the highest score the position in the middle of those will be set as the max_pos);
* _bonus_info_ for the additional information from the input .bed file.

If some information could not be defined, the point . will be written in the corresponding row and column. For example, if there was no _bonus_info_ provided from the input .bed file, the point . could be found at the output file of the footprints_extraction.

During the run, a log-file is written. The log-file is called _footprints_extraction.log_ and is saved to the working directory as well. If the user does not want to receive the messages about the run in the terminal, the parameter --silent can be used to force the script write the messages only to the log-file.

## Optional parameters
There are also optional parameters that could be set to refine the search:
* --window_length. This parameter sets the length of a sliding window. By default the length 200 bp was taken as the smallest peak of the test data was 200 bp long. 
* --step. This parameter sets the number of positions to slide the window forward. By default the step of a length 100 bp was taken.
* --percentage. By default each signal from the bigWig file is compared to the threshold which is the mean of signals within a window (or within a peak, if it is smaller than the chosen window_length). Though the user has a possibility to move this threshold using the parameter --percentage. For example --percentage 10 will add 10% of the found threshold within a window and set it as a new threshold to compare the signal to. The default percentage is set to 0.
* --min_gap. This parameter sets the minimal allowed distancee between two neighbour footprints. All footprints that are nearer to each other than this distance will be merged. By default the _min_gap_ is set to 6 bp.

Changing the optional parameters can lead to varying the number of found footprints and their length. The smaller the _--percentage_ is set, the longer the found footrpints will be, but the quality of each footprint will be lower. The length of a window has also an impact on the length of found footprints. The bigger the _--window_len_ parameter ist set, the longer the found footprints will be, though this change is not as remarkable as while changing the _--percentage_ parameter.

There is a possibility for a user to set the max allowed number of base pairs in between two footprints with the help of a parameter _--min_gap_. By default 6 bp are allowed. That means, all footprints, that have a smaller number of footprints in between, will be merged. The score for a merged footprint is calculated as a mean of scores of both original footprints. If the user doesn't want to merge the footprints, the _--min_gap_ should be set to 1.

## Example data

The files in [example](./example) are based on the Buenrostro's [research](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374986/). The file in .bigWig format represents the uncontinuous signal from the ATAC-seq and the file in .BED format contains the corresponding peaks. The run with example data can be started using following command:

```bash
python footprints_extraction.py --bigwig example/buenrostro50k_chr1_fp.bw --bed example/buenrostro50k_chr1_peaks.bed --output_file example_output/output.bed
```

The output file containing the footprints, as well as the logfile can be found in the [example_output](./example_output) folder.
