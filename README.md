# call_peaks
The script uses the uncontinuous score from a bigWig file to estimate footprints within the peaks of interest. A sliding window algorithm is used to look through the peaks. A region with significant high signal will be saved as a footprint. To estimate the significance of a footprint the signal is compared with the mean of the signals within the window. After the run the footprints are saved in a .bed file.
## Required input
There are two required input files for the peak calling. They are:
* --bigwig a bigWig file with the signal;
* --bed a corresponding .bed file with peaks of interest.

## Output files
The user can set the name for the output .bed file using the parameter --output_file. By default the output file is called _call_peaks_output.bed_ and is saved to the working directory. For each footprint an individual score is found. The score shows the quality of the footprint and is found as the mean of the signals from bigWig file for the corresponding positions. The output .bed file has 8 columns:
* _chr_ for the chromosom;
* _start_ for the start position of the footprint;
* _end_ for the end position of the footprint;
* _name_ for the unique name of the footprint;
* _score_ for the quality of the footprint;
* _len_ for the length of the footprint;
* _max_pos_ for the position where the highest signal from the bigWig file is (if there are several positions with the highest score the position in the middle of those will be set as the max_pos);
* _bonus_info_ for the additional information from the input .bed file.

During the run, a log-file is written. The log-file is called _call_peaks_log.txt_ and is saved to the working directory as well. If the user does not want to receive the messages about the run in the terminal, the parameter --silent can be used to force the script write the messages only to the log-file.

## Optional parameters
There are also optional parameters that could be set to refine the search:
* --window_length. This parameter sets the length of a sliding window. By default the length 200 bp was taken as the smallest peak of the test data was 200 bp long. 
* --step. This parameter sets the number of positions to slide the window forward. By default the step of a length 100 bp was taken.
* --percentage. By default each signal from the bigWig file is compared to the threshold which is the mean of signals within a window (or within a peak, if it is smaller than the chosen window_length). Though the user has a possibility to move this threshold using the parameter --percentage. For example --percentage 10 will add 10% of the found threshold within a window and set it as a new threshold to compare the signal to. The default percentage is set to 0.

Changing the optional parameters can lead to varying the number of found footprints and their length. 
