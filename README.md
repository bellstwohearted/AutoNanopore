# AutoNanopore
**An automated adaptive and robust method to locate translocation events in solid-state nanopore current traces.**

The article has been published in *ACS Omega*. [Click here to read it](https://pubs.acs.org/doi/10.1021/acsomega.2c02927)

**Abstract** Solid-state nanopore sequencing has shown impressive performances in several research scenarios but is still challenging, mainly due to the ultrafast speed of DNA translocation and significant noises embedded in raw signals. Hence, event detection, aiming to locate precisely these translocation events, is the fundamental step of data analysis. However, existing event detection methods use either a user-defined global threshold or an adaptive threshold determined by the data, assuming the baseline current to be stable over time. These disadvantages limit their applications in real-world application scenarios, especially considering that the results of different methods are often inconsistent. In this study, we develop an automated adaptive method called AutoNanopore, for fast and accurate event detection in current traces. The method consists of three consecutive steps: current trace segmentation, current amplitude outlier identification by straightforward statistical analyses, and event characterization. Then we propose ideas/metrics on how to quantitatively evaluate the performance of an event detection method, followed by comparing the performance of AutoNanopore against two state-of-the-art methods, OpenNanopore and EventPro. Finally, we examine if one method can detect the overlapping events detected by the other two, demonstrating that AutoNanopore has the highest coverage ratio. Moreover, AutoNanopore also performs well in detecting challenging events: e.g., those with significantly varying baselines.

------------------

We performed the analysis by running the python scirpt, a sample abf file and the corresponding output are also given.

```Python
python AutoNanopore.py --file_path dataset/1.abf --output_path dataset/output/
```

The output will be

```Bash
Processing started
sampling rate is 250000
total number of data point is 75000000
total length of the file is 300.0s
10000 slices generated
9975 valid slices selected
Optimizing ...
Processing finished within 13.242s, 80 events detected.
```

Finally, a file named `1.csv` will be generated, including the characterizations of all the events detected.

There are five parameters:

* `--file_path`, the path of the raw signal `.abf` file, it is required for the program to run.

* `--output_path`, the path for the output `.csv` file to be saved,  it is required for the program to run.

* `--signal_direction`, the direction of event detection, the default is `0`, meaning that the events are detected in positive-going, set it to `1` if the events are detected in negative-going direction.

* `--theta`, the multiple of IQR in amplitude outlier identification, the default is `1.5`. Note that this is the ***core*** parameter of the method. While increasing the value of `theta` would lead to detection of high-quality events, this may also result in missing of some events, so `theta = 1.5` is suitable, according to the statistical principles. Moreover, we strongly suggest that the value should not be smaller than `1.5`, or the risk of false detection may significantly increase.

* `--window_size`, the length (*ms* in time) of each segment, the default is `30`. In our practice (300-second current traces), this is a suitable value, for the reason that smaller window sizes increase the computation complexity, while larger ones may lead to missing of events. We suggest that if tha length of the data is smaller, this parameter should be smaller so that sufficient slices are generated, which favor the follow-up statistical analyses.

*Please feel free to contact me if you have any questions.*
