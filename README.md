# AutoNanopore
**An automated adaptive and robust method to locate translocation events in solid-state nanopore current traces.**

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

* `--window_size`, the length (*ms* in time) of each segment, the default is `30`. In our practice, this is a suitable value, for the reason that smaller window sizes increase the computation complexity, while larger ones may lead to missing of events.

*Please feel free to contact me if you have any questions.*
