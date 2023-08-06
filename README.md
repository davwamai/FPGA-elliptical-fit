# FPGA-elliptical-fit
an FPGA approach for finding the best-fit for an elliptical point cloud. Utilizes Fitzgibbon's method and Cardano's method for eigensystem computation.

unfinished skeleton code. just here to give me a place to store my changes and allow me to update from work/home, but cool enough to show off its development.

# USAGE
This is a testing environment to ensure that the algorithm works as intended before porting to HLS. The HLS GUI is an unfinished nightmare and makes testing difficult, so this is a good way to emulate how data will be passed into the HLS module without having to define the hls_stream in a memory file. 

To run the code, clone the repository and navigate to the testing directory in a terminal, then simply run the make command:

```bash
make
```

Then run the executable like this:

```bash
./fitz
```
