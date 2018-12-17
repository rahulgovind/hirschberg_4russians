# Hirschberg and the four russians

This project is an implementation of the Hirschberg algorithm for global alignment. We have also integrated the four
russians technique into our Hirschberg algorithm which lets us compute optimal edit distance alignments in
subquadratic time and linear memory.

# Setup

## Create the Makefile using cmake

```cmake .```

## Compile the source code using the generated makefile

```make```

# Usage

- For running standard Hirscberg, run ```./project -m hirschberg```
- For running Hirschberg with the Four Russians technique, run ```./project -m russians```

You will be prompted to enter the two sequences through standard input

If you wish to directly input files, you can run ```./project -m [method] -f /path/to/file```. The first line of the
file must be the first sequence and the second line must be the second sequence.