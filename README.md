# rmdup

`rmdup` is a python script which can be used to remove dumplicated and unreasonable alignments in the mapping file. `rmdup` is the upgraded version of [unique-sam](https://github.com/dlmeduLi/unique-sam). It use [pysam](http://pysam.readthedocs.org/en/latest/) as the underlying library to sort and parse the sam/bam file.

## Install

`rmdup` is denpended on [pysam](http://pysam.readthedocs.org/en/latest/) library, you should install `pysam` library first:
```shell
pip install pysam
```
`rmdup.py` is a python script, simply run the script in the command line:

```python
rmdup.py [options] in.bam -o out.bam 
```

## Usage

The input file of `rmdup.py` must be in `.bam` format. 

### parameters

```shell
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -o OUTPUTFILE, --output-file=OUTPUTFILE
                        write the result to output file
  -s, --sort            sort the input SAM file before further processing
  -k REGKEY, --reg-key=REGKEY
                        regexp to extract the key part of the read qname
  -n REGNUM, --reg-num=REGNUM
                        regexp to extract the number part of the read qname
  -u UPPERLIMIT, --upper-limit=UPPERLIMIT
                        the upper limit of the paired reads length
  -l LOWERLIMIT, --lower-limit=LOWERLIMIT
                        the lower limit of the paired reads length
  -t, --strict-mode     strict mode
```

## Unique Strategy

### Normal Mode (normal mode)

1. Keep alignment pairs with the highest score. (when more than 1 alignments with highest score appear, randomly keep one)
2. Reads are properly mapped (on the same chromsome & on the opposite strand)
3. Reads length should be in a proper range

### Strict Mode (strict mode `-t`)

1. Keep alignment pairs with the highest score. (when more than 1 alignments with highest score appear, **remove them all**)
2. Reads are properly mapped (on the same chromsome & on the opposite strand)
3. Reads length should be in a proper range

## Log File Specification

| Symbol | Description |
| :-: | :-: |
| !	| Error lines |
| <	| Low score alignments |
| =	| Pairs with more than one best score |
| ~	| Read pair mapped on the same strand |
| ? | Segment length too short |
| @ | Mapped positions overlapped |	

## Suggestions

If you have any problems or suggestions, please contact [dlmeduli@163.com](mailto:dlmeduli@163.com).

