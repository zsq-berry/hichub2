# hichub2

This is a python package to call collective chromatin interactions in one state and dynamics across multiple states from Hi-C data.


Call single state hichub:
```
python3 hichub_single_state.py
Options:
  -h, --help            show this help message and exit
  -p <file>, --path=<file>
                        Input/output path
  -n <file>, --name=<file>
                        Hic file name
  -r <int>, --resolution=<int>
                        Resolution for the output hic file.

```

Call multiple states hichub:
```
python3 multiple.py
Options:
  -h, --help            show this help message and exit
  -p <file>, --path=<file>
                        Input/output path
  -n <file>, --names=<file>
                        HiC file names separated by comma
  -N <file>, --norm=<file>
                        HiC matrix normalization (NONE/KR)
  -c <int>, --nclusters=<int>
                        Number of clusters
  -r <int>, --resolution=<int>
                        Resolution fof hic matrix
```
