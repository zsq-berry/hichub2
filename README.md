# hichub2

This is a python package to call collective chromatin interactions in one state and dynamics across multiple states from Hi-C data.

Call single state hichub:

![alt text](https://github.com/zsq-berry/hichub2/blob/4084af25db298dc2ee6adb0671eb11c3ffa95925/Picture1.png)

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

![alt text](https://github.com/zsq-berry/hichub2/blob/4d83ccb67a13424b5308eb256475401215453b44/Picture2.png)

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
