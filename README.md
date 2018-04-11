# DTTrendMonitoring

CMS DT efficiency plots Vs *Instantaneous Luminosity, Integrated Luminosity* and other variables for trend monitoring.

More info about the code could be found in the pdf presentation stored in the `doc/` folder

[Wiki page](https://github.com/clacaputo/DTTrendMonitoring/wiki)

If in AFS you can set up the environment with 

```
source Setup.sh
```

- To create the csv file with integrated lumi per RUN number using brilcalc use ``` python  createJSONs.py ``` and follow the instruction therein written. The normal usage is;


```
python createJSONs.py -e fileName.root -n refName
```

For example

```
python createJSONs.py -e /eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/Express/DTTree_Run295463.root -n Run295463
```

- refName must have the format Run[number] in order to be used in the next step.



There are two type of run: 
- As a function of variables such integrated luminosity or time where you add new bins in order to have a trending plot.
  - The code is run_type1.C
- As a function of variables such instantaneous luminosity or PU where the variable range is fixed and you can increase the statistic using new files. 
  - The code is run_type2.C

To run the code for the first option, please follow these steps:

```
root 

.L run_type1.C+

```

Then to start from scratch with a new file

```
run_type1("refName","storingName","","file.root")

```

Where refName correspond to the reference name and it must be the same one used in the previous step for the integrate luminosity. 
storingName is the name of the root file where the results are stored.
The storage located in data/results . refName correspond also to the directory name that is saved in the web page with all the plots. 

As example:

```
run_type1("Run295463","295463","","/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/Express/DTTree_Run295463.root")

run_type1("Run2016B","test","","/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016BZMu23Sep2016-v1.root")
```

If there is already a stored file and you want to update the old plots with the new file please do:

```
run_type1("refName","storingName","storedName","file.root")

``` 

where storedName correspond to the file name in data/results the previous results are stored.

For the second option the procedure is identical.

It is also possible to run more then one file at the same time:

```
run_type1("refName","storingName","","file1.root","file2.root","file3.root")
```