# DTTrendMonitoring

CMS DT efficiency plots Vs *Instantaneous Luminosity, Integrated Luminosity* and other variables for trend monitoring.

More info about the code could be found in the pdf presentation stored in the `doc/` folder

[Wiki page](https://github.com/clacaputo/DTTrendMonitoring/wiki)


An envoirment from cmssw is needed to have the right libraries. So you need to type cmsenv from a release > 8_0_X

Check-out the code from this repository. 

```
git clone https://github.com/cmsdtoffline/DTTrendMonitoring.git DTTrendMonitoring  
```
```
cd DTTrendMonitoring/
```

```
git checkout plotterClass
```


To set up brilcalc use:

```
source Setup.sh
```

- To create the csv file with integrated lumi per RUN number using brilcalc use ``` python  createJSONs.py ``` and follow the instruction therein written. The normal usage is;


```
python createJSONs.py -e fileName.root -n refName -y year
```

For example

```
python createJSONs.py -e /eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/Express/DTTree_Run295463.root -y 2017
```

- year is needed to compute the integrated luminosity from the begining of the first run of that year.

There are two type of run: 
- As a function of variables such integrated luminosity or time where you add new bins in order to have a trending plot.
  - The code is run_IncreasingRange.C
- As a function of variables such instantaneous luminosity or PU where the variable range is fixed and you can increase the statistic using new files. 
  - The code is run_FixedRange.C


To save the plots in a web page incrCont.webFolder must be modifed with the right path and the file index.php that is in the folder named plot has to be copied in the same path.

To run the code for the first option, please follow these steps:

```
root 

.L run_IncreasingRange.C+

```

Then to start from scratch with a new file

```
run_IncreasingRange("refName","storingName","","file.root")

```

Where refName correspond to the reference name and it must be the same one used in the previous step for the integrate luminosity. 
storingName is the name of the root file where the results are stored.
The storage located in data/results . refName correspond also to the directory name that is saved in the web page with all the plots. 

As example:

```
run_IncreasingRange("Run295463","295463","","/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/Express/DTTree_Run295463.root")

run_IncreasingRange("Run2016B","test","","/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016BZMu23Sep2016-v1.root")
```

If there is already a stored file and you want to update the old plots with the new file please do:

```
run_IncreasingRange("refName","storingName","storedName","file.root")

``` 

where storedName correspond to the file name in data/results the previous results are stored.

For the second option the procedure is identical.

It is also possible to run more then one file at the same time:

```
run_IncreasingRange("refName","storingName","","file1.root","file2.root","file3.root")
```