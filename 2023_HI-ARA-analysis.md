---
title: "2023 HI ARA Analysis"
author: "Angel Avalos"
date: "2023-09-05"
output: 
  html_document: 
    keep_md: yes
---

# 2023 HI ARA Analysis

### Import packages and set working directory.

##### R

```r
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
library(reticulate)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
```

##### Python

```python
import pandas as pd
import numpy as np
import os
```

### ARA Data

```python
# Importing file
data = pd.read_csv("data/ara-reports.csv")
data = data.replace(r'^\s*$', np.nan, regex=True)
data = data.rename(columns={"Date":"total-date"})
# Extracting date
for i,v in enumerate(data["total-date"]):
  if type(v) == str:
    name = v[:-12]
    data.at[i,"total-date"] = name
  
data["date"] = pd.to_datetime(data["total-date"])
data = data.sort_values("date")

data.to_csv("output/ara-report-dates.csv",index=False)
```

### 2023 ARAs

```python
# in the future, convert to ethylene ppm based on standard curve
toptop=pd.DataFrame()
for i in os.listdir("data/"):
  if os.path.isdir("data/"+i):
    path = "data/"+i
    top = pd.DataFrame()
    for j in os.listdir(path):
      name = j.split("_",1)[1].replace("_rep1_MS.csv","").replace("_rep2_MS.csv","")
      data = pd.read_csv(path+"/"+j, header=3)
      data = data.iloc[:,1:24]
      # Note that these criteria are based on manual inspection of values, subject to change.
      data = data[data["RT"].between(2.61,2.68)]
      data.insert(loc=0,column="ID",value=name)
      top = pd.concat([top,data], axis=0)
    blank = top[top["ID"]=="uninoc_pos"]["Area"].mean()
    top = top[top["ID"]!="uninoc_pos"]
    top["Area"] = top["Area"]-blank
    top.reset_index(drop=True,inplace=True)
    top.insert(loc=24,column="date",value=i)
    toptop = pd.concat([top,toptop], axis=0)
toptop.reset_index(drop=True,inplace=True)
#data = pd.read_csv("data/20230902/20230902_BCW202070_pos_rep1_MS.csv", header=3)
#data = data.iloc[:,1:24]
#data = data[data["RT"].between(2.55,2.75)]
#name = "data/20230902/20230902_BCW202070_pos_rep1_MS.csv".split("_",1)[1].replace("_rep1_MS.csv","").replace("_rep2_MS.csv","")
#data.insert(loc=0,column="ID",value=name)
#blank = top[top["ID"]=="uninoc_pos"]["Area"].mean()
#top[top["ID"]!="uninoc_pos"]["Area"]-blank
```


```r
data = py$toptop

ggplot(data=data, aes(x=ID, y=Area)) +
  geom_bar(stat="summary",fun="mean", aes(fill=date)) + 
  geom_point() +
  ylab("Ethylene Peak Area") +
  theme(axis.text.x = element_text(angle=90, vjust=0.3))
```

![](2023_HI-ARA-analysis_files/figure-html/plot-1.png)<!-- -->
