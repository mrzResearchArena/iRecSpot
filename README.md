## iRecSpot-EF: Effective Sequence Based Features for Recombination Hotspot Prediction

### Installation Process :
    Required Python Packages:

    Install: python (version >= 3.5)
    Install: sklearn (version >= 0.19.1)
    Install: numpy (version >= 1.13.1)
    Install: pandas (version >= 0.21.0)
    Install: matplotlib (version >= 2.1.0)
    Install: pickle

    pip install < package name >
    example: pip install sklearn

    or
    We can download from anaconda cloud.


### Code Description :
- **Features Extraction :**
  ```console
  user@machine:~$ python extractionFeatures.py
  ```
  Note #1: It will provide a dataset named **fullDataset.csv**.

  Note #2: **readXY.py** ( This file will fetch **hotSpot.fasta** and **coldSpot.fasta**. )


- **Features Selection :**
  ```console
  user@machine:~$ python selectionFeatures.py
  ```

  Note #1: It will provide a dataset named **selectedDataset.csv**.


- **Load Model :**
  ``` console
  user@machine:~$ python modelDump.py
  ```
  Note #1: It will provide a model named **iRecSpotModel.pkl**.

- **Test Dataset :**
  ``` console
  user@machine:~$ python evaluation.py
  ```
  Note #1: It will provide test result from **testFASTA.fasta** file.

- **Run Machine Learning Classifiers :**
  ``` console
  user@machine:~$ python classifiers.py
  ```
  Note #1: It will provides:
  - Classification results from **dataset.csv** file.
  - ROC Curve (**auROC.png**)
  - Accuracy comparison via boxPlot. (**AccuracyBoxPlot.png**)


