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
  Note #1: It will provide a dataset named **fullDataset.csv** from FASTA sequences.
  ([fullDataset.csv](https://drive.google.com/file/d/1-DHKnHMcVDZATUYZ8BwQzdLUAWQJRxyg/))

  Note #2: **readXY.py** ( This file will fetch data from **hotSpot.fasta** and **coldSpot.fasta** files. )


- **Features Selection :**
  ```console
  user@machine:~$ python selectionFeatures.py
  ```
  Note #1: It will provide a dataset named **selectedDataset.csv** from **fullDataset.csv**.


- **Training Model :**
  ``` console
  user@machine:~$ python modelDump.py
  ```
  Note #1: It will provide a model named **iRecSpotModel.pkl**.

- **Test Dataset :**
  ``` console
  user@machine:~$ python evaluation.py
  ```
  Or, if we want to write results in a **.TXT** file then type:
  ``` console
  user@machine:~$ python evaluation.py > results.txt
  200
  # Input the stepSize (Here 200 is the dummy value.)
  # We can use any name instead of results.txt
  ```
  Note #1: In this portion, the programme will ask the user the **stepSize**.
  - **stepSize** : Segmentation length for each testing FASTA.
  - **Recommended stepSize >= 200.**
  
  Note #2: It will provide test result from **testFASTA.fasta** file.
  
  Note #3: readiRecX.py ( This file will fetch data from **testFASTA.fasta** file. )
  
  Note #3: If probability < 0.70 then the predicted type is assigned to 'None'. Otherwise, assigned to 'hotspot' or 'coldspot' based on the training model.

- **Run Machine Learning Classifiers :**
  ``` console
  user@machine:~$ python classifiers.py
  ```
  Or, if we want to write results in a **.TXT** file then type:
  ``` console
  user@machine:~$ python classifiers.py > classifiersResults.txt
  # We can use any name instead of classifiersResults.txt
  ```
  
  Note #1: It will provides:
  - Classification results from **selectedDataset.csv** file.
  - Generate a ROC Curve (**auROC.png**)
  - Generate an accuracy comparison via boxPlot. (**AccuracyBoxPlot.png**)

