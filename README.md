# sabatier-catalyst-detection

Please cite the following work when using this code and data:

 ```
@article{ifandi2022,
   author={Elena Ifandi and Daphne Teck Ching Lai and Stavros Kalaitzidis and Muhammad Saifullah Abu Bakar and Tassos Grammatikopoulos and Chun-Kit Lai and Basilios Tsikouras},
   title={Detecting naturally occurring Sabatier catalyst(s) in a multivariate rock system using stochastic, machine-learning algorithms},
   journal={Under review},
   year ={2022}
}
```


Instructions to use code

Step 1:
Run RF_methane(30runs).R with methanenew.csv to get the Top20 features.
Extract Top20 features and save as Top20.csv

Step 2:
For t-SNE visualisation, run tsne_noClusters.R on Top20.csv
For mBIC clustering with t-SNE visualisation, run mbic+tsne.R on Top20.csv
