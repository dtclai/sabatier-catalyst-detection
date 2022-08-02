# sabatier-catalyst-detection

Step 1:
Run RF_methane(30runs).R with methanenew.csv to get the Top20 features.
Extract Top20 features and save as Top20.csv

Step 2:
For t-SNE visualisation, run tsne_noClusters.R on Top20.csv
For mBIC clustering with t-SNE visualisation, run mbic+tsne.R on Top20.csv
