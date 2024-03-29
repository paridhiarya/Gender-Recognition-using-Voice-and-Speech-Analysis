#import subprocess
import pandas as pd
import xgboost as xgb

#subprocess.call ("/usr/bin/Rscript --vanilla C:/Users/Paridhi/Dropbox/PC/Documents/Projects/Gender Recognition using Voice and Speech Analysis/Model/Feature_Extraction.R", shell=True)

features = ["meanfreq","sd","median","Q25","Q75","IQR","skew","kurt","sp.ent","sfm","mode","centroid","meanfun","minfun","maxfun","meandom","mindom","maxdom","dfrange","modindx"]

test_df = pd.read_csv(r'C:\Users\Paridhi\Dropbox\PC\Documents\Projects\Gender Recognition using Voice and Speech Analysis\Data\Brian-Acoustics.csv')

test_X = test_df[features]
xgtest = xgb.DMatrix(test_X)

model = xgb.Booster({'nthread':4})
model.load_model(r'C:\Users\Paridhi\Dropbox\PC\Documents\Projects\Gender Recognition using Voice and Speech Analysis\voice-gender.model')

pred_test_y = model.predict(xgtest)

if pred_test_y >= 0.5:
    print("Prediction: Male")

else:
    print("Prediction: Female")
